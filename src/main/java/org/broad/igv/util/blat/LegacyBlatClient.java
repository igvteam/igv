/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.util.blat;

import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.PSLCodec;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.HttpUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Port of perl script blatPlot.pl   http://genomewiki.cse.ucsc.edu/index.php/Blat_Scripts
 *
 * @author jrobinso
 * Date: 11/21/12
 * Time: 8:28 AM
 */
public class LegacyBlatClient {

    private static Logger log = LogManager.getLogger(LegacyBlatClient.class);

    static int sleepTime = 15 * 1000;  //	#	milli seconds to wait between requests

    static String hgsid;  // cached, not sure what this is for but apparently its best to reuse it.
    static long lastQueryTime = 0;


     static List<PSLRecord> blat(final String userSeq) throws IOException {

        String serverType = PreferencesManager.getPreferences().get(Constants.BLAT_SERVER_TYPE);
        List<String> tokensList;
        if (serverType.equalsIgnoreCase("web_blat")) {
            tokensList = webBlat(userSeq);

        } else {

            Genome genome = IGV.hasInstance() ? GenomeManager.getInstance().getCurrentGenome() : null;
            String db = genome.getId();
            String species = genome.getSpecies();
            if (species == null) {
                throw new RuntimeException("Cannot determine species name for genome: " + genome.getDisplayName());
            }
            tokensList = ucscBlat(species, db, userSeq);
        }

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        PSLCodec codec = new PSLCodec(genome, true);
        List<PSLRecord> features = new ArrayList<>(tokensList.size());
        for (String tokens : tokensList) {
            PSLRecord f = codec.decode(tokens);
            if (f != null) {
                features.add(f);
            }
        }
        return features;

    }


    private static List<String> ucscBlat(String org, String db, String userSeq) throws IOException {

        String searchType = "DNA";
        String sortOrder = "query,score";
        String outputType = "psl";

        String $url = PreferencesManager.getPreferences().get(Constants.BLAT_URL).trim();
        String serverType = PreferencesManager.getPreferences().get(Constants.BLAT_SERVER_TYPE);

        String result;
        if (serverType.equalsIgnoreCase("web_blat")) {
            String urlString = ($url + "?&wb_qtype=" + searchType + "&wb_sort=" + sortOrder +
                    "&wb_output=" + outputType + "&wb_seq=" + userSeq); // + "&hgsid=" + hgsid);
            //log.warn("BLAT: " + urlString);
            result = HttpUtils.getInstance().getContentsAsString(new URL(urlString));

        } else {
            String urlString = ($url + "?org=" + org + "&db=" + db + "&type=" + searchType + "&sort=" + sortOrder +
                    "&output=" + outputType); // + "&hgsid=" + hgsid);

            //if an hgsid was obtained from the output of the first batch then resuse.  I'm not sure what this is all about.
            if (hgsid != null) {
                urlString += "&hgsid=" + hgsid;
            }

            URL url = HttpUtils.createURL(urlString);
            long dt = System.currentTimeMillis() - lastQueryTime;
            if (dt < sleepTime) {
                try {
                    Thread.sleep(dt);
                } catch (InterruptedException e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
            lastQueryTime = System.currentTimeMillis();
            Map<String, String> params = new HashMap();
            params.put("userSeq", userSeq);
            result = HttpUtils.getInstance().doPost(url, params);

        }

        return parseResult(result);
    }

    public static List<String> webBlat(String userSeq) throws IOException {

        String searchType = "DNA";
        String sortOrder = "query,score";
        String outputType = "psl";

        String $url = PreferencesManager.getPreferences().get(Constants.BLAT_URL).trim();

        String result;
        String urlString = ($url + "?&wb_qtype=" + searchType + "&wb_sort=" + sortOrder +
                "&wb_output=" + outputType + "&wb_seq=" + userSeq); // + "&hgsid=" + hgsid);
        log.info("BLAT: " + urlString);
        result = HttpUtils.getInstance().getContentsAsString(new URL(urlString));

        List<String> records = parseResult(result);
        return fixWebBlat(records);
    }

    static List<String> fixWebBlat(List<String> records) {
        // Hack -- weblat appends the filename to sequenc names.  Strip it
        List<String> fixed = new ArrayList<>(records.size());
        for (String line : records) {
            if (line.startsWith("#")) {
                fixed.add(line);
                continue;
            }

            String fixedLine = "";
            String[] tokens = Globals.singleTabMultiSpacePattern.split(line);
            for (int i = 0; i < tokens.length; i++) {
                if (i > 0) {
                    fixedLine += "\t";
                }
                String t = tokens[i];
                if (i == 13) {
                    int idx = t.indexOf(":");
                    if (idx > 0) {
                        t = t.substring(idx + 1);
                    }
                }
                fixedLine += t;
            }
            fixed.add(fixedLine);
        }
        return fixed;
    }

    /**
     * Return the parsed results as an array of PSL records, where each record is simply an array of tokens.
     *
     * @param result
     * @return
     * @throws IOException
     */
    static List<String> parseResult(String result) throws IOException {

        List<String> records = new ArrayList<>();

        BufferedReader br = new BufferedReader(new StringReader(result));
        String l;
        boolean pslSectionFound = false;
        boolean pslHeaderFound = false;
        while ((l = br.readLine()) != null) {

            String line = l.trim();
            String lowerCase = line.toLowerCase();

            if (pslHeaderFound) {

                if (lowerCase.contains("</tt>")) {
                    break;
                }

                String[] tokens = Globals.whitespacePattern.split(line);
                if (tokens.length != 21) {
                    // PSL record section over
                    // Error?
                } else {
                    records.add(line);
                }
            }

            if (lowerCase.contains("<tt>") && lowerCase.contains("<pre>") && lowerCase.contains("pslayout")) {
                pslSectionFound = true;
                continue;
            }

            if (pslSectionFound) {
                if (lowerCase.startsWith("-----------------------------")) {
                    pslHeaderFound = true;
                }
            }
        }

        return records;
    }

}

