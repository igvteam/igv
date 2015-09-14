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

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.PSLCodec;
import org.broad.igv.track.FeatureCollectionSource;
import org.broad.igv.track.FeatureSource;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

/**
 * Port of perl script blatPlot.pl   http://genomewiki.cse.ucsc.edu/index.php/Blat_Scripts
 *
 * @author jrobinso
 *         Date: 11/21/12
 *         Time: 8:28 AM
 */
public class BlatClient {

    static int sleepTime = 15 * 1000;  //	#	milli seconds to wait between requests

    static String hgsid;  // cached, not sure what this is for but apparently its best to reuse it.
    static long lastQueryTime = 0;



    public static void main(String[] args) throws IOException {

        if (args.length != 6) {
            Usage();
            System.exit(255);
        }

        String org = args[0];
        String db = args[1];
        String searchType = args[2];
        String sortOrder = args[3];
        String outputType = args[4];
        String userSeq = args[5];

        blat(org, db, searchType, sortOrder, outputType, userSeq);

    }

    static void Usage() {
        System.out.println("usage: BlatBot <organism> <db> <searchType> <sortOrder>");
        System.out.println(" <outputType> <querySequence>");
        System.out.println("\tSpecify organism using the common name with first letter");
        System.out.println("capitalized.");
        System.out.println("\te.g. Human, Mouse, Rat etc.");
        System.out.println("\tDb is database or assembly name e.g hg17, mm5, rn3 etc.");
        System.out.println("\tsearchType can be BLATGuess, DNA, RNA, transDNA or transRNA");
        System.out.println("\tsortOrder can be query,score; query,start; chrom,score");
        System.out.println("\tchrom,start; score.");
        System.out.println("\toutputType can be pslNoHeader, psl or hyperlink.");
        System.out.println("\tblats will be run in groups of $batchCount sequences, all");
    }

    public static List<String> blat(String org, String db, String userSeq) throws IOException {
        String searchType = "DNA";
        String sortOrder = "query,score";
        String outputType = "psl";

        List<String> blatRecords = blat(org, db, searchType, sortOrder, outputType, userSeq);
        return blatRecords;


    }

    public static List<String> blat(String org, String db, String searchType, String sortOrder, String outputType, String userSeq) throws IOException {
        if (searchType.equals("BLATGuess")) {
            searchType = "Blat's Guess";
        } else if (searchType.equals("transDNA")) {
            searchType = "translated DNA";
        } else if (searchType.equals("transRNA")) {
            searchType = "translated RNA";
        } else if (searchType.equals("DNA") || (searchType.equals("RNA"))) {
        } else {
            System.out.println("ERROR: have not specified an acceptable search type - it should be BLATGuess, transDNA, transRNA, DNA or RNA.");
            Usage();
            System.exit(255);
        }
        if (outputType.equals("pslNoHeader")) {
            outputType = "psl no header";
        } else if (outputType.equals("psl") || outputType.equals("hyperlink")) {
        } else {
            System.out.println("ERROR: have not specified an acceptable output type - it should be pslNoHeader, psl or hyperlink.");
            Usage();
            System.exit(255);
        }

        //$response;
        String $url = PreferenceManager.getInstance().get(PreferenceManager.BLAT_URL);

        //if an hgsid was obtained from the output of the first batch
        //then use this.

        String urlString = ($url + "?org=" + org + "&db=" + db + "&type=" + searchType + "&sort=" + sortOrder +
                "&output=" + outputType + "&userSeq=" + userSeq); // + "&hgsid=" + hgsid);
        if (hgsid != null) {
            urlString += "&hgsid=" + hgsid;
        }

        long dt = System.currentTimeMillis() - lastQueryTime;
        if (dt < sleepTime) {
            try {
                Thread.sleep(dt);
            } catch (InterruptedException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

        lastQueryTime = System.currentTimeMillis();

        String result = HttpUtils.getInstance().getContentsAsString(new URL(urlString));

        return parseResult(result);
    }

    /**
     * Return the parsed results as an array of PSL records, where each record is simply an array of tokens.
     *
     * @param result
     * @return
     * @throws IOException
     */
    static List<String> parseResult(String result) throws IOException {

        ArrayList<String> records = new ArrayList<String>();

        BufferedReader br = new BufferedReader(new StringReader(result));
        String line;
        int headerLineCount = 0;
        boolean header = false;
        boolean hgsidFound = false;
        boolean pslSectionFound = false;
        while ((line = br.readLine()) != null) {

            if (line.contains("hgsid=") && !hgsidFound) {
                int startIDX = line.indexOf("hgsid=") + 6;
                String sub = line.substring(startIDX);
                int endIDX = sub.indexOf("\"");
                if (endIDX < 0) endIDX = sub.indexOf("&");
                if (endIDX > 0) {
                    hgsid = sub.substring(0, endIDX);
                    hgsidFound = true;
                }
            }
            if (line.trim().startsWith("<TT><PRE>")) {
                pslSectionFound = true;
                if (line.contains("psLayout") && line.contains("version")) {
                    header = true;
                    headerLineCount++;
                }
            } else if (line.trim().startsWith("</PRE></TT>")) {
                break;
            }

            if (pslSectionFound) {
                if (header && headerLineCount < 6) {
                    headerLineCount++;
                    continue;
                }

                String[] tokens = Globals.whitespacePattern.split(line);
                if (tokens.length != 21) {
                    System.err.println("Unexpected number of fields (" + tokens.length + ")");
                    System.err.println(line);
                } else {
                    records.add(line);
                }
            }
        }
        return records;
    }

    public static void doBlatQuery(final String chr, final int start, final int end) {

        if((end - start) > 8000) {
            MessageUtils.showMessage("BLAT searches are limited to 8kb.  Please try a shorter sequence.");
            return;
        }

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        final byte[] seqBytes = genome.getSequence(chr, start, end);
        String userSeq = new String(seqBytes);

        doBlatQuery(userSeq);
    }

    public static void doBlatQuery(final String userSeq) {

        LongRunningTask.submit(new NamedRunnable() {
            public String getName() {
                return "Blat sequence";
            }

            public void run() {
                try {

                    Genome genome = IGV.hasInstance() ? GenomeManager.getInstance().getCurrentGenome() : null;
                    PSLCodec codec = new PSLCodec(genome, true);

                    // TODO -- something better than this!
                    String db = genome.getId();
                    String species = genome.getSpecies();
                    if (species == null) {
                        MessageUtils.showMessage("Cannot determine species name for genome: " + genome.getDisplayName());
                        return;
                    }

                    List<String> tokensList = blat(species, db, userSeq);

                    // Convert tokens to features
                    List<PSLRecord> features = new ArrayList<PSLRecord>(tokensList.size());
                    for (String tokens : tokensList) {
                        PSLRecord f = codec.decode(tokens);
                        if (f != null) {
                            features.add(f);
                        }
                    }

                    if (features.isEmpty()) {
                        MessageUtils.showMessage("No features found");
                    } else {

                        FeatureSource<PSLRecord> source = new FeatureCollectionSource(features, genome);
                        FeatureTrack newTrack = new FeatureTrack("Blat", "Blat", source);
                        newTrack.setUseScore(true);
                        newTrack.setDisplayMode(Track.DisplayMode.SQUISHED);
                        IGV.getInstance().getTrackPanel(IGV.FEATURE_PANEL_NAME).addTrack(newTrack);


                        BlatQueryWindow win = new BlatQueryWindow(IGV.getMainFrame(), userSeq, features);
                        win.setVisible(true);

                    }
                } catch (IOException e1) {

                    MessageUtils.showErrorMessage("Error running blat", e1);
                }
            }
        });
    }
}

