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

import com.google.gson.*;
import org.broad.igv.logging.*;
import org.broad.igv.feature.PSLRecord;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.PSLCodec;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.BlatTrack;
import org.broad.igv.track.SequenceTrack;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.*;

import java.io.IOException;
import java.net.URL;
import java.net.URLEncoder;
import java.util.*;

/**
 * Port of perl script blatPlot.pl   http://genomewiki.cse.ucsc.edu/index.php/Blat_Scripts
 *
 * @author jrobinso
 * Date: 11/21/12
 * Time: 8:28 AM
 */
public class BlatClient {

    private static Logger log = LogManager.getLogger(BlatClient.class);
    public static final int MINIMUM_BLAT_LENGTH = 20;


    public static List<PSLRecord> blat(String db, String userSeq) throws IOException {

        String serverType = PreferencesManager.getPreferences().get(Constants.BLAT_SERVER_TYPE, "");
        String urlpref = PreferencesManager.getPreferences().get(Constants.BLAT_URL);
        if (urlpref == null || urlpref.trim().length() == 0) {
            MessageUtils.showMessage("BLAT url is not configured");
            return Collections.EMPTY_LIST;
        }
        urlpref = urlpref.trim();

        if (serverType.equalsIgnoreCase("web_blat") || !urlpref.contains("$SEQUENCE")) {
            return LegacyBlatClient.blat(userSeq);

        } else {

            String dbEncoded = URLEncoder.encode(db, "UTF-8");
            urlpref = urlpref.replace("$SEQUENCE", userSeq).replace("$DB", dbEncoded);

            //Strip leading "file://" protocol, if any
            if (urlpref.startsWith("file://")) {
                urlpref = urlpref.substring("file://".length() + 1);
            }

            String jsonString = null;

            try {

                // If urlpref is not a URL, assume it is a command line program.  An example might be
                // blat.sh $SEQUENCE $DB

                if (URLUtils.isURL(urlpref)) {
                    jsonString = HttpUtils.getInstance().getContentsAsJSON(new URL(urlpref));
                } else {
                    jsonString = RuntimeUtils.exec(urlpref);
                }

                JsonObject obj = (new JsonParser()).parse(jsonString).getAsJsonObject();

                // If response includes "genome" property, verify against DB argument
                if (obj.has("genome")) {
                    String responseGenome = obj.get("genome").getAsString();
                    if (!(responseGenome.equalsIgnoreCase(db) || responseGenome.equalsIgnoreCase(dbEncoded))) {
                        throw new BlatException("Genome '" + db + "' not supported by BLAT server.");
                    }
                }

                // Collect PSL lines, stripping quotes from individual tokens in the process.
                JsonArray arr = obj.get("blat").getAsJsonArray();
                Iterator<JsonElement> iter = arr.iterator();
                List<String[]> results = new ArrayList<>();
                while (iter.hasNext()) {
                    JsonArray row = iter.next().getAsJsonArray();
                    String[] tokens = new String[row.size()];
                    for (int i = 0; i < row.size(); i++) {
                        String tmp = row.get(i).getAsString();
                        tokens[i] = StringUtils.stripQuotes(tmp);
                    }
                    results.add(tokens);
                }

                // Parse PSL lines into features.  The genome, if defined, is used to substitute chromosome aliases
                Genome genome = IGV.hasInstance() ? GenomeManager.getInstance().getCurrentGenome() : null;
                List<PSLRecord> features = new ArrayList<>(results.size());
                for (String[] tokens : results) {
                    features.add(PSLCodec.getPslRecord(tokens, genome));
                }

                return features;
            } catch (JsonSyntaxException e) {
                // This might be an html page, extract body text.
                String error = htmlToString(jsonString);
                throw new BlatException(error);
            } catch (InterruptedException e) {
                // Thrown from command line option
                throw new BlatException("Error executing blat command: " + e.getMessage(), e);
            }
        }
    }

    public static void doBlatQuery(final String userSeq, final String trackLabel) {

        try {
            Genome genome = GenomeManager.getInstance().getCurrentGenome();
            String db = genome.getBlatDB();
            if (db == null) {
                // If this is a known UCSC genome, or a custom blat server is configured, use the genome ID
                if (knownUCSCGenomes.contains(genome.getId()) ||
                        !PreferencesManager.getPreferences().get(Constants.BLAT_URL).equals(PreferencesManager.getDefault(Constants.BLAT_URL))) {
                    db = genome.getId();
                }
            }

            List<PSLRecord> features = blat(db, userSeq);
            if (features.isEmpty()) {
                MessageUtils.showMessage("No features found");
            } else {
                BlatTrack newTrack = new BlatTrack(db, userSeq, features, trackLabel); //species, userSeq, db, genome, trackLabel);
                IGV.getInstance().getTrackPanel(IGV.FEATURE_PANEL_NAME).addTrack(newTrack);
                IGV.getInstance().repaint();
                BlatQueryWindow win = new BlatQueryWindow(IGV.getInstance().getMainFrame(), userSeq, newTrack.getFeatures());
                win.setVisible(true);
            }
        } catch (Exception e1) {
            log.error("BLAT error.", e1);
            MessageUtils.showMessage(e1.getLocalizedMessage());
        }
    }


    public static void doBlatQueryFromRegion(final String chr, final int start, final int end, Strand strand) {

        if ((end - start) > 8000) {
            MessageUtils.showMessage("BLAT searches are limited to 8kb.  Please try a shorter sequence.");
            return;
        }

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        final byte[] seqBytes = genome.getSequence(chr, start, end);
        String userSeq = new String(seqBytes);

        if (strand == Strand.NEGATIVE) {
            userSeq = SequenceTrack.getReverseComplement(userSeq);
        }

        doBlatQuery(userSeq, "BLAT");
    }

    private static String htmlToString(String str) {

        // Remove script tags
        int idx1;
        int count = 0;
        while ((idx1 = str.indexOf("<script")) >= 0 && count < 10) {
            int idx2 = str.indexOf("</script>", idx1);
            if (idx1 > 0 && idx2 > idx1) {
                str = str.substring(0, idx1) + str.substring(idx2 + 9);
            }
            count++;
        }
        str = str.replaceAll("\\<.*?>", "").replace('\n', ' ');
        return str;

    }

    static Set knownUCSCGenomes = new HashSet<>(Arrays.asList("hg18", "hg19", "hg38", "ce10", "ce11", "galGal4",
            "galGal5", "galGal6", "panTro3", "panTro4", "panTro5", "panTro6", "felCat6", "bosTau7", "bosTau8",
            "bosTau9", "dm3", "dm6", "canFam3", "canFam5", "mm8", "mm9", "mm10", "mm39", "rn5", "rn6", "rn7",
            "rheMac3", "rheMac8", "rheMac10", "sacCer3", "susScr3", "xenTro9", "danRer10", "danRer11",
            "gorGor4", "gorGor6", "panPan2", "susScr11", "strPur2"));

}

