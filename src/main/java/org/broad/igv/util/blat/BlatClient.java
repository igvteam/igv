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
import org.apache.log4j.Logger;
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
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * Port of perl script blatPlot.pl   http://genomewiki.cse.ucsc.edu/index.php/Blat_Scripts
 *
 * @author jrobinso
 * Date: 11/21/12
 * Time: 8:28 AM
 */
public class BlatClient {

    private static Logger log = Logger.getLogger(BlatClient.class);
    public static final int MINIMUM_BLAT_LENGTH = 20;


    public static List<PSLRecord> blat(String db, String userSeq) throws IOException {

        String serverType = PreferencesManager.getPreferences().get(Constants.BLAT_SERVER_TYPE, "");
        String urlpref = PreferencesManager.getPreferences().get(
                Constants.BLAT_URL,
                "https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq=$SEQUENCE&type=DNA&db=$DB&output=json")
                .trim();

        if (serverType.equalsIgnoreCase("web_blat") || !urlpref.contains("$SEQUENCE")) {
            return LegacyBlatClient.blat(userSeq);

        } else {

            urlpref = urlpref.replace("$SEQUENCE", userSeq).replace("$DB", db);

            //Strip leading "file://" protocol, if any
            if(urlpref.startsWith("file://")) {
                urlpref = urlpref.substring("file://".length() + 1);
            }

            String jsonString = null;

            try {

                // If urlpref is not a URL, assume it is a command line program.  An example might be
                // blat.sh $SEQUENCE $DB

                if(URLUtils.isURL(urlpref)) {
                    jsonString = HttpUtils.getInstance().getContentsAsJSON(new URL(urlpref));
                } else {
                    jsonString = RuntimeUtils.exec(urlpref);
                }


                JsonObject obj = (new JsonParser()).parse(jsonString).getAsJsonObject();
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

                Genome genome = IGV.hasInstance() ? GenomeManager.getInstance().getCurrentGenome() : null;
                List<PSLRecord> features = new ArrayList<>(results.size());
                for (String[] tokens : results) {
                    features.add(PSLCodec.getPslRecord(tokens, genome));
                }

                return features;
            } catch (JsonSyntaxException e) {
                // This might be an html page, extract body as message.
                String error = jsonString;
                int idx1 = jsonString.indexOf("<BODY");
                if (idx1 > 0) {
                    int idx2 = jsonString.indexOf(">", idx1);
                    int idx3 = jsonString.indexOf("</BODY", idx2);
                    if (idx3 > 0) {
                        error = jsonString.substring(idx2 + 1, idx3);
                    }
                }
                throw new RuntimeException(error);
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        }
    }

    public static void doBlatQuery(final String userSeq, final String trackLabel) {

        LongRunningTask.submit(new NamedRunnable() {

            public String getName() {
                return "Blat sequence";
            }

            public void run() {
                try {
                    Genome genome = IGV.hasInstance() ? GenomeManager.getInstance().getCurrentGenome() : null;
                    String db = genome == null ? "hg19" : genome.getUCSCId();
                    List<PSLRecord> features = blat(db, userSeq);


                    if (features.isEmpty()) {
                        MessageUtils.showMessage("No features found");
                    } else {
                        BlatTrack newTrack = new BlatTrack(userSeq, features, trackLabel); //species, userSeq, db, genome, trackLabel);
                        IGV.getInstance().getTrackPanel(IGV.FEATURE_PANEL_NAME).addTrack(newTrack);
                        IGV.getInstance().repaint();
                        BlatQueryWindow win = new BlatQueryWindow(IGV.getMainFrame(), userSeq, newTrack.getFeatures());
                        win.setVisible(true);
                    }
                } catch (Exception e1) {
                    MessageUtils.showErrorMessage("Error running blat", e1);
                }
            }
        });

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

        doBlatQuery(userSeq, "Blat");
    }


}

