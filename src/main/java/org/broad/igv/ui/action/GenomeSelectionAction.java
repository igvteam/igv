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

package org.broad.igv.ui.action;

import org.broad.igv.Globals;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.genome.GenomeSelectionDialog;
import org.broad.igv.ui.genome.GenomeTableModel;
import org.broad.igv.ui.genome.GenomeTableRecord;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;

import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

import static org.broad.igv.prefs.Constants.BACKUP_GENOMES_SERVER_URL;
import static org.broad.igv.prefs.Constants.GENOMES_SERVER_URL;

/**
 * @author jrobinso
 * Date: 11/2/13
 * Time: 6:39 PM
 */
public class GenomeSelectionAction extends MenuAction {

    private static Logger log = LogManager.getLogger(GenomeSelectionAction.class);

    private static String [] defaultHeaders = {
                    "common name",
                    "scientific name",
                    "assembly",
                    "accession",
                    "taxonId"
    };

    private static String[] legacyColumns = {
            "common name",
            "url",
            "id"
    };

    private static Set<String> ignoredFields = new HashSet<>(Arrays.asList(
            "url",
            "GenArk clade",
            "id",
            "_source"
    ));


    IGV igv;
    GenomeSelectionDialog genomeSelectionDialog;

    public GenomeSelectionAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }


    @Override
    public void actionPerformed(ActionEvent event) {
        final GenomeSelectionDialog dlg = getGenomeSelectionDialog();
        dlg.setVisible(true);
    }

    private GenomeSelectionDialog getGenomeSelectionDialog() {

        if (genomeSelectionDialog == null) {
            List<GenomeTableRecord> records = new ArrayList<>();
            final IGVPreferences preferences = PreferencesManager.getPreferences();
            final String genomesServerURL = preferences.get(GENOMES_SERVER_URL);
            final boolean isCustom =
                    !(genomesServerURL.equals(preferences.getDefault(GENOMES_SERVER_URL)) ||
                            genomesServerURL.equals(preferences.getDefault(BACKUP_GENOMES_SERVER_URL)));
            final String genarkURL = "https://hgdownload.soe.ucsc.edu/hubs/UCSC_GI.assemblyHubList.txt";

            List<String> headers;
            if (isCustom) {
                headers = new ArrayList<>();
            } else {
                headers = (List.of(defaultHeaders));
            }

            List<String> errors = new ArrayList<>();

            try {
                String [] fields = addRecords(genomesServerURL, records, "IGV");
                if(isCustom) {
                    addUniqueFields(fields, headers);
                }
            } catch (Exception e) {
                log.error("Error loading genome list from: " + genomesServerURL, e);
                if (isCustom) {
                    MessageUtils.showMessage("Error loading genome list from: " + genomesServerURL + "   (" + e.getMessage() + ")");
                    return null;
                } else {
                    errors.add("Error loading genome list from: " + genomesServerURL + "   (" + e.getMessage() + ")");
                    try {
                        addRecords(preferences.get(BACKUP_GENOMES_SERVER_URL), records, "IGV");
                    } catch (Exception ex) {
                        errors.add("Error loading genome list from: " + preferences.get(BACKUP_GENOMES_SERVER_URL) + "   (" + ex.getMessage() + ")");
                    }
                }
            }


            if (!isCustom) {
                try {
                    addRecords(genarkURL, records, "Genark");
                } catch (Exception e) {
                    log.error("Error connecting to UCSC Genark server URL: " + e.getMessage());
                    errors.add("Error connecting to UCSC Genark server: " + genarkURL + "  (" + e.getMessage() + ")");
                }
            }

            if (errors.size() > 0) {
                StringBuilder sb = new StringBuilder();
                for (String error : errors) {
                    sb.append(error).append("\n");
                }
                MessageUtils.showMessage(sb.toString());
            }

            GenomeTableModel model = new GenomeTableModel(headers, records);
            genomeSelectionDialog = new GenomeSelectionDialog(IGV.getInstance().getMainFrame(), model, !isCustom);
            genomeSelectionDialog.setTitle("Genomes");
        }
        return genomeSelectionDialog;
    }

    private void addUniqueFields(String[] fields, List<String> headers) {
        for (String field : fields) {
            if (!ignoredFields.contains(field) && !headers.contains(field)) {
                headers.add(field);
            }
        }
    }

    private static String[] addRecords(String url, List<GenomeTableRecord> records, String source) {
        String[] headers = null;
        try (BufferedReader reader = ParsingUtils.openBufferedReader(url)) {
            boolean firstRow = true;
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) {
                    headers = Globals.tabPattern.split(line.substring(1).trim());
                } else {
                    String[] values = Globals.tabPattern.split(line.trim());
                    if (firstRow) {
                        if (values.length == 3 && headers.length != 3) {
                            headers = legacyColumns;
                        }
                        firstRow = false;
                    }
                    Map<String, String> attributes = new HashMap<>();
                    for (int i = 0; i < headers.length; i++) {
                        attributes.put(headers[i], values[i]);
                    }
                    attributes.put("_source", source);
                    records.add(new GenomeTableRecord(attributes));
                }
            }
        } catch (IOException e) {
            log.error(e);
            MessageUtils.showMessage("Error loading genome list: " + url + "  " + e.getMessage());
        }
        return headers;

    }
}
