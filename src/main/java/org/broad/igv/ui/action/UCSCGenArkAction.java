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

import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.load.GenomeConfig;
import org.broad.igv.feature.genome.load.HubGenomeLoader;
import org.broad.igv.feature.genome.load.JsonGenomeLoader;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.genome.GenomeSelectionDialog;
import org.broad.igv.ui.genome.GenomeTableModel;
import org.broad.igv.ui.genome.GenomeTableRecord;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ParsingUtils;

import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * @author jrobinso
 * Date: 11/2/13
 * Time: 6:39 PM
 */
public class UCSCGenArkAction extends MenuAction {

    private static Logger log = LogManager.getLogger(UCSCGenArkAction.class);
    IGV igv;
    GenomeSelectionDialog genomeSelectionDialog;

    public UCSCGenArkAction(String label, int mnemonic, IGV igv) {
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
            List<String> headers = new ArrayList<>(List.of(new String[]{
                    "accession",
                    "assembly",
                    "scientific name",
                    "common name",
                    "taxonId"
            }));
            List<GenomeTableRecord> records = new ArrayList<>();

//            try (BufferedReader reader = ParsingUtils.openBufferedReader(PreferencesManager.getPreferences().getGenomeListURL())) {
//                String genomeRecord;
//                while ((genomeRecord = reader.readLine()) != null) {
//                    if (genomeRecord.startsWith("<") || genomeRecord.startsWith("#")) {
//                        continue;
//                    }
//                    genomeRecord = genomeRecord.trim();
//                    String[] fields = genomeRecord.split("\t");
//                    //# accession	assembly	scientific name	common name	taxonId	GenArk clade
//                    if (fields.length >= 3) {
//                        Map<String, String> attributes = Map.of(
//                                "common name", fields[0],
//                                "url", fields[1],
//                                "accession", fields[2],
//                                "source", "IGV"
//                        );
//                        records.add(new GenomeTableRecord(attributes));
//                    } else {
//                        log.error("Found invalid server genome list record: " + genomeRecord);
//                    }
//                }
//            } catch (IOException e) {
//                log.error(e);
//                MessageUtils.showMessage("Error loading genome list: " + e.getMessage());
//            }

            String[] fields = addRecords(PreferencesManager.getPreferences().getGenomeListURL(), records);
            addUniqueFields(fields, headers);

            fields = addRecords("https://hgdownload.soe.ucsc.edu/hubs/UCSC_GI.assemblyHubList.txt", records);
            addUniqueFields(fields, headers);

            GenomeTableModel model = new GenomeTableModel(headers, records);
            genomeSelectionDialog = new GenomeSelectionDialog(null, model);
            genomeSelectionDialog.setTitle("Genomes");
        }
        return genomeSelectionDialog;
    }

    private Set<String> ignoredFields = new HashSet<>(Arrays.asList(
            "url",
            "GenArk clade"
    ));

    private void addUniqueFields(String[] fields, List<String> headers) {
        for (String field : fields) {
            if (!ignoredFields.contains(field) && !headers.contains(field)) {
                headers.add(field);
            }
        }
    }

    private static String[] addRecords(String url, List<GenomeTableRecord> records) {
        String[] headers = null;
        try (BufferedReader reader = ParsingUtils.openBufferedReader(url)) {

            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) {
                    headers = Globals.tabPattern.split(line.substring(1).trim());
                } else {
                    System.out.println(line);
                    String[] values = Globals.tabPattern.split(line.trim());
                    Map<String, String> attributes = new HashMap<>();
                    for (int i = 0; i < headers.length; i++) {
                        attributes.put(headers[i], values[i]);
                    }
                    records.add(new GenomeTableRecord(attributes));
                }
            }
        } catch (IOException e) {
            log.error(e);
            MessageUtils.showMessage("Error loading genark list: " + e.getMessage());
        }
        return headers;

    }
}
