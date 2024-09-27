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
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.load.HubGenomeLoader;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.commandbar.GenomeListManager;
import org.broad.igv.ui.table.SearchableTableDialog;
import org.broad.igv.ui.table.SearchableTableModel;
import org.broad.igv.ui.table.SearchableTableRecord;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.encode.EncodeFileBrowser;
import org.broad.igv.util.encode.EncodeFileRecord;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 * Date: 11/2/13
 * Time: 6:39 PM
 */
public class UCSCGenArkAction extends MenuAction {

    private static Logger log = LogManager.getLogger(UCSCGenArkAction.class);
    IGV igv;

    public UCSCGenArkAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }


    @Override
    public void actionPerformed(ActionEvent event) {

        try (BufferedReader reader = ParsingUtils.openBufferedReader("https://hgdownload.soe.ucsc.edu/hubs/UCSC_GI.assemblyHubList.txt")) {

            String[] headers = null;
            List<SearchableTableRecord> records = new ArrayList<>();
            String headerLine = null;
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) {
                    headerLine = line;

                } else {
                    if (headers == null) {
                        headers = Globals.tabPattern.split(headerLine.substring(1).trim());
                    }
                    String[] values = Globals.tabPattern.split(line.trim());
                    Map<String, String> attributes = new HashMap<>();
                    for (int i = 0; i < headers.length; i++) {
                        attributes.put(headers[i], values[i]);
                    }
                    records.add(new SearchableTableRecord(attributes));
                }
            }
            SearchableTableModel model = new SearchableTableModel(headers, records);
            SearchableTableDialog dlg = new SearchableTableDialog(null, model);
            dlg.setTitle("UCSC GenArk");
            dlg.setVisible(true);

            SearchableTableRecord rec = dlg.getSelectedRecord();
            if (!dlg.isCanceled() && rec != null) {
                String accession = rec.getAttributeValue("accession");
                String url = HubGenomeLoader.convertToHubURL(accession);
                GenomeManager.getInstance().loadGenome(url);
            }
        } catch (IOException e) {
            log.error(e);
            MessageUtils.showMessage("Error loading GenArk hub: " + e.getMessage());
        }
    }
}
