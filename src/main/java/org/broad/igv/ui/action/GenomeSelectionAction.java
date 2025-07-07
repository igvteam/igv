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

import org.broad.igv.feature.genome.HostedGenomes;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.genome.GenomeSelectionDialog;
import org.broad.igv.ui.genome.GenomeTableModel;
import org.broad.igv.ui.genome.GenomeDescriptor;

import java.awt.event.ActionEvent;
import java.util.Arrays;
import java.util.List;

/**
 * @author jrobinso
 * Date: 11/2/13
 * Time: 6:39 PM
 */
public class GenomeSelectionAction extends MenuAction {

    private static Logger log = LogManager.getLogger(GenomeSelectionAction.class);

        private static String[] defaultHeaders = {
            "common name",
            "scientific name",
            "assembly",
            "accession",
            "taxonId"
    };

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
            List<GenomeDescriptor> records = HostedGenomes.getRecords();
            List<String> headers = Arrays.asList(defaultHeaders);
            GenomeTableModel model = new GenomeTableModel(headers, records);
            genomeSelectionDialog = new GenomeSelectionDialog(IGV.getInstance().getMainFrame(), model);
            genomeSelectionDialog.setTitle("Genomes");
        }
        return genomeSelectionDialog;
    }

}
