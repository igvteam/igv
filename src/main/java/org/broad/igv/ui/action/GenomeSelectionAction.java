package org.broad.igv.ui.action;

import org.broad.igv.feature.genome.HostedGenomes;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.genome.GenomeSelectionDialog;
import org.broad.igv.ui.genome.GenomeTableModel;
import org.broad.igv.ui.genome.GenomeListItem;

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
            "accession"
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
            List<GenomeListItem> records = HostedGenomes.getRecords();
            String [] headers = defaultHeaders;
            GenomeTableModel model = new GenomeTableModel(headers, records);
            genomeSelectionDialog = new GenomeSelectionDialog(IGV.getInstance().getMainFrame(), model);
            genomeSelectionDialog.setTitle("Genomes");
        }
        return genomeSelectionDialog;
    }

}
