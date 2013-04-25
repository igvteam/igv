package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.gs.GSFileBrowser;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.IGVMainFrame;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.event.ActionEvent;
import java.util.Arrays;

/**
 * @author jrobinso
 *         Date: 4/24/13
 *         Time: 4:26 PM
 */
public class LoadGenomeFromGSMenuAction extends MenuAction {

    static Logger log = Logger.getLogger(LoadGenomeFromGSMenuAction.class);
    IGV igv;

    public LoadGenomeFromGSMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
        setToolTipText("Load genome from GenomeSpace");
    }

    @Override
    public void actionPerformed(ActionEvent e) {


        try {
            GSFileBrowser dlg = new GSFileBrowser(IGV.getMainFrame());
            dlg.setVisible(true);

            String url = dlg.getFileURL();
            if (url != null) {
                igv.loadGenome(url, null);
            }
        } catch (Exception e1) {
            log.error("Error fetching directory listing on GenomeSpace server.", e1);
            MessageUtils.showMessage("Error fetching directory listing on GenomeSpace server: " + e1.getMessage());
        }

    }

}
