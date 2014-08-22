package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.ga4gh.Ga4ghAPIHelper;
import org.broad.igv.ui.IGV;

import java.awt.event.ActionEvent;
import java.io.IOException;

/**
 * Created by jrobinso on 8/16/14.
 */
public class BrowseGoogleAPIAction  extends MenuAction {

    private static Logger log = Logger.getLogger(BrowseGoogleAPIAction.class);


    IGV igv;

    public BrowseGoogleAPIAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }


    @Override
    public void actionPerformed(ActionEvent event) {

        try {
            Ga4ghAPIHelper.openLoadDialog(this.igv, IGV.getMainFrame());

        } catch (IOException e) {
            log.error("Error opening Encode browser", e);
        }
    }

}
