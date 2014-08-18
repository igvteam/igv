package org.broad.igv.ui.action;

import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ga4gh.GoogleAPIHelper;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.encode.EncodeFileBrowser;
import org.broad.igv.util.encode.EncodeFileRecord;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.util.*;

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
            GoogleAPIHelper.openLoadDialog(this.igv, IGV.getMainFrame() );

        } catch (IOException e) {
            log.error("Error opening Encode browser", e);
        }
    }

}
