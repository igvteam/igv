/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.broad.igv.logging.*;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackGroup;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.panel.DataPanelContainer;
import org.broad.igv.ui.panel.TrackPanel;

import java.awt.event.ActionEvent;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 */
public class FitDataToWindowMenuAction extends MenuAction {

    static Logger log = LogManager.getLogger(FitDataToWindowMenuAction.class);
    IGV igv;

    public FitDataToWindowMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    /**
     * The action method. A swing worker is used, so "invoke later" and explicit
     * threads are not neccessary.
     *
     */
    public void actionPerformed(ActionEvent e) {
        igv.fitTracksToPanel();
    }



}
