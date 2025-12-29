/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.action;

import org.igv.logging.*;
import org.igv.track.Track;
import org.igv.track.TrackGroup;
import org.igv.ui.IGV;
import org.igv.ui.UIConstants;
import org.igv.ui.panel.DataPanelContainer;
import org.igv.ui.panel.TrackPanel;

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
