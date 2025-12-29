/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.logging.*;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.LongRunningTask;

import java.awt.event.ActionEvent;

/**
 * @author jrobinso
 */
public class ReloadSessionMenuAction extends MenuAction {

    static Logger log = LogManager.getLogger(SaveSessionMenuAction.class);
    IGV igv;

    /**
     * @param label
     * @param mnemonic
     * @param igv
     */
    public ReloadSessionMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    /**
     * Method description
     *
     * @param e
     */
    @Override
    public void actionPerformed(ActionEvent e) {
        String currentSessionFilePath = igv.getSession().getPath();
        if (currentSessionFilePath != null) {
            LongRunningTask.submit(() -> this.igv.loadSession(currentSessionFilePath, null));
        }
    }
}
