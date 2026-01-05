package org.igv.ui.action;


import org.igv.logging.*;
import org.igv.ui.IGV;
import org.igv.util.LongRunningTask;

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
