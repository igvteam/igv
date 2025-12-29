/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.ui.action;

import org.igv.ui.IGV;

import java.awt.event.ActionEvent;

/**
 * @author jrobinso
 */
public class NewSessionMenuAction extends MenuAction {

    IGV igv;

    public NewSessionMenuAction(String label, int mnemonic, IGV igv) {
        super(label, null, mnemonic);
        this.igv = igv;
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        igv.newSession();
    }
}
