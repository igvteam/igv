/*
 * MenuAction.java
 *
 * Created on November 7, 2007, 2:07 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package org.broad.igv.ui.action;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * @author eflakes
 */
public class MenuAction extends AbstractAction {

    private String toolTipText;

    /**
     * Creates a new instance of MenuAction
     */
    public MenuAction(String name, Icon icon, int mnemonic) {
        super(name, icon);
        if(mnemonic >= 0) {
            putValue(MNEMONIC_KEY, mnemonic);
        }
    }

    public MenuAction(String name, Icon icon) {
        super(name, icon);
    }

    public MenuAction(String name) {
        super(name, null);
    }

    public void actionPerformed(ActionEvent event) {
        JOptionPane.showMessageDialog(null, "Functionality not implemented!");
    }

    public String getToolTipText() {
        return toolTipText;
    }

    public void setToolTipText(String toolTipText) {
        this.toolTipText = toolTipText;
    }

}
