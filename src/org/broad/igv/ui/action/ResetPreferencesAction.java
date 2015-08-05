/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui.action;

import org.broad.igv.ui.IGV;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.event.ActionEvent;

/**
 * @author jrobinso
 */
public class ResetPreferencesAction extends MenuAction {

    // TODO -- The batch referenceFrame is likely to be used by many actions. Move this
    // member to a base class ?
    IGV mainFrame;

    public ResetPreferencesAction(String label, IGV mainFrame) {
        super(label);
        this.mainFrame = mainFrame;
        setToolTipText(UIConstants.RESET_FACTORY_TOOLTIP);
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                execute();

            }
        });
    }

    public void execute() {

        int status =
                JOptionPane.showConfirmDialog(
                        mainFrame.getMainFrame(),
                        "All application preferences will be reset to" +
                                " their default settings." +
                                "\n Continue?",
                        null,
                        JOptionPane.OK_CANCEL_OPTION,
                        JOptionPane.PLAIN_MESSAGE,
                        null);

        if (status == JOptionPane.CANCEL_OPTION ||
                status == JOptionPane.CLOSED_OPTION) {
            return;
        }
        mainFrame.resetToFactorySettings();
    }
}
