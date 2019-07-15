
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

package org.broad.igv.ui;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.util.BrowserLauncher;

import java.awt.*;
import java.awt.event.*;
import java.io.IOException;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

import static org.broad.igv.ui.UIConstants.SERVER_BASE_URL;

/**
 * @author Jim Robinson
 */
public class VersionUpdateDialog extends JDialog {

    private static Logger log = Logger.getLogger(VersionUpdateDialog.class);


    boolean skipVersion;

    public VersionUpdateDialog(String versionString) {
        initComponents();
        messagePane.setEditorKit(JEditorPane.createEditorKitForContentType("text/html"));
        messagePane.setEditable(false);
        messagePane.setText("<strong>A later version of IGV is available  (" + versionString + ").<p/>Download from " +
                "<a href=" + Globals.downloadURL + ">" + Globals.downloadURL + "</a></strong>");
    }


    public boolean isSkipVersion() {
        return skipVersion;
    }

    private void okButtonActionPerformed(ActionEvent e) {
        setVisible(false);
        skipVersion = skipCB.isSelected();
    }

    private void messagePaneHyperlinkUpdate2(HyperlinkEvent e) {

        try {
            if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
                BrowserLauncher.openURL(e.getURL());
            }
        } catch (IOException ex) {
            log.error("Error opening browser", ex);
        }
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        skipCB = new JCheckBox();
        messagePane = new JEditorPane();
        buttonBar = new JPanel();
        okButton = new JButton();

        //======== this ========
        setModalityType(Dialog.ModalityType.APPLICATION_MODAL);
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BorderLayout(0, 10));

                //---- skipCB ----
                skipCB.setText("Skip this update");
                contentPanel.add(skipCB, BorderLayout.SOUTH);

                //---- messagePane ----
                messagePane.setFont(new Font("Lucida Grande", Font.BOLD, 14));
                messagePane.addHyperlinkListener(new HyperlinkListener() {
                    @Override
                    public void hyperlinkUpdate(HyperlinkEvent e) {
                        messagePaneHyperlinkUpdate2(e);
                    }
                });
                contentPanel.add(messagePane, BorderLayout.CENTER);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(540, 260);
        setLocationRelativeTo(null);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JCheckBox skipCB;
    private JEditorPane messagePane;
    private JPanel buttonBar;
    private JButton okButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
