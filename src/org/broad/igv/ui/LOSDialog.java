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
 * Created by JFormDesigner on Sat Sep 19 20:46:36 PDT 2015
 */

package org.broad.igv.ui;

import javax.swing.event.*;

import org.broad.igv.PreferenceManager;
import org.broad.igv.util.HttpUtils;

import java.awt.*;
import java.awt.event.*;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author James Robinson
 */
public class LOSDialog extends JDialog {

    public LOSDialog(Frame owner) {
        super(owner);
        initComponents();
        initContent();
        remindMeCB.setSelected(true);

        contentPane.setEditorKit(JEditorPane.createEditorKitForContentType("text/html"));
        contentPane.setEditable(false);

        contentPane.addHyperlinkListener(new HyperlinkListener() {
            public void hyperlinkUpdate(HyperlinkEvent e) {
                if(e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
                    try {
                        final URI mailtoURI = e.getURL().toURI();
                        Desktop.getDesktop().mail(mailtoURI);
                        setVisible(false);
                        PreferenceManager.getInstance().put(PreferenceManager.SHOW_LOS, false);
                        dispose();
                    } catch (Exception e1) {
                        e1.printStackTrace();
                    }
                }
            }
        });
    }

    private void initContent() {
        try {
            contentPane.setPage(new URL("http://www.broadinstitute.org/igv/projects/los.html"));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void okButtonActionPerformed(ActionEvent e) {

        if(!remindMeCB.isSelected()) {
            PreferenceManager.getInstance().put(PreferenceManager.SHOW_LOS, false);
        }

        this.setVisible(false);
        this.dispose();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        contentPane = new JEditorPane();
        remindMeCB = new JCheckBox();
        buttonBar = new JPanel();
        okButton = new JButton();

        //======== this ========
        Container contentPane2 = getContentPane();
        contentPane2.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BorderLayout());

                //======== scrollPane1 ========
                {

                    //---- contentPane ----
                    contentPane.setPreferredSize(new Dimension(580, 400));
                    contentPane.setContentType("text/html");
                    scrollPane1.setViewportView(contentPane);
                }
                contentPanel.add(scrollPane1, BorderLayout.CENTER);

                //---- remindMeCB ----
                remindMeCB.setText("Remind me later");
                contentPanel.add(remindMeCB, BorderLayout.SOUTH);
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
        contentPane2.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JScrollPane scrollPane1;
    private JEditorPane contentPane;
    private JCheckBox remindMeCB;
    private JPanel buttonBar;
    private JButton okButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
