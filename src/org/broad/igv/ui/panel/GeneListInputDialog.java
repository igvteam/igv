/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
 * Created by JFormDesigner on Fri Sep 17 11:22:26 EDT 2010
 */

package org.broad.igv.ui.panel;

import org.broad.igv.Globals;
import org.broad.igv.ui.util.MessageUtils;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.net.URLEncoder;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Jim Robinson
 */
public class GeneListInputDialog extends JDialog {

    private String[] genes;

    public GeneListInputDialog(Frame owner) {
        super(owner);
        initComponents();
    }

    public GeneListInputDialog(Dialog owner) {
        super(owner);
        initComponents();
    }

    private void okButtonActionPerformed(ActionEvent e) {
        if (listNameField.getText().length() == 0) {
            MessageUtils.showMessage("Name is required");
            genes = null;
        } else {
            parseGenes(genesField.getText());
            if(genes.length == 0) {
                MessageUtils.showMessage("Lists must contain at least 1 locus");
                genes = null;
            }
        }
        if(genes != null) {
            saveGeneList();
        }
        setVisible(false);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        genes = null;
        setVisible(false);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        label1 = new JLabel();
        listNameField = new JTextField();
        label2 = new JLabel();
        scrollPane1 = new JScrollPane();
        genesField = new JTextArea();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setModal(true);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(null);

                //---- label1 ----
                label1.setText("Gene list name: ");
                contentPanel.add(label1);
                label1.setBounds(new Rectangle(new Point(10, 21), label1.getPreferredSize()));
                contentPanel.add(listNameField);
                listNameField.setBounds(120, 15, 315, listNameField.getPreferredSize().height);

                //---- label2 ----
                label2.setText("Enter or paste gene list below");
                contentPanel.add(label2);
                label2.setBounds(10, 70, 240, label2.getPreferredSize().height);

                //======== scrollPane1 ========
                {
                    scrollPane1.setViewportView(genesField);
                }
                contentPanel.add(scrollPane1);
                scrollPane1.setBounds(10, 100, 425, 250);

                { // compute preferred size
                    Dimension preferredSize = new Dimension();
                    for (int i = 0; i < contentPanel.getComponentCount(); i++) {
                        Rectangle bounds = contentPanel.getComponent(i).getBounds();
                        preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                        preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                    }
                    Insets insets = contentPanel.getInsets();
                    preferredSize.width += insets.right;
                    preferredSize.height += insets.bottom;
                    contentPanel.setMinimumSize(preferredSize);
                    contentPanel.setPreferredSize(preferredSize);
                }
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 5), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    private void parseGenes(String text) {
        genes = text.trim().split("\\s+");
    }

    public String[] getGenes() {
        return genes;
    }

    public String getGeneListName() {
        return listNameField.getText();
    }
    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JLabel label1;
    private JTextField listNameField;
    private JLabel label2;
    private JScrollPane scrollPane1;
    private JTextArea genesField;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


    private void saveGeneList() {
        File file = null;
        try {
            final String listName = getGeneListName();
            file = new File(Globals.getGeneListDirectory(), getLegalFilename(listName));
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
            pw.println("#name=" + listName);
            for (String s : genes) {
                pw.println(s);
            }
            pw.close();
        } catch (IOException e) {
            if (file != null) {
                MessageUtils.showMessage("Error writing gene list file: " + file.getAbsolutePath() + " " + e.getMessage());
            }
            e.printStackTrace();
        }

    }


    private String getLegalFilename(String s) {
        try {
            return URLEncoder.encode(s, "UTF-8") + ".list.txt";
        } catch (UnsupportedEncodingException e) {
            return s;
        }
    }


}
