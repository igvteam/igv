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
 * Created by JFormDesigner on Fri Sep 17 11:22:26 EDT 2010
 */

package org.broad.igv.lists;

import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.util.MessageUtils;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Jim Robinson
 */
public class GeneListEditDialog extends JDialog {

    private GeneList geneList;
    private boolean canceled = true;
    private boolean bedOptionChanged = false;


    public GeneListEditDialog(Dialog owner, GeneList geneList) {
        super(owner);
        initComponents();
        this.geneList = geneList;

        String name = geneList.getName();
        if (name != null) {
            listNameField.setText(name);
        }

        StringBuffer buf = new StringBuffer();
        java.util.List<String> loci = geneList.getLoci();
        if (loci != null) {
            for (String gene : geneList.getLoci()) {
                buf.append(gene);
                buf.append("\n");
            }
            genesField.setText(buf.toString());
        }

        bedCB.setSelected(PreferenceManager.getInstance().getAsBoolean(PreferenceManager.GENE_LIST_BED_FORMAT));
    }


    private String[] parseGenes(String text) {
        return text.trim().split("\\s+");

    }


    private void okButtonActionPerformed(ActionEvent e) {
        String name = listNameField.getText();
        if (name == null || name.length() == 0) {
            MessageUtils.showMessage("Name is required");
            return;
        } else if (bedCB.isSelected()) {
            saveGeneList(name.trim(), parseBed(genesField.getText()));
        } else {
            String[] genes = parseGenes(genesField.getText());
            if (genes != null & genes.length == 0) {
                MessageUtils.showMessage("Lists must contain at least 1 locus");
                return;
            }
            saveGeneList(name.trim(), Arrays.asList(genes));
        }

        if(this.bedOptionChanged) {
            PreferenceManager.getInstance().put(PreferenceManager.GENE_LIST_BED_FORMAT, bedCB.isSelected());
        }
        
        setVisible(false);
    }

    private java.util.List<String> parseBed(String string) {

        java.util.List<String> loci = new ArrayList<String>();
        BufferedReader br = new BufferedReader(new StringReader(string));
        String nextLine;
        try {
            while ((nextLine = br.readLine()) != null) {
                String[] tokens = nextLine.split("\\s+");
                if (tokens.length > 2) {
                    loci.add(tokens[0] + ":" + (Integer.parseInt(tokens[1]) + 1) + "-" + tokens[2]);
                }
            }
        } catch (IOException e) {
            MessageUtils.showErrorMessage("Error parsing bed data", e);
        }
        return loci;
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void bedCBActionPerformed(ActionEvent e) {
        this.bedOptionChanged = true;
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
        scrollPane2 = new JScrollPane();
        descriptionField = new JTextArea();
        label3 = new JLabel();
        bedCB = new JCheckBox();
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
                contentPanel.setPreferredSize(new Dimension(435, 600));
                contentPanel.setLayout(null);

                //---- label1 ----
                label1.setText("Name: ");
                contentPanel.add(label1);
                label1.setBounds(5, 15, label1.getPreferredSize().width, 15);

                //---- listNameField ----
                listNameField.setPreferredSize(new Dimension(520, 28));
                contentPanel.add(listNameField);
                listNameField.setBounds(new Rectangle(new Point(65, 9), listNameField.getPreferredSize()));

                //---- label2 ----
                label2.setText("<html>Enter or paste genes or loci below &nbsp;&nbsp;<i>(e.g EGFR or chr1:1000-2000)");
                contentPanel.add(label2);
                label2.setBounds(5, 155, 425, 37);

                //======== scrollPane1 ========
                {
                    scrollPane1.setViewportView(genesField);
                }
                contentPanel.add(scrollPane1);
                scrollPane1.setBounds(6, 195, 579, 415);

                //======== scrollPane2 ========
                {
                    scrollPane2.setViewportView(descriptionField);
                }
                contentPanel.add(scrollPane2);
                scrollPane2.setBounds(5, 75, 579, 70);

                //---- label3 ----
                label3.setText("Description: ");
                contentPanel.add(label3);
                label3.setBounds(new Rectangle(new Point(5, 50), label3.getPreferredSize()));

                //---- bedCB ----
                bedCB.setText("BED format");
                bedCB.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        bedCBActionPerformed(e);
                    }
                });
                contentPanel.add(bedCB);
                bedCB.setBounds(new Rectangle(new Point(475, 160), bedCB.getPreferredSize()));

                { // compute preferred size
                    Dimension preferredSize = new Dimension();
                    for(int i = 0; i < contentPanel.getComponentCount(); i++) {
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
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

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
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    @Override
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
        setSize(650, 710);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
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
    private JScrollPane scrollPane2;
    private JTextArea descriptionField;
    private JLabel label3;
    private JCheckBox bedCB;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


    private void saveGeneList(String name, java.util.List<String> genes) {
        canceled = false;
        geneList.setName(name);
        geneList.setLoci(genes);
        geneList.setDescription(descriptionField.getText().trim());
        GeneListManager.getInstance().saveGeneList(geneList);
    }


    public boolean isCanceled() {
        return canceled;
    }
}
