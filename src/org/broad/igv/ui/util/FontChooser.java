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
 * Created by JFormDesigner on Thu May 19 21:44:40 EDT 2011
 */

package org.broad.igv.ui.util;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;

import org.broad.igv.ui.FontManager;
import org.jdesktop.layout.GroupLayout;
import org.jdesktop.layout.LayoutStyle;

/**
 * @author Jim Robinson
 */
public class FontChooser extends JDialog {

    private Font selectedFont;
    boolean canceled = false;

    public FontChooser(Dialog owner, Font font) {
        super(owner);
        initComponents();
        setTitle("Font Chooser");
        init(font);
    }

    private void init(Font font) {
        selectedFont = font;
        String[] fontFamilies = GraphicsEnvironment.getLocalGraphicsEnvironment().getAvailableFontFamilyNames();
        this.fontList.setListData(fontFamilies);
        String family = font.getFamily();
        fontList.setSelectedValue(font.getFamily(), true);
        sizeComboBox.setSelectedItem(String.valueOf(font.getSize()));
        exampleLabel.setFont(font);
    }

    private void fontListValueChanged(ListSelectionEvent e) {
        updateFont();
    }

    private void sizeComboBoxActionPerformed(ActionEvent e) {
        updateFont();
    }

    private void boldCBActionPerformed(ActionEvent e) {
        updateFont();
    }

    private void italicCBActionPerformed(ActionEvent e) {
        updateFont();
    }

    private void okButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
    }

    private void updateFont() {
        String fontName = (String) fontList.getSelectedValue();
        int size = Integer.parseInt((String) sizeComboBox.getSelectedItem());
        boolean isBold = boldCB.isSelected();
        boolean isItalic = italicCB.isSelected();
        int attrs = Font.PLAIN;
        if (isBold) attrs = Font.BOLD;
        if (isItalic) attrs |= Font.ITALIC;
        selectedFont = new Font(fontName, attrs, size);
        this.exampleLabel.setFont(selectedFont);
    }

    public Font getSelectedFont() {
        return selectedFont;
    }

    public boolean isCanceled() {
        return canceled;
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        sizeComboBox = new JComboBox();
        label1 = new JLabel();
        boldCB = new JCheckBox();
        italicCB = new JCheckBox();
        exampleLabel = new JTextPane();
        scrollPane1 = new JScrollPane();
        fontList = new JList();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {

                //---- sizeComboBox ----
                sizeComboBox.setModel(new DefaultComboBoxModel(new String[]{
                        "6",
                        "8",
                        "9",
                        "10",
                        "11",
                        "12",
                        "14",
                        "16",
                        "20",
                        "24"
                }));
                sizeComboBox.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        sizeComboBoxActionPerformed(e);
                    }
                });

                //---- label1 ----
                label1.setText("Size");
                label1.setLabelFor(sizeComboBox);

                //---- boldCB ----
                boldCB.setText("<html><b>Bold");
                boldCB.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        boldCBActionPerformed(e);
                    }
                });

                //---- italicCB ----
                italicCB.setText("<html><i>Italic");
                italicCB.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        italicCBActionPerformed(e);
                    }
                });

                //---- exampleLabel ----
                exampleLabel.setText("Lorem ipsum dolor sit amet, consectetur adipiscing elit.");

                //======== scrollPane1 ========
                {

                    //---- fontList ----
                    fontList.addListSelectionListener(new ListSelectionListener() {
                        public void valueChanged(ListSelectionEvent e) {
                            fontListValueChanged(e);
                        }
                    });
                    scrollPane1.setViewportView(fontList);
                }

                GroupLayout contentPanelLayout = new GroupLayout(contentPanel);
                contentPanel.setLayout(contentPanelLayout);
                contentPanelLayout.setHorizontalGroup(
                        contentPanelLayout.createParallelGroup()
                                .add(contentPanelLayout.createSequentialGroup()
                                .addContainerGap()
                                .add(scrollPane1, GroupLayout.PREFERRED_SIZE, 274, GroupLayout.PREFERRED_SIZE)
                                .add(22, 22, 22)
                                .add(contentPanelLayout.createParallelGroup(GroupLayout.TRAILING)
                                .add(contentPanelLayout.createSequentialGroup()
                                        .add(contentPanelLayout.createParallelGroup()
                                                .add(contentPanelLayout.createSequentialGroup()
                                                        .add(label1, GroupLayout.PREFERRED_SIZE, 52, GroupLayout.PREFERRED_SIZE)
                                                        .addPreferredGap(LayoutStyle.RELATED)
                                                        .add(sizeComboBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                                .add(boldCB, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                .add(italicCB, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                        .add(66, 66, 66))
                                .add(contentPanelLayout.createSequentialGroup()
                                .add(exampleLabel, GroupLayout.PREFERRED_SIZE, 198, GroupLayout.PREFERRED_SIZE)
                                .addContainerGap())))
                );
                contentPanelLayout.setVerticalGroup(
                        contentPanelLayout.createParallelGroup()
                                .add(contentPanelLayout.createSequentialGroup()
                                .addContainerGap()
                                .add(contentPanelLayout.createParallelGroup()
                                        .add(scrollPane1, GroupLayout.DEFAULT_SIZE, 264, Short.MAX_VALUE)
                                        .add(contentPanelLayout.createSequentialGroup()
                                        .add(contentPanelLayout.createParallelGroup(GroupLayout.BASELINE)
                                                .add(label1, GroupLayout.PREFERRED_SIZE, 22, GroupLayout.PREFERRED_SIZE)
                                                .add(sizeComboBox, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                        .addPreferredGap(LayoutStyle.RELATED)
                                        .add(boldCB, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .add(5, 5, 5)
                                        .add(italicCB, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(LayoutStyle.RELATED, 57, Short.MAX_VALUE)
                                        .add(exampleLabel, GroupLayout.PREFERRED_SIZE, 123, GroupLayout.PREFERRED_SIZE)))
                                .addContainerGap())
                );
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
        setSize(530, 365);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents


    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JComboBox sizeComboBox;
    private JLabel label1;
    private JCheckBox boldCB;
    private JCheckBox italicCB;
    private JTextPane exampleLabel;
    private JScrollPane scrollPane1;
    private JList fontList;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


}
