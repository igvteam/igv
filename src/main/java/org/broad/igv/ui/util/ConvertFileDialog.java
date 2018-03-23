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
 * Created by JFormDesigner on Tue Nov 15 14:22:59 EST 2016
 */

package org.broad.igv.ui.util;

import org.broad.igv.ui.IGV;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;

/**
 * Before converting .shape, .map, .ct, .db, or .dp files (associated with individual RNAs)
 * to IGV-loadable formats, allow the user to select the applicable chromosome,
 * coordinate offset, and strand, since that information is not present in those files.
 *
 * @author sbusan
 */
public class ConvertFileDialog extends JDialog {

    ConvertOptions opts = new ConvertOptions();

    private ConvertFileDialog(Frame owner, String message, java.util.List<String> chromosomes) {
        super(owner);
        this.setModal(true);
        initComponents();
        label.setText("<html>" + message + "</html>");
        okButton.setText("Continue");

        // FIXME: limit input to empty string or integer
        /*NumberFormat format = NumberFormat.getInstance();
        NumberFormatter formatter = new NumberFormatter(format);
        formatter.setValueClass(Integer.class);
        formatter.setMinimum(0);
        formatter.setMaximum(Integer.MAX_VALUE);
        formatter.setAllowsInvalid(false);
        DefaultFormatterFactory factory = new DefaultFormatterFactory(formatter);
        startTextField.setFormatterFactory(factory);*/

        DefaultComboBoxModel boxModel = new DefaultComboBoxModel();
        for (String chrom : chromosomes){
            boxModel.addElement(chrom);
        }
        chromBox.setModel(boxModel);

        getRootPane().setDefaultButton(okButton);
    }

    public static ConvertOptions showConvertFileDialog(String message) {
        ConvertFileDialog dlg = new ConvertFileDialog(IGV.getMainFrame(),
                                                      message,
                                                      IGV.getInstance().getGenomeManager().getCurrentGenome().getAllChromosomeNames());
        dlg.setVisible(true);
        dlg.opts.chrom = dlg.chromBox.getSelectedItem().toString();
        dlg.opts.start = Integer.parseInt(dlg.startTextField.getText());
        if (dlg.reverseRadio.isSelected()) dlg.opts.strand = "-";
        return dlg.opts;
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        opts.doConvert = false;
        setVisible(false);
    }

    private void okButtonActionPerformed(ActionEvent e) {
        opts.doConvert = true;
        setVisible(false);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        label = new JLabel();
        forwardRadio = new JRadioButton();
        reverseRadio = new JRadioButton();
        label1 = new JLabel();
        label2 = new JLabel();
        label3 = new JLabel();
        startTextField = new JFormattedTextField();
        chromBox = new JComboBox();
        buttonBar = new JPanel();
        cancelButton = new JButton();
        okButton = new JButton();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setResizable(false);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(null);

                //---- label ----
                label.setText("text");
                label.setFont(label.getFont().deriveFont(label.getFont().getStyle() & ~Font.BOLD));
                contentPanel.add(label);
                label.setBounds(0, 0, 415, 140);

                //---- forwardRadio ----
                forwardRadio.setText("Forward");
                forwardRadio.setSelected(true);
                contentPanel.add(forwardRadio);
                forwardRadio.setBounds(65, 180, forwardRadio.getPreferredSize().width, 23);

                //---- reverseRadio ----
                reverseRadio.setText("Reverse");
                contentPanel.add(reverseRadio);
                reverseRadio.setBounds(145, 180, reverseRadio.getPreferredSize().width, 23);

                //---- label1 ----
                label1.setText("Strand:");
                label1.setHorizontalAlignment(SwingConstants.RIGHT);
                contentPanel.add(label1);
                label1.setBounds(5, 185, label1.getPreferredSize().width, 15);

                //---- label2 ----
                label2.setText("Chr:");
                label2.setHorizontalAlignment(SwingConstants.RIGHT);
                contentPanel.add(label2);
                label2.setBounds(10, 155, 42, 15);

                //---- label3 ----
                label3.setText("Start:");
                label3.setHorizontalAlignment(SwingConstants.RIGHT);
                contentPanel.add(label3);
                label3.setBounds(5, 215, 49, 15);

                //---- startTextField ----
                startTextField.setText("1");
                contentPanel.add(startTextField);
                startTextField.setBounds(65, 210, 165, startTextField.getPreferredSize().height);

                //---- chromBox ----
                chromBox.setMaximumRowCount(100);
                contentPanel.add(chromBox);
                chromBox.setBounds(65, 150, 165, chromBox.getPreferredSize().height);

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
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 0, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(e -> cancelButtonActionPerformed(e));
                buttonBar.add(cancelButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 5), 0, 0));

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(e -> okButtonActionPerformed(e));
                buttonBar.add(okButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 5, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(450, 355);
        setLocationRelativeTo(getOwner());

        //---- buttonGroup1 ----
        ButtonGroup buttonGroup1 = new ButtonGroup();
        buttonGroup1.add(forwardRadio);
        buttonGroup1.add(reverseRadio);
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JLabel label;
    private JRadioButton forwardRadio;
    private JRadioButton reverseRadio;
    private JLabel label1;
    private JLabel label2;
    private JLabel label3;
    private JFormattedTextField startTextField;
    private JComboBox chromBox;
    private JPanel buttonBar;
    private JButton cancelButton;
    private JButton okButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
