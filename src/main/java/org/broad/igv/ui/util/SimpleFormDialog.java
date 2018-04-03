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
 * Created by JFormDesigner on Fri Feb 10 08:19:01 EST 2012
 */

package org.broad.igv.ui.util;

import java.awt.*;
import java.awt.event.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import javax.swing.*;
import javax.swing.border.*;

/**
 * Constructs a simple multi-field text entry form from a list of labels.   Example use:
 *
 *    List labels = Arrays.asList("The first field", "F2", "The third field");
 *    SimpleFormDialog dlg = new SimpleFormDialog((Frame) null, "Prompt", labels);
 *    dlg.setVisible(true);
 *    if(!dlg.isCanceled) {
 *        String firstValue = dlg.getValue("The first field");
 *    }
 *
 * @author Jim Robinson
 */
public class SimpleFormDialog extends JDialog {


    boolean canceled = false;

    /**
     * Container to map labels to associated input boxes
     */
    public Map<String, JTextField> textComponents;

    public SimpleFormDialog(Frame owner, String prompt, java.util.List<String> labels) {
        super(owner);
        initComponents();
        initContentPanel(labels);
        promptLabel.setText(prompt);
    }

    public SimpleFormDialog(Dialog owner, String prompt, java.util.List<String> labels) {
        super(owner);
        initComponents();
        initContentPanel(labels);
        promptLabel.setText(prompt);
    }

    private void initContentPanel(java.util.List<String> labels) {

        textComponents = new HashMap<String, JTextField>(labels.size());
        contentPanel.setLayout(new GridLayout(labels.size(), 2, 5, 5));
        int width = 0;
        int height = 0;
        for (String l : labels) {
            JTextField field = new JTextField();
            textComponents.put(l, field);
            final JLabel jlabel = new JLabel(l);
            contentPanel.add(jlabel);
            contentPanel.add(field);

            Dimension d = jlabel.getPreferredSize();
            width = Math.max(width, (int) d.getWidth());

            double h2 = (int) field.getPreferredSize().getHeight();
            height += (int) Math.max(h2, d.getHeight()) + 5;
        }
        // Why is setting size explicitly neccessary?
        width = width*2 + 20;
        height += 20;
        contentPanel.setSize(width, height);
        contentPanel.setPreferredSize(new Dimension(width, height));
        pack();
    }

    public String getValue(String label) {
        return textComponents.get(label).getText();
    }

    public boolean isCanceled() {
        return canceled;
    }

    private void okButtonActionPerformed(ActionEvent e) {
        canceled = false;
        setVisible(false);
        dispose();
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
        dispose();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();
        panel1 = new JPanel();
        promptLabel = new JLabel();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(null);

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

            //======== panel1 ========
            {
                panel1.setBorder(new EmptyBorder(5, 0, 15, 5));
                panel1.setLayout(new BorderLayout());

                //---- promptLabel ----
                promptLabel.setFont(new Font("Lucida Grande", Font.BOLD, 13));
                panel1.add(promptLabel, BorderLayout.CENTER);
            }
            dialogPane.add(panel1, BorderLayout.NORTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    private JPanel panel1;
    private JLabel promptLabel;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


    // Example usage
    public static void main(String[] args) {

        java.util.List labels = Arrays.asList("The first field", "F2", "The third field");
        SimpleFormDialog dlg = new SimpleFormDialog((Frame) null, "Prompt", labels);
        dlg.setVisible(true);

    }
}
