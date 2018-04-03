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
 * Created by JFormDesigner on Mon Aug 30 19:08:05 EDT 2010
 */

package org.broad.igv.util;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Jim Robinson
 */
public class UserPasswordDialog extends JDialog {

    boolean canceled = true;

    public UserPasswordDialog(Frame owner, String user, String host) {
        super(owner);
        setModal(true);
        initComponents();

        if (user != null) {
            userField.setText(user);
        }
        if(host != null) {
            hostLabel.setText(host);
        }
    }


    private void okButtonActionPerformed(ActionEvent e) {
        canceled = false;
        setVisible(false);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();
        contentPanel2 = new JPanel();
        label1 = new JLabel();
        hostLabel = new JLabel();
        label3 = new JLabel();
        label4 = new JLabel();
        userField = new JTextField();
        passwordField = new JPasswordField();
        panel1 = new JPanel();
        label2 = new JLabel();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

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

            //======== contentPanel2 ========
            {
                contentPanel2.setLayout(null);

                //---- label1 ----
                label1.setText("Enter user name and password for:");
                contentPanel2.add(label1);
                label1.setBounds(new Rectangle(new Point(15, 10), label1.getPreferredSize()));

                //---- hostLabel ----
                hostLabel.setText("ftp://ftp........................................");
                contentPanel2.add(hostLabel);
                hostLabel.setBounds(new Rectangle(new Point(30, 35), hostLabel.getPreferredSize()));

                //---- label3 ----
                label3.setText("User name:");
                contentPanel2.add(label3);
                label3.setBounds(new Rectangle(new Point(30, 125), label3.getPreferredSize()));

                //---- label4 ----
                label4.setText("Password");
                contentPanel2.add(label4);
                label4.setBounds(new Rectangle(new Point(30, 160), label4.getPreferredSize()));
                contentPanel2.add(userField);
                userField.setBounds(125, 120, 305, userField.getPreferredSize().height);
                contentPanel2.add(passwordField);
                passwordField.setBounds(125, 155, 305, passwordField.getPreferredSize().height);

                //======== panel1 ========
                {
                    panel1.setLayout(null);

                    { // compute preferred size
                        Dimension preferredSize = new Dimension();
                        for (int i = 0; i < panel1.getComponentCount(); i++) {
                            Rectangle bounds = panel1.getComponent(i).getBounds();
                            preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                            preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                        }
                        Insets insets = panel1.getInsets();
                        preferredSize.width += insets.right;
                        preferredSize.height += insets.bottom;
                        panel1.setMinimumSize(preferredSize);
                        panel1.setPreferredSize(preferredSize);
                    }
                }
                contentPanel2.add(panel1);
                panel1.setBounds(new Rectangle(new Point(125, 130), panel1.getPreferredSize()));

                //---- label2 ----
                label2.setText("<html>If this is a public server try \"anonymous\" for user name and your<br>\nemail address for password.");
                contentPanel2.add(label2);
                label2.setBounds(new Rectangle(new Point(15, 65), label2.getPreferredSize()));

                { // compute preferred size
                    Dimension preferredSize = new Dimension();
                    for (int i = 0; i < contentPanel2.getComponentCount(); i++) {
                        Rectangle bounds = contentPanel2.getComponent(i).getBounds();
                        preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                        preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                    }
                    Insets insets = contentPanel2.getInsets();
                    preferredSize.width += insets.right;
                    preferredSize.height += insets.bottom;
                    contentPanel2.setMinimumSize(preferredSize);
                    contentPanel2.setPreferredSize(preferredSize);
                }
            }
            dialogPane.add(contentPanel2, BorderLayout.NORTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    private JPanel contentPanel2;
    private JLabel label1;
    private JLabel hostLabel;
    private JLabel label3;
    private JLabel label4;
    private JTextField userField;
    private JPasswordField passwordField;
    private JPanel panel1;
    private JLabel label2;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public String getUser() {
        return userField.getText();
    }

    public String getPassword() {
        return new String(passwordField.getPassword());
    }

    public boolean isCanceled() {
        return canceled;
    }
}
