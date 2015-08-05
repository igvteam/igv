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
 * Created by JFormDesigner on Sun Jun 05 16:03:30 EDT 2011
 */

package org.broad.igv.util;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Jim Robinson
 */
public class LoginDialog extends JDialog {

    boolean canceled = false;

    public LoginDialog(Frame owner) {
        this(owner, false, "", false);
    }

    public LoginDialog(Frame owner, boolean isGenomeSpace, String resource, boolean proxyChallenge) {
        super(owner);
        initComponents();
        getRootPane().setDefaultButton(okButton);
        if (isGenomeSpace) {
            promptLabel.setText("Please login to your GenomeSpace account");
            iconLabel.setVisible(true);
        } else {
            if (proxyChallenge) {
                promptLabel.setText("<html>Please enter username and password for your Proxy server to access<br>" + resource);
            } else {
                promptLabel.setText("<html>Please enter username and password to access<br>" + resource);
            }
            iconLabel.setVisible(false);
        }
    }

    public String getUsername() {
        return usernameField.getText();
    }

    public char[] getPassword() {
        return passwordField.getPassword();
    }

    public boolean isCanceled() {
        return canceled;
    }

    private void okButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        canceled = true;
        setVisible(false);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        iconLabel = new JLabel();
        promptLabel = new JLabel();
        passwordField = new JPasswordField();
        label2 = new JLabel();
        usernameField = new JTextField();
        label1 = new JLabel();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setModalityType(Dialog.ModalityType.APPLICATION_MODAL);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(null);

                //---- iconLabel ----
                iconLabel.setIcon(new ImageIcon(getClass().getResource("/images/genomespacelogo.png")));
                contentPanel.add(iconLabel);
                iconLabel.setBounds(20, 160, 320, iconLabel.getPreferredSize().height);

                //---- promptLabel ----
                promptLabel.setText("Enter username and password to access this resource");
                promptLabel.setFont(new Font("Arial", Font.PLAIN, 14));
                contentPanel.add(promptLabel);
                promptLabel.setBounds(5, 0, 395, 65);
                contentPanel.add(passwordField);
                passwordField.setBounds(125, 115, 220, 32);

                //---- label2 ----
                label2.setText("Password:");
                label2.setFont(new Font("Arial", Font.PLAIN, 14));
                contentPanel.add(label2);
                label2.setBounds(new Rectangle(new Point(25, 123), label2.getPreferredSize()));
                contentPanel.add(usernameField);
                usernameField.setBounds(125, 70, 220, 32);

                //---- label1 ----
                label1.setText("Username:");
                label1.setFont(new Font("Arial", Font.PLAIN, 14));
                contentPanel.add(label1);
                label1.setBounds(25, 72, 90, 28);

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

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JLabel iconLabel;
    private JLabel promptLabel;
    private JPasswordField passwordField;
    private JLabel label2;
    private JTextField usernameField;
    private JLabel label1;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
