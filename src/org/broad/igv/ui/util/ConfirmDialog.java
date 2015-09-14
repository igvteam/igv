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
 * Created by JFormDesigner on Sun Aug 15 22:36:06 EDT 2010
 */

package org.broad.igv.ui.util;

import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * A confirm dialog with a "do not show again" checkbox.  For release 1.5rc2 this class is more or less hardocded
 * for the "check size" function of TrackLoader.
 *
 * @author Jim Robinson
 */
public class ConfirmDialog extends JDialog {

    String key;
    boolean okPressed = false;

    private ConfirmDialog(Frame owner, String message, String key) {
        super(owner);
        this.setModal(true);
        initComponents();
        label.setText("<html>" + message + "</html>");
        okButton.setText("Continue");
        this.key = key;
        getRootPane().setDefaultButton(cancelButton);
    }


    private void okButtonActionPerformed(ActionEvent e) {
        if (doNotShowCheckbox.isSelected()) {
            PreferenceManager.getInstance().put(key, "false");
        }
        okPressed = true;
        setVisible(false);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        if (doNotShowCheckbox.isSelected()) {
            PreferenceManager.getInstance().put(key, "false");
        }
        okPressed = false;
        setVisible(false);
    }


    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        label = new JLabel();
        buttonBar = new JPanel();
        doNotShowCheckbox = new JCheckBox();
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
                contentPanel.add(label);
                label.setBounds(10, 0, 565, 195);

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
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 0, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

                //---- doNotShowCheckbox ----
                doNotShowCheckbox.setText("Do not show this message again.");
                buttonBar.add(doNotShowCheckbox, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 5, 5), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(cancelButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 5, 5), 0, 0));

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 5, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(600, 300);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    public static void optionallyShowInfoDialog(String message, String key) {

        boolean show = PreferenceManager.getInstance().getAsBoolean(key);
        show &= ( !Globals.isHeadless() && !Globals.isSuppressMessages() );
        if (show) {
            ConfirmDialog dlg = new ConfirmDialog(IGV.getMainFrame(), message, key);
            dlg.okButton.setText("OK");
            dlg.cancelButton.setVisible(false);
            dlg.setVisible(true);
        }
    }

    public static boolean optionallyShowConfirmDialog(String message, String key, boolean defaultValue) {

        boolean show = PreferenceManager.getInstance().getAsBoolean(key);
        show &= !(Globals.isSuppressMessages() || Globals.isHeadless() || Globals.isTesting());
        if (show) {
            ConfirmDialog dlg = new ConfirmDialog(IGV.getMainFrame(), message, key);

            dlg.setVisible(true);
            return dlg.okPressed;

        }
        return defaultValue;
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JLabel label;
    private JPanel buttonBar;
    private JCheckBox doNotShowCheckbox;
    private JButton cancelButton;
    private JButton okButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


}
