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
 * Created by JFormDesigner on Wed Mar 13 11:24:25 EDT 2013
 */

package org.broad.igv.ui.util;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

/**
 * @author Jacob Silterra
 */
public class CancellableProgressDialog extends JDialog {
    public CancellableProgressDialog(Frame owner) {
        super(owner);
        initComponents();
    }

    public CancellableProgressDialog(Dialog owner) {
        super(owner);
        initComponents();
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        permText = new JLabel();
        vSpacer1 = new JPanel(null);
        statusText = new JLabel();
        progressBar = new JProgressBar();
        buttonBar = new JPanel();
        hSpacer1 = new JPanel(null);
        button = new JButton();
        hSpacer2 = new JPanel(null);

        //======== this ========
        setAlwaysOnTop(true);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));
                contentPanel.add(permText);

                //---- vSpacer1 ----
                vSpacer1.setMaximumSize(new Dimension(12, 10));
                contentPanel.add(vSpacer1);

                //---- statusText ----
                statusText.setText("...");
                contentPanel.add(statusText);
                contentPanel.add(progressBar);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridLayout(1, 3));
                buttonBar.add(hSpacer1);

                //---- button ----
                button.setText("Cancel");
                button.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(button);
                buttonBar.add(hSpacer2);
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
    private JLabel permText;
    private JPanel vSpacer1;
    private JLabel statusText;
    private JProgressBar progressBar;
    private JPanel buttonBar;
    private JPanel hSpacer1;
    private JButton button;
    private JPanel hSpacer2;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    public JProgressBar getProgressBar() {
        return progressBar;
    }

    public void addButtonActionListener(ActionListener cancelActionListener) {
        button.addActionListener(cancelActionListener);
    }

    public void setStatus(final String status) {
        Runnable updater = new Runnable() {
            @Override
            public void run() {
                statusText.setText(status);
            }
        };
        UIUtilities.invokeOnEventThread(updater);
    }

    public void setPermText(final String perm){
        Runnable updater = new Runnable() {
            @Override
            public void run() {
                permText.setText(perm);
            }
        };
        UIUtilities.invokeOnEventThread(updater);
    }


    /**
     * Create a show a progress dialog with a single button (default text "Cancel")
     * @param dialogsParent
     * @param title
     * @param buttonActionListener The {@code ActionListener} to be called when the  button is pressed.
     * @param autoClose Whether to automatically close the dialog when it's finished
     * @param monitor Optional (may be null). Status text is updated based on monitor.updateStatus
     * @return
     */
    public static CancellableProgressDialog showCancellableProgressDialog(Frame dialogsParent, String title, final ActionListener buttonActionListener, final boolean autoClose, ProgressMonitor monitor){
        final CancellableProgressDialog progressDialog = new CancellableProgressDialog(dialogsParent);

        progressDialog.setTitle(title);
        progressDialog.addButtonActionListener(buttonActionListener);

        if(monitor != null && monitor instanceof IndefiniteProgressMonitor) progressDialog.getProgressBar().setIndeterminate(true);

        if(monitor != null){
            monitor.addPropertyChangeListener(new PropertyChangeListener() {
                @Override
                public void propertyChange(PropertyChangeEvent evt) {
                    if (evt.getPropertyName().equals(ProgressMonitor.STATUS_PROPERTY)) {
                        progressDialog.setStatus("" + evt.getNewValue());
                    } else if (evt.getPropertyName().equals(ProgressMonitor.PROGRESS_PROPERTY) && (Integer) evt.getNewValue() >= 100) {
                        progressDialog.button.setText("Done");
                        if(autoClose){
                            progressDialog.button.doClick(1);
                        }
                    }else if (evt.getPropertyName().equals(ProgressMonitor.PROGRESS_PROPERTY)) {
                        progressDialog.getProgressBar().setValue((Integer) evt.getNewValue());
                    }
                }
            });
        }

        progressDialog.setVisible(true);
        progressDialog.toFront();

        return progressDialog;
    }

}
