/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

package org.broad.igv.ui.util;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

/**
 * @author eflakes
 */
public class ProgressBar extends JPanel
        implements PropertyChangeListener {

    private boolean isReady = false;
    private JProgressBar progressBar;
    private boolean closeOnCompletion = false;
    private Window progressParentWindow;
    private ProgressMonitor monitor;

    public ProgressBar(int minimumProgress, int maximumProgress) {
        this(null, minimumProgress, maximumProgress, false, null);
    }

    public ProgressBar(Window progressParentWindow, int minimumProgress, int maximumProgress, boolean closeOnCompletion, ProgressMonitor monitor) {
        this.progressParentWindow = progressParentWindow;
        this.closeOnCompletion = closeOnCompletion;
        this.monitor = monitor;
        setLayout(new BorderLayout());
        progressBar = new JProgressBar(minimumProgress, maximumProgress);
        this.add(progressBar);
        setReady(true);
    }

    public void propertyChange(PropertyChangeEvent evt) {

        if (isReady()) {

            if (ProgressMonitor.PROGRESS_PROPERTY.equalsIgnoreCase(evt.getPropertyName())) {

                int value = (Integer) evt.getNewValue();
                final int progress;

                if (value > progressBar.getMaximum()) {
                    progress = 100;
                } else {
                    progress = value;
                }

                progressBar.setValue(progress);


                if (progress > 99 && closeOnCompletion) {
                    setReady(false); // Accept no more input
                    progressParentWindow.setVisible(false);
                }
            }
        } else {
            Object source = evt.getSource();
            if (source instanceof ProgressMonitor) {
                ((ProgressMonitor) source).setReady(false);
            }
        }
    }

    static public class ProgressDialog extends JDialog {

        public static boolean isAlreadyShowing = false;

        public ProgressDialog() {
            super();
        }

        public ProgressDialog(Frame frame) {
            super(frame);
        }
    }

    public static ProgressBar showProgressDialog(Frame parent, ProgressMonitor monitor, boolean closeOnCompletion) {
        return showProgressDialog(parent, "", monitor, closeOnCompletion);
    }

    public static ProgressBar showProgressDialog(Frame dialogsParent, String title, ProgressMonitor monitor, boolean closeOnCompletion) {

        ProgressDialog.isAlreadyShowing = true; // To prevent mutiple dialogs at same time
        ProgressDialog progressDialog = null;

        if (dialogsParent == null) {

            progressDialog = new ProgressDialog();
            progressDialog.setSize(500, 25);

            // Center on screen
            UIUtilities.centerWindow(progressDialog);
        } else {
        }
        progressDialog = new ProgressDialog(dialogsParent);
        progressDialog.setSize(500, 25);
        if (dialogsParent.isVisible()) {
            progressDialog.setLocationRelativeTo(dialogsParent);
        } else {

            // Center on screen
            UIUtilities.centerWindow(progressDialog);
        }

        final ProgressBar bar = new ProgressBar(progressDialog, 0, 100, closeOnCompletion, monitor);
        bar.setSize(500, 25);
        bar.setPreferredSize(bar.getSize());
        progressDialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        progressDialog.addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(WindowEvent e) {
                bar.setReady(false);
                ProgressDialog.isAlreadyShowing = false;
            }
        });
        progressDialog.setModal(false);
        progressDialog.setTitle(title);
        progressDialog.getContentPane().add(bar);
        progressDialog.pack();
        monitor.addPropertyChangeListener(bar);
        progressDialog.setVisible(true);

        return bar;
    }

    /**
     * Set the progress bars value.
     *
     * @param value The new value set on the progress bar.
     */
    public void setValue(int value) {
        progressBar.setValue(value);

        if (value > 99 && closeOnCompletion) {
            setReady(false); // Accept no more input
            progressParentWindow.setVisible(false);
        }
    }

    /**
     * Get the current progress amount from the progress bar.
     *
     * @return
     */
    public int getValue() {
        return progressBar.getValue();
    }

    /**
     * Reset the progress bar to 0.
     */
    public void reset() {
        progressBar.setValue(0);
    }

    /**
     * Sets whether or not this class will respond to progress requests.
     *
     * @param ready
     */
    public void setReady(boolean ready) {
        isReady = ready;
        if (monitor != null) {
            monitor.setReady(ready);
        }
    }

    /**
     * returns whether or not this class will respond to progress requests.
     *
     * @return
     */
    public boolean isReady() {
        return isReady;
    }


    public void close() {
        setReady(false);
        if (progressParentWindow != null) {
            setValue(100);
            progressParentWindow.setVisible(false);
            ProgressDialog.isAlreadyShowing = false;
        }
    }
}