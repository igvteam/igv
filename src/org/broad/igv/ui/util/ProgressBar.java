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

package org.broad.igv.ui.util;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

/**
 * Panel showing a progress bar, which can optionally close parent window when finished
 * @author eflakes
 */
public class ProgressBar extends JPanel
        implements PropertyChangeListener {

    private boolean isReady = false;
    private JProgressBar progressBar;
    private boolean closeOnCompletion = false;
    private Window progressParentWindow;
    private ProgressMonitor monitor;

    public ProgressBar(Window progressParentWindow, int minimumProgress, int maximumProgress, boolean closeOnCompletion, ProgressMonitor monitor) {
        this.progressParentWindow = progressParentWindow;
        this.closeOnCompletion = closeOnCompletion;
        this.monitor = monitor;
        setLayout(new BorderLayout());
        progressBar = new JProgressBar(minimumProgress, maximumProgress);
        progressBar.setIndeterminate(true); // Indeterminate by default
        this.add(progressBar);
        setReady(true);
    }

    public void setIndeterminate(boolean value) {
        progressBar.setIndeterminate(value);
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

        private ProgressBar progressBar = null;

        public ProgressDialog() {
            super();

        }

        public ProgressDialog(Frame frame) {
            super(frame);
        }

        public ProgressBar getProgressBar() {
            return progressBar;
        }

        public void setProgressBar(ProgressBar progressBar) {
            if(this.progressBar != null) throw new IllegalStateException("ProgressBar already set");
            this.progressBar = progressBar;
            getContentPane().add(progressBar);
        }
    }

    /**
     * Initialize a ProgressDialog but do not show it
     * @param dialogsParent
     * @param title
     * @param monitor
     * @param closeOnCompletion
     * @return
     */
    private static ProgressDialog createProgressDialog(Frame dialogsParent, String title, ProgressMonitor monitor, boolean closeOnCompletion) {

        ProgressDialog progressDialog = null;

        if (dialogsParent == null) {

            progressDialog = new ProgressDialog();
            progressDialog.setSize(500, 25);

            // Center on screen
            UIUtilities.centerWindow(progressDialog);
        } else {
            progressDialog = new ProgressDialog(dialogsParent);
            progressDialog.setSize(500, 25);
            if (dialogsParent.isVisible()) {
                progressDialog.setLocationRelativeTo(dialogsParent);
            } else {
                // Center on screen
                UIUtilities.centerWindow(progressDialog);
            }
        }

        final ProgressBar bar = new ProgressBar(progressDialog, 0, 100, closeOnCompletion, monitor);
        bar.setSize(500, 25);
        bar.setPreferredSize(bar.getSize());
        progressDialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        progressDialog.addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(WindowEvent e) {
                bar.setReady(false);
            }
        });
        progressDialog.setModal(false);
        progressDialog.setTitle(title);
        progressDialog.setProgressBar(bar);

        monitor.addPropertyChangeListener(bar);
        progressDialog.setModalExclusionType(Dialog.ModalExclusionType.APPLICATION_EXCLUDE);

        return progressDialog;
    }

    /**
     * Create and show a ProgressDialog
     * @param dialogsParent
     * @param title
     * @param monitor
     * @param closeOnCompletion
     * @return
     */
    public static ProgressDialog showProgressDialog(Frame dialogsParent, String title, ProgressMonitor monitor, boolean closeOnCompletion){
        ProgressDialog progressDialog = createProgressDialog(dialogsParent, title, monitor, closeOnCompletion);

        progressDialog.pack();
        progressDialog.setVisible(true);
        progressDialog.toFront();

        return progressDialog;
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
            if(progressParentWindow != null){
                progressParentWindow.setVisible(false);
            }
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

}