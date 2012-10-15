/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * FileChooserDialog.java
 *
 * Created on November 1, 2007, 3:06 PM
 */
package org.broad.igv.ui.util;

import org.broad.igv.Globals;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;

/**
 * @author eflakes
 * @deprecated -- use FileDialogUtils
 */
public class FileChooserDialog extends javax.swing.JDialog {

    private File lastProgramaticallySelectedFile;

    private File previousFile;

    private File lastDirectory;

    private boolean isCanceled = true;

    private boolean autoRememberLastFile = false;

    private boolean autoRememberLastDirectory = false;

    private boolean isDisposableOnOk = false;

    private boolean isDisposableOnCancel = false;

    // Constants
    final public static int FILES_AND_DIRECTORIES = JFileChooser.FILES_AND_DIRECTORIES;

    final public static int FILES_ONLY = JFileChooser.FILES_ONLY;

    final public static int DIRECTORIES_ONLY = JFileChooser.DIRECTORIES_ONLY;

    /**
     * Creates new form FileChooserDialog
     */
    public FileChooserDialog(java.awt.Frame parent, boolean modal) {

        super(parent, modal);
        initComponents();
        setSize(725, 475);
        setLocationRelativeTo(parent);

        fileChooser.addPropertyChangeListener(new PropertyChangeListener() {

            public void propertyChange(final PropertyChangeEvent e) {

                if (e.getPropertyName().equals(JFileChooser.FILE_FILTER_CHANGED_PROPERTY)) {

                    File selectedFile = lastProgramaticallySelectedFile;
                    if (selectedFile == null) {
                        return;
                    }

                    Object filter = e.getNewValue();
                    if (filter instanceof SnapshotFileChooser.SnapshotFileFilter) {
                        SnapshotFileChooser.SnapshotFileFilter snapshotFilter =
                                (SnapshotFileChooser.SnapshotFileFilter) filter;

                        File newFile = changeFileExtension(selectedFile, snapshotFilter.getExtension());
                        fileChooser.setSelectedFile(newFile);
                    }
                }
            }
        });
    }

    @Override
    public void setVisible(boolean value) {

        // MAC has a bug in rescanCurrentDirectory() which is supposes to get 
        // fixed in Java version 6
        if (!Globals.IS_MAC) {
            fileChooser.rescanCurrentDirectory();
        }
        super.setVisible(value);
    }

    private File changeFileExtension(File file, String newExtension) {

        if (file == null ||
                file.isDirectory() ||
                newExtension == null ||
                newExtension.trim().equals("") ||
                newExtension.trim().equals(".")) {
            return file;
        }

        String filePath = file.getAbsolutePath();
        int index = filePath.lastIndexOf(".");

        if (index > 0) {
            filePath = filePath.substring(0, index);
        }

        return new File((filePath + newExtension));
    }

    public void setFileSelectionMode(int mode) {
        fileChooser.setFileSelectionMode(mode);
    }

    public File getSelectedFile() {
        return fileChooser.getSelectedFile();
    }

    public File[] getSelectedFiles() {
        return fileChooser.getSelectedFiles();
    }

    public boolean isCanceled() {
        return isCanceled;
    }

    public void setDisposeOnOk(boolean isDisposableOnOk) {
        this.isDisposableOnOk = isDisposableOnOk;
    }

    public void setDisposeOnCancel(boolean isDisposableOnCancel) {
        this.isDisposableOnCancel = isDisposableOnCancel;
    }

    public void addChoosableFileFilter(FileFilter filter) {
        fileChooser.addChoosableFileFilter(filter);
    }

    public FileFilter[] getChoosableFileFilters() {
        return fileChooser.getChoosableFileFilters();
    }

    public void removeChoosableFileFilter(FileFilter filter) {
        fileChooser.removeChoosableFileFilter(filter);
    }


    public void setFileFilter(FileFilter filter) {
        fileChooser.setFileFilter(filter);
    }

    public void setMultiSelectionEnabled(boolean isEnabled) {
        fileChooser.setMultiSelectionEnabled(true);
    }

    public void setSelectedFile(File file) {

        fileChooser.setSelectedFile(file);
        lastProgramaticallySelectedFile = file;

        if (autoRememberLastFile) {
            setPreviousFile(file);
        }
    }

    /**
     * Determines if a call to setCurrentFile will automatically remember
     * the file passed as the last file used. This will not remember
     * choices made via the UI - that must be set manually.
     *
     * @param value false if it should not remember last file (default = true)
     */
    public void setAutoRememberLastFile(boolean value) {
        autoRememberLastFile = value;
    }

    public void setPreviousFile(File lastFile) {
        this.previousFile = lastFile;
    }

    public File getPreviousFile() {
        return previousFile;
    }

    public void setCurrentDirectory(File file) {

        fileChooser.setCurrentDirectory(file);

        if (autoRememberLastDirectory) {
            setLastDirectory(file);
        }

        fileChooser.rescanCurrentDirectory();
    }

    public File getCurrentDirectory() {
        return fileChooser.getCurrentDirectory();
    }

    /**
     * Determines if a call to setCurrentDirectory will automatically remember
     * the file passed as the last directory used. This will not remember
     * choices made via the UI - that must be set manually.
     *
     * @param value false if it should not remember last directory (default = true)
     */
    public void setAutoRememberLastDirectory(boolean value) {
        autoRememberLastDirectory = value;
    }

    public void setLastDirectory(File lastDirectory) {
        this.lastDirectory = lastDirectory;
    }

    public File getLastDirectory() {

        if (lastDirectory == null) {
            lastDirectory = fileChooser.getCurrentDirectory();
        }

        return lastDirectory;
    }

    public void addFileChooserPropertyChangeListener(PropertyChangeListener listener) {
        fileChooser.addPropertyChangeListener(listener);
    }

    /**
     * This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc=" Generated Code ">//GEN-BEGIN:initComponents
    private void initComponents() {
        fileChooser = new javax.swing.JFileChooser();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        fileChooser.setApproveButtonText("Ok");
        fileChooser.setDialogTitle("Select a file or folder");
        fileChooser.setDialogType(javax.swing.JFileChooser.CUSTOM_DIALOG);
        fileChooser.setFileSelectionMode(javax.swing.JFileChooser.FILES_AND_DIRECTORIES);
        fileChooser.setMultiSelectionEnabled(true);
        fileChooser.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                fileChooserActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
                layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                        .add(layout.createSequentialGroup()
                                .addContainerGap()
                                .add(fileChooser, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
                layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
                                .addContainerGap()
                                .add(fileChooser, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 402, Short.MAX_VALUE))
        );
        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void fileChooserActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fileChooserActionPerformed

        final String action = evt.getActionCommand();

        java.awt.EventQueue.invokeLater(new Runnable() {

            public void run() {

                if (action.equals(javax.swing.JFileChooser.CANCEL_SELECTION)) {
                    isCanceled = true;
                    FileChooserDialog.this.setVisible(false);

                    if (isDisposableOnCancel) {
                        FileChooserDialog.this.dispose();
                    }
                } else if (action.equals(javax.swing.JFileChooser.APPROVE_SELECTION)) {
                    isCanceled = false;
                    FileChooserDialog.this.setVisible(false);

                    if (isDisposableOnOk) {
                        FileChooserDialog.this.dispose();
                    }
                }
            }
        });
    }//GEN-LAST:event_fileChooserActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        java.awt.EventQueue.invokeLater(new Runnable() {

            public void run() {
                new FileChooserDialog(new javax.swing.JFrame(), true).setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JFileChooser fileChooser;
    // End of variables declaration//GEN-END:variables
}
