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
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui;

import com.jidesoft.plaf.LookAndFeelFactory;
import org.broad.igv.ui.util.FileDialogUtils;

import java.io.File;

/**
 * @author eflakes
 */
public class ResourceFileBuilderDialog {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        LookAndFeelFactory.installDefaultLookAndFeel();

        File file = FileDialogUtils.chooseDirectory("Resource Folder Selection", null);

//        FileChooserDialog dialog = new FileChooserDialog(null, true);
//        dialog.setFileSelectionMode(FileChooserDialog.DIRECTORIES_ONLY);
//        dialog.setTitle("Resource Folder Selection");
//        dialog.setMultiSelectionEnabled(true);
//        dialog.setSelectedFile(null);
//        dialog.setDefaultCloseOperation(JDialog.EXIT_ON_CLOSE);
//        dialog.setVisible(true);

        if (file != null) {
            (new ResourceFileBuilder()).process(file);
        }
    }

}
