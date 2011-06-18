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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.ui;

import com.jidesoft.plaf.LookAndFeelFactory;
import org.broad.igv.ui.util.FileChooserDialog;

import javax.swing.*;
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

        FileChooserDialog dialog = new FileChooserDialog(null, true);
        dialog.setFileSelectionMode(FileChooserDialog.DIRECTORIES_ONLY);
        dialog.setTitle("Resource Folder Selection");
        dialog.setMultiSelectionEnabled(true);
        dialog.setSelectedFile(null);
        dialog.setDefaultCloseOperation(JDialog.EXIT_ON_CLOSE);
        dialog.setVisible(true);

        if (!dialog.isCanceled()) {

            File file = dialog.getSelectedFile();
            (new ResourceFileBuilder()).process(file);

        }
    }

}
