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
