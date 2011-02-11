package org.broad.igv.ui.filefilters;

import org.apache.commons.io.FilenameUtils;

import javax.swing.filechooser.FileFilter;
import java.io.File;

/**
 * @author Fabien Campagne
 *         Date: Feb 11, 2011
 *         Time: 11:47:28 AM
 */
public class AlignmentFileFilter extends FileFilter {

    public boolean accept(File file) {
        String extension = FilenameUtils.getExtension(file.getName());
        if (extension != null) {

            return extension.equalsIgnoreCase("entries") ||
                    extension.equalsIgnoreCase("bam");
        } else {
            return false;
        }

    }

    public String getDescription() {
        return "All supported alignment files (*.bam, *.entries).";
    }

}

