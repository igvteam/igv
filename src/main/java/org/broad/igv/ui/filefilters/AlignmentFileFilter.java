package org.broad.igv.ui.filefilters;

import org.apache.commons.io.FilenameUtils;

import javax.swing.filechooser.FileFilter;
import java.io.File;

/**
 * Filters filenames to only show alignment files. This is useful when users store a
 * number of Goby alignments in the same directory: only one file for each distinct alignment will be
 * available for selection.
 *
 * @author Fabien Campagne
 *         Date: Feb 11, 2011
 *         Time: 11:47:28 AM
 */
public class AlignmentFileFilter extends FileFilter implements java.io.FileFilter {

    public boolean accept(File file) {
        String extension = FilenameUtils.getExtension(file.getName());
        if (extension != null) {

            return extension.equalsIgnoreCase("entries") ||
                    extension.equalsIgnoreCase("bam") || extension.equals("sam") ||
                    extension.equals("aligned");
        } else {
            return false;
        }

    }

    public String getDescription() {
        return "All supported alignment files";
    }

}

