package org.igv.ui.filefilters;

import org.apache.commons.io.FilenameUtils;

import javax.swing.filechooser.FileFilter;
import java.io.File;

/**
 * Filters file names to only show coverage files (.tdf, .counts). This is useful when users store a
 * number of distinct type of files in the same directory: only coverage data becomes
 * available for selection.
 *
 * @author Fabien Campagne
 *         Date: Jun 11, 2011
 *         Time: 10:25:42 AM
 */
public class CoverageFileFilter extends FileFilter implements java.io.FileFilter {

    public boolean accept(File file) {
        String extension = FilenameUtils.getExtension(file.getName());
        if (extension != null) {

            return extension.equalsIgnoreCase("counts") ||
                    extension.equalsIgnoreCase("tdf") ;
        } else {
            return false;
        }

    }

    public String getDescription() {
        return "All supported coverage files";
    }

}

