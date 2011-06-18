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
                    file.getName().endsWith(".sorted.txt") || extension.equals("aligned");
        } else {
            return false;
        }

    }

    public String getDescription() {
        return "All supported alignment files";
    }

}

