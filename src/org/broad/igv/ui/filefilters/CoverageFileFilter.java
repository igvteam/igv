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

package org.broad.igv.ui.filefilters;

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

