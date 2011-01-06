/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
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
package org.broad.igv.sam.reader;

import net.sf.samtools.SAMFileReader;
import org.apache.log4j.Logger;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.goby.GobyAlignmentQueryReader;
import org.broad.igv.util.IGVHttpUtils;
import org.broad.igv.util.ResourceLocator;
//import org.broad.igv.goby.GobyAlignmentQueryReader;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * @author jrobinso
 *         TODO consider renaming this class AlignmentQueryReaderFactory now that two alignment formats are supported.
 */
public class SamQueryReaderFactory {
    private static Logger log = Logger.getLogger(SamQueryReaderFactory.class);

    static {
        SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.SILENT);
    }

    public static AlignmentQueryReader getReader(String path, boolean requireIndex) throws IOException {
        return getReader(new ResourceLocator(path), requireIndex);
    }

    public static AlignmentQueryReader getReader(ResourceLocator locator) throws IOException {
        return getReader(locator, true);
    }

    public static AlignmentQueryReader getReader(ResourceLocator locator, boolean requireIndex) throws IOException {

        String path = locator.getPath().toLowerCase();


        AlignmentQueryReader reader = null;

        String samFile = locator.getPath();
        if (path.endsWith(".sam")) {
            reader = new SamQueryTextReader(samFile, requireIndex);

        } else if (path.endsWith("sorted.txt")
                || path.endsWith(".aligned")
                || path.endsWith(".aligned.txt")
                || path.endsWith("bedz")
                || path.endsWith("bed")
                || path.endsWith("psl")
                || path.endsWith("pslx")) {
            reader = new GeraldQueryReader(samFile, requireIndex);
        } else if (path.endsWith(".bam")) {
            if (locator.isLocal()) {
                reader = new BAMQueryReader(new File(samFile));
            } else if (IGVHttpUtils.isURL(locator.getPath().toLowerCase())) {
                try {
                    URL url = new URL(locator.getPath());
                    reader = new BAMHttpQueryReader(url, requireIndex);
                }
                catch (MalformedURLException e) {
                    log.error("", e);
                    throw new DataLoadException("Error loading BAM file: " + e.toString(), locator.getPath());
                }

            } else {
                reader = new BAMRemoteQueryReader(locator);
            }
        } else if (path.endsWith(".bam.list")) {
            reader = new BAMRemoteQueryReader(locator);
        } else if (locator.isLocal() && GobyAlignmentQueryReader.supportsFileType(locator.getPath())) {
            try {
                reader = new GobyAlignmentQueryReader(locator.getPath());
            } catch (IOException e) {
                throw new RuntimeException("Cannot load Goby alignment " + locator.getPath(), e);

            }
        } else {
            throw new RuntimeException("Cannot find reader for aligment file: " + locator.getPath());
        }

        return reader;
    }
}

