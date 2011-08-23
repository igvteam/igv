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

package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;
import org.broad.tribble.util.SeekableStream;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Represents a sequence database composed of plain text files with no white space, one per chromosome, in a directory.
 * This is the original IGV "sequence" format, replaced in favor if indexed fasta files.
 *
 * @author jrobinso
 * @Date 8/8/11
 */

public class IGVSequence implements Sequence {

    private static Logger log = Logger.getLogger(IGVSequence.class);

    private String dirPath;
    private Map<String, String> chrFileNameCache = new HashMap();

    public IGVSequence(String dirPath) {
        if (!dirPath.endsWith("/")) {
            dirPath = dirPath + "/";
        }
        this.dirPath = dirPath;
    }

    public byte[] readSequence(String chr, int start, int end) {

        String fn = chr + ".txt";
        fn = getChrFileName(fn);
        String seqFile = dirPath + fn;

        SeekableStream is = null;
        try {


            is = IGVSeekableStreamFactory.getStreamFor(seqFile);

            byte[] bytes = new byte[end - start];
            is.seek(start);
            is.read(bytes);
            return bytes;

        } catch (Exception ex) {
            log.error("Error reading genome sequence from: " + seqFile, ex);
            return null;
        } finally {
            if (is != null) {
                try {
                    is.close();
                } catch (IOException ex) {
                    log.error("Error closing sequence file.", ex);
                }
            }
        }
    }


    private String getChrFileName(String fn) {
        String chrFN = chrFileNameCache.get(fn);
        if (chrFN == null) {
            chrFN = FileUtils.legalFileName(fn);
            chrFileNameCache.put(fn, chrFN);
        }
        return chrFN;
    }
}
