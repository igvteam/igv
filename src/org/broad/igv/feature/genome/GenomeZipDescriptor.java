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

import java.io.IOException;
import java.io.InputStream;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;


public class GenomeZipDescriptor extends GenomeDescriptor {
    private Map<String, ZipEntry> zipEntries;
    private ZipFile genomeZipFile;

    public GenomeZipDescriptor(String name,
                               int version,
                               boolean chrNamesAltered,
                               String id,
                               String cytoBandFileName,
                               String geneFileName,
                               String chrAliasFileName,
                               String geneTrackName,
                               String sequenceLocation,
                               ZipFile genomeZipFile,
                               Map<String, ZipEntry> zipEntries,
                               boolean chromosomesAreOrdered) {
        super(name, version, chrNamesAltered, id, cytoBandFileName, geneFileName, chrAliasFileName, geneTrackName,
                sequenceLocation, chromosomesAreOrdered);
        this.zipEntries = zipEntries;
        this.genomeZipFile = genomeZipFile;

    }

    public InputStream getCytoBandStream()
            throws IOException {

        String fileName = cytoBandFileName;
        if (fileName == null) {
            return null;
        }

        boolean  isGZipped = fileName.toLowerCase().endsWith((".gz"));
        InputStream is = genomeZipFile.getInputStream(zipEntries.get(fileName));
        return isGZipped ? new GZIPInputStream(is) : is;
    }

    public InputStream getGeneStream()
            throws IOException {
        if (geneFileName == null) {
            return null;
        }
        InputStream is = genomeZipFile.getInputStream(zipEntries.get(geneFileName));
        return (geneFileName.endsWith(".gz") ? new GZIPInputStream(is) : is);
    }

    @Override
    public InputStream getChrAliasStream() throws IOException {

        String fileName = chrAliasFileName;
        if (fileName == null || !zipEntries.containsKey(fileName)) {
            return null;
        }
        return genomeZipFile.getInputStream(zipEntries.get(fileName));
    }

}
