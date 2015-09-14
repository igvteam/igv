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

package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;

import java.io.IOException;
import java.io.InputStream;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;


public class GenomeZipDescriptor extends GenomeDescriptor {

    private static Logger log = Logger.getLogger(GenomeZipDescriptor.class);

    private Map<String, ZipEntry> zipEntries;
    private ZipFile genomeZipFile;

    public GenomeZipDescriptor(String name,
                               boolean chrNamesAltered,
                               String id,
                               String cytoBandFileName,
                               String geneFileName,
                               String chrAliasFileName,
                               String geneTrackName,
                               String sequenceLocation,
                               boolean hasCustomSequenceLocation,
                               ZipFile genomeZipFile,
                               Map<String, ZipEntry> zipEntries,
                               boolean chromosomesAreOrdered,
                               boolean fasta,
                               boolean fastaDirectory,
                               String fastaFileNameString) {
        super(name, chrNamesAltered, id, cytoBandFileName, geneFileName, chrAliasFileName, geneTrackName,
                sequenceLocation, hasCustomSequenceLocation, chromosomesAreOrdered, fasta, fastaDirectory, fastaFileNameString);
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


    public void close() {
        try {
            genomeZipFile.close();
        } catch (IOException e) {
            log.error("Error closing genomeZipFile", e);
        }
    }

}
