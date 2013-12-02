/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
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
