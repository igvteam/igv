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

package org.broad.igv.feature.genome.load;

import org.apache.log4j.Logger;
import org.broad.igv.util.HttpUtils;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipInputStream;


public class GenomeDescriptor {

    public final static String GENOME_ARCHIVE_VERSION_KEY = "version";
    public final static String GENOME_ARCHIVE_PROPERTY_FILE_NAME = "property.txt";
    public final static String GENOME_ARCHIVE_ID_KEY = "id";
    public final static String GENOME_ARCHIVE_NAME_KEY = "name";
    public final static String GENOME_ORDERED_KEY = "ordered";
    public final static String GENOME_GENETRACK_NAME = "geneTrackName";
    public final static String GENOME_URL_KEY = "url";
    public final static String GENOME_ARCHIVE_CYTOBAND_FILE_KEY = "cytobandFile";
    public final static String GENOME_ARCHIVE_GENE_FILE_KEY = "geneFile";
    public final static String GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY = "sequenceLocation";
    public final static String COMPRESSED_SEQUENCE_PATH = "compressedSequencePath";
    public static final String GENOME_CHR_ALIAS_FILE_KEY = "chrAliasFile";
    public static final String SEQUENCE_MAP_FILE = "sequenceMap.txt";
    private static Logger log = Logger.getLogger(GenomeDescriptor.class);

    private Map<String, ZipEntry> zipEntries;
    private ZipFile genomeZipFile;

    private String id;
    private String name;
    protected String cytoBandFileName;
    protected String geneFileName;
    protected String chrAliasFileName;
    private String geneTrackName;
    private String url;
    private String sequencePath;
    private String compressedSequencePath;

    /**
     * Unzips a ".genome" archive and extracts metadata from property.txt.   Also has methods to create streams for
     * the embedded files.
     */
    public static GenomeDescriptor parseGenomeArchiveFile(File archiveFile)
            throws IOException {

        if (!archiveFile.exists()) {
            throw new FileNotFoundException("Genome file: " + archiveFile.getAbsolutePath() + " does not exist.");
        }

        GenomeDescriptor genomeDescriptor = null;
        Map<String, ZipEntry> zipEntries = new HashMap();
        ZipFile zipFile = new ZipFile(archiveFile);
        FileInputStream fileInputStream = null;
        try {
            fileInputStream = new FileInputStream(archiveFile);
            ZipInputStream zipInputStream = new ZipInputStream(fileInputStream);
            ZipEntry zipEntry = zipInputStream.getNextEntry();

            while (zipEntry != null) {
                String zipEntryName = zipEntry.getName();
                zipEntries.put(zipEntryName, zipEntry);
                zipEntry = zipInputStream.getNextEntry();
            }

            zipEntry = zipEntries.get(GENOME_ARCHIVE_PROPERTY_FILE_NAME);
            if (zipEntry != null) {
                InputStream inputStream = zipFile.getInputStream(zipEntry);
                Properties properties = new Properties();
                properties.load(inputStream);

                String sequencePath = properties.getProperty(GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY);
                String compressedSequencePath = properties.getProperty(COMPRESSED_SEQUENCE_PATH);

                if ((sequencePath != null) && !HttpUtils.isRemoteURL(sequencePath)) {
                    sequencePath = getFullPath(archiveFile, sequencePath);
                }

                if ((compressedSequencePath != null) && !HttpUtils.isRemoteURL(compressedSequencePath)) {
                    compressedSequencePath = getFullPath(archiveFile, sequencePath);
                }
                // The new descriptor
                genomeDescriptor = new GenomeDescriptor();
                genomeDescriptor.name = properties.getProperty(GENOME_ARCHIVE_NAME_KEY);
                genomeDescriptor.id = properties.getProperty(GENOME_ARCHIVE_ID_KEY);
                genomeDescriptor.cytoBandFileName = properties.getProperty(GENOME_ARCHIVE_CYTOBAND_FILE_KEY);
                genomeDescriptor.geneFileName = properties.getProperty(GENOME_ARCHIVE_GENE_FILE_KEY);
                genomeDescriptor.chrAliasFileName = properties.getProperty(GENOME_CHR_ALIAS_FILE_KEY);
                genomeDescriptor.geneTrackName = properties.getProperty(GENOME_GENETRACK_NAME, "Gene");
                genomeDescriptor.sequencePath = sequencePath;
                genomeDescriptor.compressedSequencePath = compressedSequencePath;
                genomeDescriptor.zipEntries = zipEntries;
                genomeDescriptor.genomeZipFile = zipFile;
                genomeDescriptor.url = properties.getProperty(GENOME_URL_KEY);
            }

        } finally {
            try {
                if (fileInputStream != null) {
                    fileInputStream.close();
                }
            } catch (IOException ex) {
                log.warn("Error closing imported genome zip stream!", ex);
            }
        }
        return genomeDescriptor;
    }

    static String getFullPath(File f, String sequencePath) throws IOException {
        File sequenceFolder = null;
        sequenceFolder = new File(sequencePath);
        boolean isAbsolutePath = sequenceFolder.isAbsolute() ||
                sequencePath.startsWith("/") || sequencePath.startsWith("\\");
        if (!isAbsolutePath) {
            sequenceFolder = new File(f.getParent(), sequencePath);
        }
        sequencePath = sequenceFolder.getCanonicalPath();
        sequencePath.replace('\\', '/');
        return sequencePath;
    }

    public String getName() {
        return name;
    }

    public String getId() {
        return id;
    }

    public String getGeneFileName() {
        return geneFileName;
    }

    public String getGeneTrackName() {
        return geneTrackName;
    }

    public String getSequencePath() {
        return compressedSequencePath == null ? sequencePath : compressedSequencePath;
    }

    public String getUrl() {
        return url;
    }

    public boolean hasCytobands() {
        return cytoBandFileName != null && cytoBandFileName.length() > 0;
    }

    public InputStream getCytoBandStream()
            throws IOException {

        String fileName = cytoBandFileName;
        if (fileName == null) {
            return null;
        }

        boolean isGZipped = fileName.toLowerCase().endsWith((".gz"));
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
