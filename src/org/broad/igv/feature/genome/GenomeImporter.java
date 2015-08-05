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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.HttpUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 * /**
 *
 * @author jrobinso
 */
public class GenomeImporter {
    public static final int MAX_CONTIGS = 1500000;

    static Logger log = Logger.getLogger(GenomeImporter.class);
    public static final Pattern SEQUENCE_NAME_SPLITTER = Pattern.compile("\\s+");


    /**
     * Create a zip containing all the information and data required to load a
     * genome. All file/directory validation is assume to have been done by validation
     * outside of this method.
     *
     * @param genomeFile
     * @param genomeId          Id of the genome.
     * @param genomeDisplayName The genome name that is user-friendly.
     * @param fastaFile         The location of a fasta file, or directory of fasta files
     * @param geneAnnotFile
     * @param cytobandFile
     * @return The newly created genome archive file.
     */
    public File createGenomeArchive(File genomeFile,
                                    String genomeId,
                                    String genomeDisplayName,
                                    String fastaFile,
                                    File geneAnnotFile,
                                    File cytobandFile,
                                    File chrAliasFile) throws IOException {

        if ((genomeFile == null) || (genomeId == null) || (genomeDisplayName == null)) {

            log.error("Invalid input for genome creation: ");
            log.error("\tGenome file=" + genomeFile.getAbsolutePath());
            log.error("\tGenome Id=" + genomeId);
            log.error("\tGenome Name" + genomeDisplayName);
            return null;
        }


        File propertyFile = null;
        FileWriter propertyFileWriter = null;
        try {
            boolean fastaDirectory = false;
            List<String> fastaFileNames = new ArrayList<String>();

            if (!FileUtils.resourceExists(fastaFile)) {
                String msg = "File not found: " + fastaFile;
                throw new GenomeException(msg);
            }
            if (fastaFile.toLowerCase().endsWith(Globals.ZIP_EXTENSION)) {
                String msg = "Error.  Zip archives are not supported.  Please select a fasta file.";
                throw new GenomeException(msg);
            }
            if (fastaFile.toLowerCase().endsWith(Globals.GZIP_FILE_EXTENSION)) {
                String msg = "Error.  GZipped files are not supported.  Please select a non-gzipped fasta file.";
                throw new GenomeException(msg);
            }

            List<String> fastaIndexPathList = new ArrayList<String>();
            String fastaIndexPath = fastaFile + ".fai";

            File sequenceInputFile = new File(fastaFile);
            if (sequenceInputFile.exists()) {
                // Local file
                if (sequenceInputFile.isDirectory()) {
                    fastaDirectory = true;
                    List<File> files = getSequenceFiles(sequenceInputFile);
                    for (File file : files) {
                        if (file.getName().toLowerCase().endsWith(Globals.GZIP_FILE_EXTENSION)) {
                            String msg = "<html>Error.  One or more fasta files are gzipped: " + file.getName() +
                                    "<br>All fasta files must be gunzipped prior to importing.";
                            throw new GenomeException(msg);
                        }

                        File indexFile = new File(sequenceInputFile, file.getName() + ".fai");
                        if (!indexFile.exists()) {
                            FastaUtils.createIndexFile(file.getAbsolutePath(), indexFile.getAbsolutePath());
                        }
                        fastaIndexPathList.add(fastaIndexPath);
                        fastaFileNames.add(file.getName());
                    }
                } else {
                    // Index if neccessary
                    File indexFile = new File(fastaIndexPath);
                    if (!indexFile.exists()) {
                        FastaUtils.createIndexFile(fastaFile, fastaIndexPath);
                    }
                    fastaIndexPathList.add(fastaIndexPath);
                }
            } else {
                if (!FileUtils.resourceExists(fastaIndexPath)) {
                    String msg = "<html>Index file " + fastaIndexPath + " Not found. " +
                            "<br>Remote fasta files must be indexed prior to importing.";
                    throw new GenomeException(msg);
                }
            }


            fastaFile = FileUtils.getRelativePath(genomeFile.getParent(), fastaFile);

            // Create "in memory" property file
            byte[] propertyBytes = createGenomePropertyFile(genomeId, genomeDisplayName, fastaFile, geneAnnotFile,
                    cytobandFile, chrAliasFile, fastaDirectory, fastaFileNames);
            File[] inputFiles = {geneAnnotFile, cytobandFile, chrAliasFile};

            // Create archive
            createGenomeArchive(genomeFile, inputFiles, propertyBytes);

        } finally {
            if (propertyFileWriter != null) {
                try {
                    propertyFileWriter.close();
                } catch (IOException ex) {
                    log.error("Failed to close genome archive: +" + genomeFile.getAbsolutePath(), ex);
                }
            }

            if (propertyFile != null) propertyFile.delete();
        }
        return genomeFile;
    }


    private List<File> getSequenceFiles(File sequenceDir) {
        ArrayList<File> files = new ArrayList();
        for (File f : sequenceDir.listFiles()) {
            if (f.getName().startsWith(".") || f.isDirectory() || f.getName().endsWith(".fai")) {
                continue;
            } else {
                files.add(f);
            }
        }
        return files;
    }


    /**
     * This method creates the property.txt file that is stored in each
     * .genome file. This is not the user-defined genome property file
     * created by storeUserDefinedGenomeListToFile(...)
     *
     * @param genomeId
     * @param genomeDisplayName
     * @param sequenceLocation Path to nucleotide sequence. Can be absolute or relative, also local or remote
     * @param geneAnnotFile
     * @param cytobandFile
     * @param fastaFileNames
     * @return
     */
    public byte[] createGenomePropertyFile(String genomeId,
                                           String genomeDisplayName,
                                           String sequenceLocation,
                                           File geneAnnotFile,
                                           File cytobandFile,
                                           File chrAliasFile,
                                           boolean fastaDirectory,
                                           List<String> fastaFileNames) throws IOException {

        PrintWriter propertyFileWriter = null;
        try {

            ByteArrayOutputStream propertyBytes = new ByteArrayOutputStream();

            // Add the new property file to the archive
            propertyFileWriter = new PrintWriter(new OutputStreamWriter(propertyBytes));

            propertyFileWriter.println("fasta=true"); // Fasta is the only format supported now

            propertyFileWriter.println("fastaDirectory=" + fastaDirectory);

            if (fastaDirectory) {
                propertyFileWriter.print("fastaFiles=");
                for (String fif : fastaFileNames) {
                    propertyFileWriter.print(fif + ",");
                }
                propertyFileWriter.println();
            }

            propertyFileWriter.println("ordered=" + !fastaDirectory);
            if (genomeId != null) {
                propertyFileWriter.println(Globals.GENOME_ARCHIVE_ID_KEY + "=" + genomeId);
            }
            if (genomeDisplayName != null) {
                propertyFileWriter.println(Globals.GENOME_ARCHIVE_NAME_KEY + "=" + genomeDisplayName);
            }
            if (cytobandFile != null) {
                propertyFileWriter.println(Globals.GENOME_ARCHIVE_CYTOBAND_FILE_KEY + "=" + cytobandFile.getName());
            }
            if (geneAnnotFile != null) {
                propertyFileWriter.println(Globals.GENOME_ARCHIVE_GENE_FILE_KEY + "=" + geneAnnotFile.getName());
            }
            if (chrAliasFile != null) {
                propertyFileWriter.println(Globals.GENOME_CHR_ALIAS_FILE_KEY + "=" + chrAliasFile.getName());
            }
            if (sequenceLocation != null) {
                if (!HttpUtils.isRemoteURL(sequenceLocation)) {
                    sequenceLocation = sequenceLocation.replace('\\', '/');
                }
                propertyFileWriter.println(Globals.GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY + "=" + sequenceLocation);
            }

            propertyFileWriter.flush();
            return propertyBytes.toByteArray();

        } finally {
            if (propertyFileWriter != null) {
                propertyFileWriter.close();

            }
        }

    }


    final static int ZIP_ENTRY_CHUNK_SIZE = 64000;

    static public void createGenomeArchive(File zipOutputFile, File[] inputFiles, byte[] propertyBytes)
            throws FileNotFoundException, IOException {

        if (zipOutputFile == null) {
            return;
        }

        if ((inputFiles == null) || (inputFiles.length == 0)) {
            return;
        }

        ZipOutputStream zipOutputStream = null;

        try {
            zipOutputStream = new ZipOutputStream(new FileOutputStream(zipOutputFile));

            ZipEntry propertiesEntry = new ZipEntry("property.txt");
            propertiesEntry.setSize(propertyBytes.length);
            zipOutputStream.putNextEntry(propertiesEntry);
            zipOutputStream.write(propertyBytes);


            for (File file : inputFiles) {

                if (file == null) {
                    continue;
                }

                long fileLength = file.length();

                ZipEntry zipEntry = new ZipEntry(file.getName());
                zipEntry.setSize(fileLength);
                zipOutputStream.putNextEntry(zipEntry);

                BufferedInputStream bufferedInputstream = null;
                try {
                    InputStream inputStream = new FileInputStream(file);
                    bufferedInputstream = new BufferedInputStream(inputStream);

                    int bytesRead = 0;
                    byte[] data = new byte[ZIP_ENTRY_CHUNK_SIZE];

                    while ((bytesRead = bufferedInputstream.read(data)) != -1) {
                        zipOutputStream.write(data, 0, bytesRead);
                    }
                } finally {
                    if (bufferedInputstream != null) {
                        bufferedInputstream.close();
                    }
                }
            }
        } finally {
            if (zipOutputStream != null) {
                zipOutputStream.flush();
                zipOutputStream.close();
            }
        }
    }

}
