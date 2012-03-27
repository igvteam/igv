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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.feature.genome;

import org.apache.log4j.Logger;
import org.broad.igv.DirectoryManager;
import org.broad.igv.Globals;

import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.*;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

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
     * @param archiveOutputLocation
     * @param genomeFileName
     * @param genomeId                       Id of the genome.
     * @param genomeDisplayName              The genome name that is user-friendly.
     * @param sequenceLocation               The location of sequence data.
     * @param sequenceInputFile
     * @param refFlatFile                    RefFlat file.
     * @param cytobandFile                   Cytoband file.
     * @param sequenceOutputLocationOverride
     * @param monitor
     * @return The newly created genome archive file.
     */
    public File createGenomeArchive(File archiveOutputLocation,
                                    String genomeFileName,
                                    String genomeId,
                                    String genomeDisplayName,
                                    String sequenceLocation,
                                    File sequenceInputFile,
                                    File refFlatFile,
                                    File cytobandFile,
                                    File chrAliasFile,
                                    String sequenceOutputLocationOverride,
                                    ProgressMonitor monitor) throws IOException {

        if ((archiveOutputLocation == null) || (genomeFileName == null) || (genomeId == null) || (genomeDisplayName == null)) {

            log.error("Invalid input for genome creation: ");
            log.error("\tGenome Output Location=" + archiveOutputLocation);
            log.error("\tGenome filename=" + genomeFileName);
            log.error("\tGenome Id=" + genomeId);
            log.error("\tGenome Name" + genomeDisplayName);
            return null;
        }


        // Create a tmp directory for genome files
        File tmpdir = new File(DirectoryManager.getGenomeCacheDirectory(), genomeFileName + "_tmp");
        if (tmpdir.exists()) {
            tmpdir.delete();
        }
        tmpdir.mkdir();

        File propertyFile = null;

        File archive = null;
        FileWriter propertyFileWriter = null;
        try {
            boolean fastaDirectory = false;
            List<String> fastaFileNames = new ArrayList<String>();

            if (sequenceInputFile != null) {

                List<String> fastaIndexPathList = new ArrayList<String>();

                if (sequenceInputFile.isDirectory()) {
                    fastaDirectory = true;
                    List<File> files = getSequenceFiles(sequenceInputFile);
                    for (File file : files) {
                        if (file.getName().toLowerCase().endsWith(Globals.GZIP_FILE_EXTENSION)) {
                            String msg = "<html>Error.  One or more fasta files are gzipped: " + file.getName() +
                                    "<br>All fasta files must be gunzipped prior to importing.";
                            throw new GenomeException(msg);
                        }
                        String fastaPath = file.getAbsolutePath();
                        String fastaIndexPath = fastaPath + ".fai";
                        File indexFile = new File(fastaIndexPath);
                        if (!indexFile.exists()) {
                            FastaIndex.createIndexFile(fastaPath, fastaIndexPath);
                        }
                        fastaIndexPathList.add(fastaIndexPath);
                        fastaFileNames.add(file.getName());
                    }
                } else if (sequenceInputFile.getName().toLowerCase().endsWith(Globals.ZIP_EXTENSION)) {
                    String msg = "Error.  Zip archives are not supported.  Please select a fasta file.";
                    throw new GenomeException(msg);
                } else if (sequenceInputFile.getName().toLowerCase().endsWith(Globals.GZIP_FILE_EXTENSION)) {
                    String msg = "Error.  GZipped files are not supported.  Please select a non-gzipped fasta file.";
                    throw new GenomeException(msg);
                } else {
                    // Single fasta -- index
                    String fastaPath = sequenceInputFile.getAbsolutePath();
                    String fastaIndexPath = fastaPath + ".fai";
                    File indexFile = new File(fastaIndexPath);
                    if (!indexFile.exists()) {
                        FastaIndex.createIndexFile(fastaPath, fastaIndexPath);
                    }
                    fastaIndexPathList.add(fastaIndexPath);
                }

            }

            // Create Property File for genome archive
            if (sequenceOutputLocationOverride != null && sequenceOutputLocationOverride.length() > 0) {
                sequenceLocation = sequenceOutputLocationOverride;
            } else {
                sequenceLocation = FileUtils.getRelativePath(archiveOutputLocation.getAbsolutePath(), sequenceLocation);
            }


            propertyFile = createGenomePropertyFile(genomeId, genomeDisplayName, sequenceLocation, refFlatFile,
                    cytobandFile, chrAliasFile, fastaDirectory, fastaFileNames, tmpdir);
            archive = new File(archiveOutputLocation, genomeFileName);
            File[] inputFiles = {refFlatFile, cytobandFile, propertyFile, chrAliasFile};
            Utilities.createZipFile(archive, inputFiles);


        } finally {
            if (propertyFileWriter != null) {
                try {
                    propertyFileWriter.close();
                } catch (IOException ex) {
                    log.error("Failed to close genome archive: +" + archive.getAbsolutePath(), ex);
                }
            }

            if (propertyFile != null) propertyFile.delete();
            org.apache.commons.io.FileUtils.deleteDirectory(tmpdir);
        }
        return archive;
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
     * @param relativeSequenceLocation
     * @param refFlatFile
     * @param cytobandFile
     * @param fastaFileNames
     * @return
     */
    public File createGenomePropertyFile(String genomeId,
                                         String genomeDisplayName,
                                         String relativeSequenceLocation,
                                         File refFlatFile,
                                         File cytobandFile,
                                         File chrAliasFile,
                                         boolean fastaDirectory,
                                         List<String> fastaFileNames,
                                         File tmpdir) throws IOException {

        PrintWriter propertyFileWriter = null;
        try {

            File propertyFile = new File(tmpdir, "property.txt");
            propertyFile.createNewFile();

            // Add the new property file to the archive
            propertyFileWriter = new PrintWriter(new FileWriter(propertyFile));

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
            if (refFlatFile != null) {
                propertyFileWriter.println(Globals.GENOME_ARCHIVE_GENE_FILE_KEY + "=" + refFlatFile.getName());
            }
            if (chrAliasFile != null) {
                propertyFileWriter.println(Globals.GENOME_CHR_ALIAS_FILE_KEY + "=" + chrAliasFile.getName());
            }
            if (relativeSequenceLocation != null) {
                if (!HttpUtils.getInstance().isURL(relativeSequenceLocation)) {
                    relativeSequenceLocation = relativeSequenceLocation.replace('\\', '/');
                }
                propertyFileWriter.println(Globals.GENOME_ARCHIVE_SEQUENCE_FILE_LOCATION_KEY + "=" + relativeSequenceLocation);
            }
            return propertyFile;

        } finally {
            if (propertyFileWriter != null) {
                propertyFileWriter.close();

            }
        }

    }

}
