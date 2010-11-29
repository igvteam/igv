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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.feature;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.GenomeManager;

import java.util.HashMap;
import java.util.HashSet;

/**
 * @author eflakes
 */
public class CommandLineGenomeBuilder {

    private final static String GENOME_OUTPUT_LOCATION = "outputLocation"; // Mandatory
    private final static String GENOME_DISPLAY_NAME = "genomeName";//Mandatory
    private final static String SEQUENCE_FILE = "sequenceFile"; // Mandatory
    private final static String GENOME_ID = "genomeId";
    private final static String CYTOBAND_FILE = "cytobandFile";
    private final static String GENE_FILE = "geneFile";
    private final static String SEQUENCE_OUTPUT_LOCATION_OVERRIDE = "sequenceLocationOverride";
    private final static String GENOME_ARCHIVE_FILENAME_PREFIX_OVERRIDE = "archiveFilenamePrefix";

    static final HashSet<String> validCommandLineOptions = new HashSet();

    static {
        validCommandLineOptions.add(GENOME_OUTPUT_LOCATION);
        validCommandLineOptions.add(GENOME_DISPLAY_NAME);
        validCommandLineOptions.add(SEQUENCE_FILE);
        validCommandLineOptions.add(CYTOBAND_FILE);
        validCommandLineOptions.add(GENE_FILE);
        validCommandLineOptions.add(SEQUENCE_OUTPUT_LOCATION_OVERRIDE);
        validCommandLineOptions.add(GENOME_ARCHIVE_FILENAME_PREFIX_OVERRIDE);
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {

        /*
       args = new String[] {
           (GENOME_OUTPUT_LOCATION+"=C:/transfer"),
           (GENOME_DISPLAY_NAME+"=blahblah2"),
           (SEQUENCE_FILE+"=C:/Alex Test Data/input/Cgla.fsa"),
           (GENOME_ARCHIVE_FILENAME_PREFIX_OVERRIDE+"=foobar"),
           (SEQUENCE_OUTPUT_LOCATION_OVERRIDE+"=http://seq")
       };
        */


        buildGenome(args);
    }

    public static GenomeManager.GenomeListItem buildGenome(String[] args) throws Exception {

        if (args.length == 0) {
            throw new RuntimeException("No command line parameters supplied!");
        }

        HashMap<String, String> commandLineValues = new HashMap();

        for (String parameter : args) {

            String optionAndValue[] = parameter.split("=");
            if (optionAndValue.length != 2) {
                throw new RuntimeException("Invalid command line option=value set:[" + parameter + "]");
            }
            if (!validCommandLineOptions.contains(optionAndValue[0])) {
                throw new RuntimeException("Invalid command line option:[" + optionAndValue[0] + "]");
            }
            commandLineValues.put(optionAndValue[0], optionAndValue[1]);
        }

        if (!commandLineValues.containsKey(GENOME_OUTPUT_LOCATION)) {
            throw new RuntimeException("No genome output location was supplied!");
        }
        if (!commandLineValues.containsKey(GENOME_DISPLAY_NAME)) {
            throw new RuntimeException("No Genome name was supplied!");
        }
        if (!commandLineValues.containsKey(SEQUENCE_FILE)) {
            throw new RuntimeException("No Sequence file was supplied!");
        }
        String genomeZipOutputLocation = commandLineValues.get(GENOME_OUTPUT_LOCATION); // Mandatory
        String genomeDisplayName = commandLineValues.get(GENOME_DISPLAY_NAME);//Mandatory
        String fastaFileName = commandLineValues.get(SEQUENCE_FILE); // Mandatory
        String genomeId = commandLineValues.get(GENOME_ID);
        String sequenceOutputLocationOverride = commandLineValues.get(SEQUENCE_OUTPUT_LOCATION_OVERRIDE);
        String cytobandFileName = commandLineValues.get(CYTOBAND_FILE);
        String refFlatFileName = commandLineValues.get(GENE_FILE);
        String genomeArchiveFilenamePrefix = commandLineValues.get(GENOME_ARCHIVE_FILENAME_PREFIX_OVERRIDE);

        genomeDisplayName = genomeDisplayName.trim();
        genomeId = genomeId.trim();

        if (genomeId == null || genomeId.length() == 0) {
            genomeId = GenomeManager.getInstance().generateGenomeKeyFromText(genomeDisplayName);
        }

        String genomeFileName = null;
        if (genomeArchiveFilenamePrefix != null) {
            genomeFileName = genomeArchiveFilenamePrefix + Globals.GENOME_FILE_EXTENSION;
        } else {
            genomeFileName = genomeId + Globals.GENOME_FILE_EXTENSION;
        }

        String relativeSequenceLocation = null;

        // Sequence location
        if (fastaFileName != null && !fastaFileName.trim().equals("")) {

            if (genomeArchiveFilenamePrefix != null) {
                relativeSequenceLocation = "/" + genomeArchiveFilenamePrefix;
            } else {
                relativeSequenceLocation = "/" + genomeId ;
            }
        }

        GenomeManager.GenomeListItem item = GenomeManager.getInstance().defineGenome(
                genomeZipOutputLocation,
                cytobandFileName,
                refFlatFileName,
                fastaFileName,
                null,
                relativeSequenceLocation,
                genomeDisplayName,
                genomeId,
                genomeFileName,
                null,
                sequenceOutputLocationOverride);

        return item;
    }
}
