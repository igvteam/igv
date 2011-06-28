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

package org.broad.igv.util;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Jun 27, 2011
 */
public class GenomeUtils {

    static String rootDir = "/Volumes/seq_transcriptome/public/IGV/PORTALS/Schizos/genomes/";

    public static void main(String[] args) throws IOException {

        File cytoFile = new File(rootDir + "/SJ4.tmp", "Sjaponicus4_cytoband.txt");
        File aliasFile = new File(rootDir + "/SJ4.tmp", "SJ4_alias.tab");
        File seqDir =     new File("/Volumes/seq_transcriptome/public/IGV/seq/Sjaponicus4.genome_seq");
        renameContigs(aliasFile, cytoFile, seqDir);
    }


    public static void renameContigs(File aliasFile, File cytobandFile, File seqDir) throws IOException {

        BufferedReader aliasReader = null;
        BufferedReader cytoReader = null;
        PrintWriter cytoWriter = null;

        File tmpCytoFile = new File(cytobandFile.getAbsolutePath() + ".tmp");

        File backup = new File(cytobandFile.getAbsolutePath() + ".bak");
        copyFile(cytobandFile, backup);

        try {
            aliasReader = new BufferedReader(new FileReader(aliasFile));
            cytoReader = new BufferedReader(new FileReader(cytobandFile));
            cytoWriter = new PrintWriter(new BufferedWriter(new FileWriter(tmpCytoFile)));

            Map<String, String> map = new HashMap();
            String nextLine;
            while ((nextLine = aliasReader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                if (tokens.length > 1) {
                    map.put(tokens[0], tokens[1]);

                    
                } else {
                    System.out.println(nextLine);
                }

            }

            while ((nextLine = cytoReader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                String newName = map.get(tokens[0]);
                if (newName != null) {
                    String oldName = tokens[0];
                    tokens[0] = newName;

                    File seqFile = new File(seqDir, oldName + ".txt");
                    File newSeqfile = new File(seqDir, newName + ".txt");
                    copyFile(seqFile, newSeqfile);
                }
                cytoWriter.print(tokens[0]);
                for (int i = 1; i < tokens.length; i++) {
                    cytoWriter.print("\t" + tokens[i]);
                }
                cytoWriter.println();
            }

        } finally {
            if (aliasReader != null) aliasReader.close();
            if (cytoReader != null) cytoReader.close();
            if (cytoWriter != null) cytoWriter.close();
        }
        copyFile(tmpCytoFile, cytobandFile);



    }


    private static void copyFile(File fromFile, File toFile) throws IOException {

        System.out.println("cp " + fromFile.getName() + " " + toFile.getName());

        Runtime.getRuntime().exec("cp " + fromFile.getAbsolutePath() + " " + toFile.getAbsolutePath());

    }
}
