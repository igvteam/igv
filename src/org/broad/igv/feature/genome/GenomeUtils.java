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

import org.broad.igv.Globals;
import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.*;

/**
 * Static utility functions for genome data-wrangling.
 *
 * @author jrobinso
 *         Date: 4/22/13
 *         Time: 1:27 PM
 */
public class GenomeUtils {


    public static void main(String[] args) throws IOException {

//        String directory = ".";
//        if (args.length > 0) {
//            directory = args[0];
//        }
//        String genomeList = "http://igv.broadinstitute.org/genomes/genomes.txt";
//        if (args.length > 1) {
//            genomeList = args[0];
//        }
//        exportAllChromSizes(new File(directory), genomeList);

        mergeINCDCNames(
                new File("genomes/alias/hg38_alias.tab"),
                new File("/Users/jrobinso/projects/INSDC/GCF_000001405.26.assembly.txt"),
                new File("/Users/jrobinso/projects/INSDC"));

    }


    /**
     * Create .chrom.sizes file for each genome found in the {@code genomeListPath}, and write it out to
     * {@code directory}
     *
     * @param directory
     * @param genomeListPath
     * @throws IOException
     */
    public static void exportAllChromSizes(File directory, String genomeListPath) throws IOException {

        //<Server-Side Genome List>
        // Human hg19	http://igv.broadinstitute.org/genomes/hg19.genome	hg19
        BufferedReader br = null;

        try {
            br = ParsingUtils.openBufferedReader(genomeListPath);
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                if (tokens.length > 2) {
                    String genomePath = tokens[1];
                    try {
                        Genome genome = GenomeManager.getInstance().loadGenome(genomePath, null);
                        exportChromSizes(directory, genome);

                    } catch (Exception e) {
                        System.err.println(e.toString());
                    }
                }
            }
        } finally {
            if (br != null) br.close();
        }

    }


    /**
     * Export a "chrom.sizes" file for the specified genome
     *
     * @param directory output directory
     * @param genome
     * @throws FileNotFoundException
     */
    public static void exportChromSizes(File directory, Genome genome) throws FileNotFoundException {


        String fn = genome.getId() + ".chrom.sizes";
        File file = new File(directory, fn);
        PrintWriter pw = null;

        try {
            pw = new PrintWriter(file);
            for (String chr : genome.getAllChromosomeNames()) {

                Chromosome chromosome = genome.getChromosome(chr);
                pw.println(chromosome.getName() + "\t" + chromosome.getLength());

            }
        } finally {
            if (pw != null) pw.close();
        }

    }

    /**
     * Merge chromosome names from an NCBI assembly.txt file with an existing IGV alias file
     *
     * @param aliasFile
     * @param assemblyFile
     */
    public static void mergeINCDCNames(File aliasFile, File assemblyFile, File outputDirectory) throws IOException {

        Map<String, Set<String>> aliasRows = new LinkedHashMap<String, Set<String>>();

        BufferedReader br = null;
        PrintWriter pw = null;

        // Build alias dictionary
        br = new BufferedReader(new FileReader(aliasFile));
        String nextLine;
        while ((nextLine = br.readLine()) != null) {
            String[] tokens = Globals.whitespacePattern.split(nextLine);
            HashSet<String> row = new LinkedHashSet<String>(Arrays.asList(tokens));
            for (String nm : tokens) {
                aliasRows.put(nm, row);
            }
        }
        br.close();

        // Loop through assembly file
        int[] chrIndeces = {0, 4, 6, 9};
        br = new BufferedReader(new FileReader(assemblyFile));
        boolean start = false;
        List<String> newRows = new ArrayList<String>();
        while ((nextLine = br.readLine()) != null) {
            if (start) {

                String[] tokens = Globals.tabPattern.split(nextLine);
                boolean foundRow = false;
                for (int i : chrIndeces) {
                    Set<String> row = aliasRows.get(tokens[i]);
                    if (row != null) {
                        for (int j : chrIndeces) {
                            if (!"na".equals(tokens[j])) {
                                row.add(tokens[j]);
                            }
                        }
                        foundRow = true;
                        break;
                    }
                }
                if (!foundRow) {
                    String newRow = tokens[chrIndeces[0]];
                    for (int i = 1; i < chrIndeces.length; i++) {
                        String chrNm = tokens[chrIndeces[i]];
                        if (!"na".equals(chrNm)) {
                            newRow += ("\t" + chrNm);
                        }
                    }
                    newRows.add(newRow);
                    System.out.println("New alias row: " + newRow);
                }

            } else if (nextLine.startsWith("# Sequence-Name")) {
                start = true;
            }

        }
        br.close();

        pw = new PrintWriter(new BufferedWriter(new FileWriter(new File(outputDirectory, aliasFile.getName()))));
        Set<Set<String>> output = new HashSet<Set<String>>();
        for (Set<String> row : aliasRows.values()) {
            if (row.size() == 0) continue;
            if (!output.contains(row)) {
                output.add(row);
                List<String> chrNames = new ArrayList<String>(row);
                pw.print(chrNames.get(0));
                for (int i = 1; i < chrNames.size(); i++) {
                    pw.print("\t" + chrNames.get(i));
                }
                pw.println();
            }
        }
        for (String row : newRows) {
            pw.println(row);
        }
        pw.close();

    }
}
