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

import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Static utility functions for genome data-wrangling.
 *
 * @author jrobinso
 *         Date: 4/22/13
 *         Time: 1:27 PM
 */
public class GenomeUtils {


    public static void main(String[] args) throws IOException {

        String directory = ".";
        if (args.length > 0) {
            directory = args[0];
        }
        String genomeList = "http://igv.broadinstitute.org/genomes/genomes.txt";
        if (args.length > 1) {
            genomeList = args[0];
        }
        exportAllChromSizes(new File(directory), genomeList);

    }


    /**
     * Create .chrom.sizes file for each genome found in the {@code genomeListPath}, and write it out to
     * {@code directory}
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

                        // Export aliases, if any exist
                        exportChrAliases(directory, genome);
                    } catch (Exception e) {
                        System.err.println(e.toString());
                    }
                }
            }
        } finally {
            if (br != null) br.close();
        }

    }

    private static void exportChrAliases(File directory, Genome genome) throws FileNotFoundException {

        String id = genome.getId();
        Map<String, String> chrAliasMap = genome.getChrAliasTable();
        if (chrAliasMap != null) {

            // Filter frivolous  and automatic entries.  Its not critical that this be complete,
            // but reduces unnecessary entries.

            Map<String, String> tmp = new HashMap<String, String>();
            Map<String, String> autoAliases = genome.getAutoAliases();

            for (Map.Entry<String, String> entry : chrAliasMap.entrySet()) {

                final String key = entry.getKey();
                final String value = entry.getValue();
                if (!(key.equals(value) || value.equals(autoAliases.get(key)))) {
                    tmp.put(key, value);
                }
            }
            chrAliasMap = tmp;

            if (chrAliasMap.size() > 0) {
                String fn = genome.getId() + "_alias.tab";
                File file = new File(directory, fn);
                PrintWriter pw = null;

                try {
                    pw = new PrintWriter(file);
                    for (Map.Entry<String, String> entry : chrAliasMap.entrySet()) {
                        pw.println(entry.getKey() + "\t" + entry.getValue());
                    }
                } finally {
                    if (pw != null) pw.close();
                }
            }
        }
    }

    /**
     * Export a "chrom.sizes" file for the specified genome
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
}
