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

package org.broad.igv.hic.tools;

import jargs.gnu.CmdLineParser;
import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.hic.data.Chromosome;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

/**
 * @author Jim Robinson
 * @date 9/16/11
 */
public class HiCTools {


    public static void main(String[] argv) throws IOException, CmdLineParser.UnknownOptionException, CmdLineParser.IllegalOptionValueException {

        if (argv.length < 4) {
            System.out.println("Usage: hictools pre <options> <inputFile> <outputFile> <genomeID>");
            System.exit(0);
        }

        Globals.setHeadless(true);

        CommandLineParser parser = new CommandLineParser();
        parser.parse(argv);
        String[] args = parser.getRemainingArgs();

        if (args[0].equals("sort")) {
            AlignmentsSorter.sort(args[1], args[2], null);
        } else if (args[0].equals("pre")) {

            String genomeId = args[3];
            List<Chromosome> chromosomes = loadChromosomes(genomeId);

            long genomeLength = 0;
            for (Chromosome c : chromosomes) {
                if (c != null)
                    genomeLength += c.getSize();
            }
            chromosomes.set(0, new Chromosome(0, "All", (int) (genomeLength / 1000)));

            String[] tokens = args[1].split(",");
            List<String> files = new ArrayList(tokens.length);
            for (String f : tokens) {
                files.add(f);
            }

            Preprocessor preprocessor = new Preprocessor(new File(args[2]), chromosomes);

            preprocessor.setIncludedChromosomes(parser.getchromosomeOption());
            preprocessor.setCountThreshold(parser.getCountThresholdOption());
            preprocessor.setDiagonalsOnly(parser.getDiagonalsOption());

            preprocessor.preprocess(files);
        }
    }

    /**
     * Load chromosomes
     *
     * @param idOrFile
     * @return
     */
    public static List<Chromosome> loadChromosomes(String idOrFile) throws IOException {

        InputStream is = null;

        //First assume this is an ID an
        try {
            is = HiCTools.class.getResourceAsStream(idOrFile + ".chrom.sizes");
            if (is == null) {
                // Not an ID,  see if its a file
                File file = new File(idOrFile);
                if (file.exists()) {
                    is = new FileInputStream(file);
                } else {
                    throw new FileNotFoundException("Could not find chrom sizes file for: " + idOrFile);
                }

            }

            List<Chromosome> chromosomes = new ArrayList();
            chromosomes.add(0, null);   // Index 0 reserved for "whole genome" psuedo-chromosome

            Pattern pattern = Pattern.compile("\t");
            BufferedReader reader = new BufferedReader(new InputStreamReader(is));
            String nextLine;
            long genomeLength = 0;
            int idx = 1;

            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = pattern.split(nextLine);
                if (tokens.length == 2) {
                    String name = tokens[0];
                    int length = Integer.parseInt(tokens[1]);
                    genomeLength += length;
                    chromosomes.add(idx, new Chromosome(idx, name, length));
                    idx++;
                } else {
                    System.out.println("Skipping " + nextLine);
                }
            }

            // Add the "psuedo-chromosome" All, representing the whole genome.  Units are in kilo-bases
            chromosomes.set(0, new Chromosome(0, "All", (int) (genomeLength / 1000)));


            return chromosomes;
        } finally {
            if (is != null) is.close();
        }

    }


    /**
     * Convert a BAM file containing paried-end tags to the ascii "pair" format used for HiC.
     *
     * @param inputBam
     * @param outputFile
     * @throws IOException
     */
    public static void filterBam(String inputBam, String outputFile, List<Chromosome> chromosomes) throws IOException {

        CloseableIterator<Alignment> iter = null;
        AlignmentReader reader = null;
        PrintWriter pw = null;

        HashSet allChroms = new HashSet(chromosomes);

        try {
            pw = new PrintWriter(new FileWriter(outputFile));
            reader = AlignmentReaderFactory.getReader(inputBam, false);
            iter = reader.iterator();
            while (iter.hasNext()) {

                Alignment alignment = iter.next();
                ReadMate mate = alignment.getMate();

                // Filter unpaired and "normal" pairs.  Only interested in abnormals
                if (alignment.isPaired() &&
                        alignment.isMapped() &&
                        alignment.getMappingQuality() > 10 &&
                        mate != null &&
                        mate.isMapped() &&
                        allChroms.contains(alignment.getChr()) &&
                        allChroms.contains(mate.getChr()) &&
                        (!alignment.getChr().equals(mate.getChr()) || alignment.getInferredInsertSize() > 1000)) {

                    // Each pair is represented twice in the file,  keep the record with the "leftmost" coordinate
                    if (alignment.getStart() < mate.getStart()) {
                        String strand = alignment.isNegativeStrand() ? "-" : "+";
                        String mateStrand = mate.isNegativeStrand() ? "-" : "+";
                        pw.println(alignment.getReadName() + "\t" + alignment.getChr() + "\t" + alignment.getStart() +
                                "\t" + strand + "\t.\t" + mate.getChr() + "\t" + mate.getStart() + "\t" + mateStrand);
                    }
                }

            }
        } finally {
            pw.close();
            iter.close();
            reader.close();
        }
    }


    static class CommandLineParser extends CmdLineParser {
        private CmdLineParser.Option diagonalsOption = null;
        private CmdLineParser.Option chromosomeOption = null;
        private CmdLineParser.Option countThresholdOption = null;

        CommandLineParser() {
            diagonalsOption = addBooleanOption('d', "diagonals");
            chromosomeOption = addStringOption('c', "chromosomes");
            countThresholdOption = addIntegerOption('t', "countThreshold");
        }

        boolean getDiagonalsOption() {
            Object opt = getOptionValue(diagonalsOption);
            return opt == null ? false : ((Boolean) opt).booleanValue();
        }

        Set<String> getchromosomeOption() {
            Object opt = getOptionValue(chromosomeOption);
            if (opt != null) {
                String[] tokens = opt.toString().split(",");
                return new HashSet<String>(Arrays.asList(tokens));
            } else {
                return null;
            }
        }

        int getCountThresholdOption() {
            Object opt = getOptionValue(countThresholdOption);
            return opt == null ? -1 : ((Number) opt).intValue();
        }
    }


}
