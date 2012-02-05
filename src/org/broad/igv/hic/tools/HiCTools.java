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

import net.sf.samtools.util.CloseableIterator;
import org.broad.igv.Globals;
import org.broad.igv.hic.data.Chromosome;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.reader.AlignmentReader;
import org.broad.igv.sam.reader.AlignmentReaderFactory;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * @author Jim Robinson
 * @date 9/16/11
 */
public class HiCTools {


    public static void main(String[] args) throws IOException {

        if (args.length < 3) {
            System.out.println("Usage: hictools <command> <inputFile> <outputFile>");
        }

        Globals.setHeadless(true);

        if (args[0].equals("sort")) {
            AlignmentsSorter.sort(args[1], args[2], null);
        } else if (args[0].equals("pre") || args[0].equals("bamToPairs")) {
            String genomeId = args[3];
            if (genomeId.equals("hg18")) {
                chromosomes = hg18Chromosomes;
            } else if (genomeId.equals("b37")) {
                chromosomes = b37Chromosomes;
            } else if (genomeId.equals("chr14")) {
                chromosomes = chr14;
            } else if (genomeId.equals("mm9")) {
                chromosomes = mm9Chromosomes;
            }else {
                chromosomes = dmelChromosomes;
            }

            if (args[0].equals("bamToPairs")) {
                filterBam(args[1], args[2]);
            } else {
                // Preprocess
                long genomeLength = 0;
                for (Chromosome c : chromosomes) {
                    if (c != null)
                        genomeLength += c.getSize();
                }
                chromosomes[0] = new Chromosome(0, "All", (int) (genomeLength / 1000));

                chromosomeOrdinals = new Hashtable();
                for (int i = 0; i < chromosomes.length; i++) {
                    chromosomeOrdinals.put(chromosomes[i].getName(), i);
                }

                String[] tokens = args[1].split(",");
                List<String> files = new ArrayList(tokens.length);
                for (String f : tokens) {
                    files.add(f);
                }

                (new Preprocessor(new File(args[2]))).preprocess(files, genomeId);
            }
        }
    }

    public static void filterBam(String inputBam, String outputFile) throws IOException {

        CloseableIterator<Alignment> iter = null;
        AlignmentReader reader = null;
        PrintWriter pw = null;

        HashSet allChroms = new HashSet(Arrays.asList(chromosomes));

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

    public static Chromosome[] chromosomes;
    public static Map<String, Integer> chromosomeOrdinals;

    // Hardcoded chromosomes for dMel r4.2.1
    public final static Chromosome[] dmelChromosomes = new Chromosome[]{
            null,                                 // Placeholder for whole genome
            new Chromosome(1, "2L", 22407834),
            new Chromosome(2, "2R", 20766785),
            new Chromosome(3, "3L", 23771897),
            new Chromosome(4, "3R", 27905053),
            new Chromosome(5, "4", 1281640),
            new Chromosome(6, "X", 22224390)
    };

    public final static Chromosome[] hg18Chromosomes = new Chromosome[]{
            null,                               // Placeholder for whole genome
            new Chromosome(1, "1", 247249719),
            new Chromosome(2, "2", 242951149),
            new Chromosome(3, "3", 199501827),
            new Chromosome(4, "4", 191273063),
            new Chromosome(5, "5", 180857866),
            new Chromosome(6, "6", 170899992),
            new Chromosome(7, "7", 158821424),
            new Chromosome(8, "8", 146274826),
            new Chromosome(9, "9", 140273252),
            new Chromosome(10, "10", 135374737),
            new Chromosome(11, "11", 134452384),
            new Chromosome(12, "12", 132349534),
            new Chromosome(13, "13", 114142980),
            new Chromosome(14, "14", 106368585),
            new Chromosome(15, "15", 100338915),
            new Chromosome(16, "16", 88827254),
            new Chromosome(17, "17", 78774742),
            new Chromosome(18, "18", 76117153),
            new Chromosome(19, "19", 63811651),
            new Chromosome(20, "20", 62435964),
            new Chromosome(21, "21", 46944323),
            new Chromosome(22, "22", 49691432),
            new Chromosome(23, "23", 154913754),
            new Chromosome(24, "24", 57772954)
    };


    public final static Chromosome[] b37Chromosomes = new Chromosome[]{
            null,                               // Placeholder for whole genome
            new Chromosome(1, "1", 249250621),
            new Chromosome(2, "2", 243199373),
            new Chromosome(3, "3", 198022430),
            new Chromosome(4, "4", 191154276),
            new Chromosome(5, "5", 180915260),
            new Chromosome(6, "6", 171115067),
            new Chromosome(7, "7", 159138663),
            new Chromosome(8, "8", 146364022),
            new Chromosome(9, "9", 141213431),
            new Chromosome(10, "10", 135534747),
            new Chromosome(11, "11", 135006516),
            new Chromosome(12, "12", 133851895),
            new Chromosome(13, "13", 115169878),
            new Chromosome(14, "14", 107349540),
            new Chromosome(15, "15", 102531392),
            new Chromosome(16, "16", 90354753),
            new Chromosome(17, "17", 81195210),
            new Chromosome(18, "18", 78077248),
            new Chromosome(19, "19", 59128983),
            new Chromosome(20, "20", 63025520),
            new Chromosome(21, "21", 48129895),
            new Chromosome(22, "22", 51304566),
            new Chromosome(23, "X", 155270560),
            new Chromosome(24, "Y", 59373566)
    };

    public final static Chromosome[] mm9Chromosomes = new Chromosome[]{
            null,                               // Placeholder for whole genome
            new Chromosome(1, "chr1", 197195432),
            new Chromosome(2, "chr2", 181748087),
            new Chromosome(3, "chr3", 159599783),
            new Chromosome(4, "chr4", 155630120),
            new Chromosome(5, "chr5", 152537259),
            new Chromosome(6, "chr6", 149517037),
            new Chromosome(7, "chr7", 152524553),
            new Chromosome(8, "chr8", 131738871),
            new Chromosome(9, "chr9", 124076172),
            new Chromosome(10, "chr10", 129993255),
            new Chromosome(11, "chr11", 121843856),
            new Chromosome(12, "chr12", 121257530),
            new Chromosome(13, "chr13", 120284312),
            new Chromosome(14, "chr14", 125194864),
            new Chromosome(15, "chr15", 103494974),
            new Chromosome(16, "chr16", 98319150),
            new Chromosome(17, "chr17", 95272651),
            new Chromosome(18, "chr18", 90772031),
            new Chromosome(19, "chr19", 61342430),
            new Chromosome(23, "chrX", 166650296),
            new Chromosome(24, "chrY", 15902555)
    };

    // For testing
    public final static Chromosome[] chr14 = new Chromosome[]{
            null,                               // Placeholder for whole genome
            new Chromosome(1, "14", 107349540)
    };


}
