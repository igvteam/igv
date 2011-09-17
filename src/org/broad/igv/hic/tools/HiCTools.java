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

import org.broad.igv.hic.data.Chromosome;

import java.io.File;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Map;

/**
 * @author Jim Robinson
 * @date 9/16/11
 */
public class HiCTools {


    public static void main(String [] args) throws IOException {

        if(args.length < 3) {
            System.out.println("Usage: hictools <command> <inputFile> <outputFile>");
        }

        if(args[0].equals("sort")) {
            AlignmentsSorter.sort(args[1], args[2], null);
        }
        else if(args[0].equals("pre")) {
            String genomeId = args[3];
            if(genomeId.equals("hg18")) {
                chromosomes = hg18Chromosomes;
            }
            else {
                chromosomes = dmelChromosomes;
            }
            chromosomeOrdinals = new Hashtable();
            for(int i=0; i<chromosomes.length; i++) {
                chromosomeOrdinals.put(chromosomes[i].getName(), i);
            }
            (new Preprocessor(new File(args[2]))).preprocess(new File(args[1]), genomeId);
        }
    }

    public static Chromosome[] chromosomes;
    public static Map<String, Integer> chromosomeOrdinals;

    // Hardcoded chromosomes for dMel r4.2.1
    public final static Chromosome[] dmelChromosomes = new Chromosome[]{
            new Chromosome(0, "2L", 22407834),
            new Chromosome(1, "2R", 20766785),
            new Chromosome(2, "3L", 23771897),
            new Chromosome(3, "3R", 27905053),
            new Chromosome(4, "4", 1281640),
            new Chromosome(5, "X", 22224390),
            new Chromosome(6, "U", 11561901)
    };

    public final static Chromosome[] hg18Chromosomes = new Chromosome[]{
            new Chromosome(0, "1", 247249719),
            new Chromosome(1, "2", 242951149),
            new Chromosome(2, "3", 199501827),
            new Chromosome(3, "4", 191273063),
            new Chromosome(4, "5", 180857866),
            new Chromosome(5, "6", 170899992),
            new Chromosome(6, "7", 158821424),
            new Chromosome(7, "8", 146274826),
            new Chromosome(8, "9", 140273252),
            new Chromosome(9, "10", 135374737),
            new Chromosome(10, "11", 134452384),
            new Chromosome(11, "12", 132349534),
            new Chromosome(12, "13", 114142980),
            new Chromosome(13, "14", 106368585),
            new Chromosome(14, "15", 100338915),
            new Chromosome(15, "16", 88827254),
            new Chromosome(16, "17", 78774742),
            new Chromosome(17, "18", 76117153),
            new Chromosome(18, "19", 63811651),
            new Chromosome(19, "20", 62435964),
            new Chromosome(20, "21", 46944323),
            new Chromosome(21, "22", 49691432),
            new Chromosome(22, "23", 154913754),
            new Chromosome(23, "24", 57772954),
            new Chromosome(0, "0", 16571),
    };




}
