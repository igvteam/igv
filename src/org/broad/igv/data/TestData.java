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
package org.broad.igv.data;

import org.apache.log4j.Logger;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * @author eflakes
 */
public class TestData {

    private static Logger logger = Logger.getLogger(TestData.class);
    private static String genomeId = "hg18";

    public static void main(String[] args) {


        generateIGVTestFiles(new File("/Users/jrobinso/IGVTestData"));
    }

    static public void generateAllTestFiles(File fileLocation) {
        generateCNTestFiles(fileLocation);
        generateGctTestFiles(fileLocation);
    }

    static public void generateIGVTestFiles(File cnFileLocation) {
        int numberOfSamples = 100;
        int fileLength = 10000;
        Genome genome = IGV.getInstance().getGenomeManager().getGenome(genomeId);
        if (genome == null) {
            throw new RuntimeException("Unknown genome: " + genomeId);
        }
        List<String> chromosomeList =
                new ArrayList<String>(genome.getChromosomeNames());
        generateIGVTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);


    }

    /**
     * Generate SNP files
     *
     * @param cnFileLocation
     */
    static public void generateCNTestFiles(File cnFileLocation) {

        Genome genome = IGV.getInstance().getGenomeManager().getGenome(genomeId);
        if (genome == null) {
            throw new RuntimeException("Unknown genome: " + genomeId);
        }
        List<String> chromosomeList =
                new ArrayList<String>(genome.getChromosomeNames());

        // Fix rows at 1 and change sample


        for (int numberOfSamples : new int[]{10, 100, 1000}) {
            int fileLength = 1;
            generateCNTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);
        }

        // Change rows and fix samples at 1
        for (int fileLength : new int[]{10, 100, 1000, 10000, 100000, 500000}) {
            int numberOfSamples = 1;
            generateCNTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);
        }

        for (int nSamples : new int[]{10, 100}) {
            for (int fileLength : new int[]{1000, 10000}) {
                int numberOfSamples = 10;
                generateCNTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);
            }
        }

        // Fix rows at 25000 and change samples
        int numberOfSamples = 1;
        int fileLength = 25000;
        generateCNTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);
        numberOfSamples = 100;
        generateCNTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);
        numberOfSamples = 300;
        generateCNTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);

        // Fix rows at 250000 and 300 samples
        fileLength = 250000;
        numberOfSamples = 300;
        generateCNTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);

        // Fix rows at 500000 and 300 samples
        fileLength = 500000;
        numberOfSamples = 300;
        generateCNTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);

    }

    /**
     * Generate a single large SNP file to stress IGV's processing limit
     *
     * @param cnFileLocation
     */
    static public void generateVeryLargeCNTestFile(File cnFileLocation) {

        Genome genome = IGV.getInstance().getGenomeManager().getGenome(genomeId);
        if (genome == null) {
            throw new RuntimeException("Unknown genome: " + genomeId);
        }
        List<String> chromosomeList =
                new ArrayList<String>(genome.getChromosomeNames());

        int numberOfSamples = 96;
        //int fileLength = 800000;
        int fileLength = 1800000;

        generateCNTestFile(genome, chromosomeList, cnFileLocation, fileLength, numberOfSamples);
    }


    /**
     * Generate SNP files
 
     */
    static private void generateCNTestFile(Genome genome, List<String> chromosomeList,
                                           File outputLocation, int fileLength, int numberOfSamples) {

        // File length has to be atleast the number of chromosomes 
        // we are working with otherwise the spec'd algorithm won't
        // work (not enough data)
        if (fileLength < chromosomeList.size()) {
            fileLength = chromosomeList.size();
        }

        int snpIdentifier = 0;
        int chunks = fileLength / chromosomeList.size();
        if (chunks < 1) {
            chunks = 1;
        }
        int divisor = chunks - 1;
        float lowestSampleValue = 0.0f;
        float highestSampleValue = lowestSampleValue + 4.0f;
        float sampleValueIncrement =
                (highestSampleValue - lowestSampleValue) / chromosomeList.size();

        GeneManager geneManager = GeneManager.getGeneManager(genome.getId());

        geneManager.sortGeneLists();
        sortChromosomeList(chromosomeList);

        // Calculate sample value increment
        if (divisor > 0) {
            sampleValueIncrement = (highestSampleValue - lowestSampleValue) / divisor;
        }

        File file = new File(outputLocation, "igv_cn_" + fileLength + "_" + numberOfSamples + ".cn");
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(file);
            StringBuffer record = new StringBuffer();

            // Create SNP Header
            record.append("SNP");
            record.append("\t");
            record.append("Chromosome");
            record.append("\t");
            record.append("PhysicalPosition");
            for (int i = 0; i < numberOfSamples; i++) {
                record.append("\t");
                record.append("Value");
                record.append((i + 1));
            }
            record.append("\n");

            writer.write(record.toString());
            record.delete(0, record.length());

            // Process each Chromosome one at a time
            float sampleValue = lowestSampleValue;
            for (String chromosomeName : chromosomeList) {

                Chromosome chromosome = genome.getChromosome(chromosomeName);

                float chromosomeChunkIncrement = chromosome.getLength();
                if (divisor > 0) {
                    chromosomeChunkIncrement = chromosome.getLength() / divisor;
                }

                chromosomeName = chromosome.getName();
                chromosomeName = chromosomeName.substring(3, chromosomeName.length());

                for (float i = 0; i < chunks; i++) {

                    long chunkStart =
                            (int) Math.floor(i * chromosomeChunkIncrement);

                    // Create Chromosome Record
                    record.append("SNP_");
                    record.append(++snpIdentifier);
                    record.append("\t");
                    record.append(chromosomeName);
                    record.append("\t");
                    record.append(chunkStart);
                    for (int j = 0; j < numberOfSamples; j++) {
                        record.append("\t");
                        record.append(sampleValue);
                    }
                    writer.println(record.toString());
                    record.delete(0, record.length());
                    sampleValue += sampleValueIncrement;
                }

                // Reset sample value to lowest if we actually have enough 
                // row data to calculate it; otherwise, we must be working 
                // with 1 loop over each chromosome so each will only have 
                // a single color
                if (divisor > 0) {
                    sampleValue = lowestSampleValue;
                }
            }

            writer.close();
            writer = null;
        } catch (Exception e) {
            logger.error("Could not write test files", e);
        } finally {

            if (writer != null) {
                writer.close();
                writer = null;
            }
        }
    }

    static private void generateIGVTestFile(Genome genome, List<String> chromosomeList,
                                            File outputLocation, int fileLength, int numberOfSamples) {

        // File length has to be atleast the number of chromosomes 
        // we are working with otherwise the spec'd algorithm won't
        // work (not enough data)
        if (fileLength < chromosomeList.size()) {
            fileLength = chromosomeList.size();
        }

        int snpIdentifier = 0;
        int chunks = fileLength / chromosomeList.size();
        if (chunks < 1) {
            chunks = 1;
        }
        int divisor = chunks - 1;
        float lowestSampleValue = 0.0f;
        float highestSampleValue = lowestSampleValue + 4.0f;
        float sampleValueIncrement =
                (highestSampleValue - lowestSampleValue) / chromosomeList.size();

        GeneManager geneManager = GeneManager.getGeneManager(genome.getId());

        geneManager.sortGeneLists();
        sortChromosomeList(chromosomeList);

        // Calculate sample value increment
        if (divisor > 0) {
            sampleValueIncrement = (highestSampleValue - lowestSampleValue) / divisor;
        }

        File file = new File(outputLocation, "igv_" + fileLength + "_" + numberOfSamples + ".igv");
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(file);
            StringBuffer record = new StringBuffer();

            // Write a few comment lines
            writer.println("#comment 1");
            writer.println("#comment 2");

            // Create SNP Header
            record.append("Chromosome\tStart\tEnd\tName");
            for (int i = 0; i < numberOfSamples; i++) {
                record.append("\t");
                record.append("Value");
                record.append((i + 1));
            }
            record.append("\n");

            writer.write(record.toString());
            record.delete(0, record.length());

            // Process each Chromosome one at a time
            float sampleValue = lowestSampleValue;
            for (String chromosomeName : chromosomeList) {

                Chromosome chromosome = genome.getChromosome(chromosomeName);

                float chromosomeChunkIncrement = chromosome.getLength();
                if (divisor > 0) {
                    chromosomeChunkIncrement = chromosome.getLength() / divisor;
                }

                chromosomeName = chromosome.getName();
                chromosomeName = chromosomeName.substring(3, chromosomeName.length());

                for (float i = 0; i < chunks; i++) {

                    long chunkStart =
                            (int) Math.floor(i * chromosomeChunkIncrement);

                    // Create Chromosome Record
                    record.append(chromosomeName);
                    record.append("\t");
                    record.append(chunkStart);
                    record.append("\t");
                    record.append(chunkStart + 1000);
                    record.append("\t");
                    record.append("Probe_");
                    record.append(++snpIdentifier);
                    for (int j = 0; j < numberOfSamples; j++) {
                        record.append("\t");
                        record.append(sampleValue);
                    }
                    writer.println(record.toString());
                    record.delete(0, record.length());
                    sampleValue += sampleValueIncrement;
                }

                // Reset sample value to lowest if we actually have enough 
                // row data to calculate it; otherwise, we must be working 
                // with 1 loop over each chromosome so each will only have 
                // a single color
                if (divisor > 0) {
                    sampleValue = lowestSampleValue;
                }
            }

            writer.close();
            writer = null;
        } catch (Exception e) {
            logger.error("Could not write test files", e);
        } finally {

            if (writer != null) {
                writer.close();
                writer = null;
            }
        }
    }

    static public void generateGctTestFiles(File outputLocation) {

        //ReferenceFrame context = ReferenceFrame.getInstance();
        Genome genome = IGV.getInstance().getGenomeManager().getGenome(genomeId);
        if (genome == null) {
            throw new RuntimeException("Unknown genome: " + genomeId);
        }
        List<String> chromosomeList =
                new ArrayList<String>(genome.getChromosomeNames());

        generateGctGeneTestFile(genome, chromosomeList, outputLocation);
        generateGctChromosomeTestFile(genome, chromosomeList, outputLocation);
    }

    static public void generateGctTestFiles(Genome genome, List<String> chromosomeList, File outputLocation) {
        generateGctGeneTestFile(genome, chromosomeList, outputLocation);
        generateGctChromosomeTestFile(genome, chromosomeList, outputLocation);
    }

    static private void generateGctGeneTestFile(Genome genome, List<String> chromosomeList, File outputLocation) {

        float lowestSampleValue = -1.5f;
        float highestSampleValue = lowestSampleValue + 3.0f;

        GeneManager geneManager = GeneManager.getGeneManager(genome.getId());

        geneManager.sortGeneLists();
        sortChromosomeList(chromosomeList);

        int divisor = chromosomeList.size() - 1;
        float sampleValueIncrement = 0;
        if (divisor > 0) {
            sampleValueIncrement = (highestSampleValue - lowestSampleValue) / divisor;
        }

        File file = new File(outputLocation, "igv_genes.gct");
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(file);
            int lengthOfFile = countTotalGenes(genome.getId(), chromosomeList);
            int numberOfSamples = 1;
            StringBuffer record = new StringBuffer();

            // Create Gene Header
            record.append("#1.2");
            record.append("\n");
            record.append(lengthOfFile);
            record.append("\t");
            record.append(numberOfSamples);
            record.append("\n");
            record.append("NAME");
            record.append("\t");
            record.append("DESCRIPTION");
            record.append("\t");
            record.append("SAMPLE");
            record.append("\n");

            writer.write(record.toString());
            record.delete(0, record.length());

            float sampleValue = lowestSampleValue;
            for (String chromosome : chromosomeList) {

                List<org.broad.tribble.Feature> genes =
                        geneManager.getGenesForChromosome(chromosome);

                for (org.broad.tribble.Feature f : genes) {

                    IGVFeature gene = (IGVFeature) f;
                    // Create Gene Record
                    record.append(gene.getName());
                    record.append("\t");
                    record.append("|@");
                    record.append(chromosome);
                    record.append(":");
                    record.append(gene.getStart());
                    record.append("-");
                    record.append(gene.getEnd());
                    record.append("|");
                    record.append("\t");
                    record.append(sampleValue);
                    writer.println(record.toString());
                    record.delete(0, record.length());
                }
                sampleValue += sampleValueIncrement;
            }

            writer.close();
            writer = null;
        } catch (Exception e) {
            logger.error("Could not write test files", e);
        } finally {

            if (writer != null) {
                writer.close();
                writer = null;
            }
        }
    }

    static private void generateGctChromosomeTestFile(Genome genome, List<String> chromosomeList, File outputLocation) {

        int chunks = 100;
        int divisor = chunks - 1;
        if (chunks < 1) {
            chunks = 1;
        }
        float lowestSampleValue = -1.5f;
        float highestSampleValue = lowestSampleValue + 3.0f;
        float sampleValueIncrement =
                (highestSampleValue - lowestSampleValue) / chromosomeList.size();

        GeneManager geneManager = GeneManager.getGeneManager(genome.getId());

        geneManager.sortGeneLists();
        sortChromosomeList(chromosomeList);

        // Calculate sample value increment
        if (divisor > 0) {
            sampleValueIncrement = (highestSampleValue - lowestSampleValue) / divisor;
        }

        File file = new File(outputLocation, "igv_chromosome.gct");
        PrintWriter writer = null;
        try {
            writer = new PrintWriter(file);
            int lengthOfFile = chromosomeList.size() * chunks;
            int numberOfSamples = 1;
            StringBuffer record = new StringBuffer();

            // Create Gene Header
            record.append("#1.2");
            record.append("\n");
            record.append(lengthOfFile);
            record.append("\t");
            record.append(numberOfSamples);
            record.append("\n");
            record.append("NAME");
            record.append("\t");
            record.append("DESCRIPTION");
            record.append("\t");
            record.append("SAMPLE");
            record.append("\n");

            writer.write(record.toString());
            record.delete(0, record.length());

            // Process each Chromosome one at a time
            float sampleValue = lowestSampleValue;
            for (String chromosomeName : chromosomeList) {

                Chromosome chromosome = genome.getChromosome(chromosomeName);

                float chromosomeChunkIncrement = chromosome.getLength();
                if (divisor > 0) {
                    chromosomeChunkIncrement = chromosome.getLength() / divisor;
                }

                for (float i = 0; i < chunks; i++) {

                    long chunkStart =
                            (int) Math.floor(i * chromosomeChunkIncrement);
                    long chunkEnd =
                            (int) Math.ceil(chunkStart + chromosomeChunkIncrement - 1);

                    // Create Chromosome Record
                    record.append("Probe_" + i);
                    record.append("\t");
                    record.append("|@");
                    record.append(chromosome.getName());
                    record.append(":");
                    record.append(chunkStart);
                    record.append("-");
                    record.append(chunkEnd);
                    record.append("|");
                    record.append("\t");
                    record.append(sampleValue);
                    writer.println(record.toString());
                    record.delete(0, record.length());
                    sampleValue += sampleValueIncrement;
                }

                // Reset sample value to lowest if we actually have enough 
                // row data to calculate it; otherwise, we must be working 
                // with 1 loop over each chromosome so each will only have 
                // a single color
                if (divisor > 0) {
                    sampleValue = lowestSampleValue;
                }
            }

            writer.close();
            writer = null;
        } catch (Exception e) {
            logger.error("Could not write test files", e);
        } finally {

            if (writer != null) {
                writer.close();
                writer = null;
            }
        }
    }

    static public int countTotalGenes(String genome, List<String> chromosomeList) {

        int totalGenes = 0;
        GeneManager geneManager = GeneManager.getGeneManager(genome);

        for (String chromosome : chromosomeList) {
            List<org.broad.tribble.Feature> genes = geneManager.getGenesForChromosome(chromosome);
            for (org.broad.tribble.Feature gene : genes) {
                ++totalGenes;
            }
        }
        return totalGenes;
    }

    static private void sortChromosomeList(List<String> chromosomeList) {

        Collections.sort(chromosomeList, new Comparator() {

            public int compare(Object arg0, Object arg1) {

                String name1 = (String) arg0;
                String name2 = (String) arg1;
                if (name1.startsWith("chr") && name2.startsWith("chr")) {

                    try {
                        name1 = name1.substring(3).trim();
                        name2 = name2.substring(3).trim();
                        Integer integer1 = Integer.parseInt(name1);
                        Integer integer2 = Integer.parseInt(name2);
                        int returnValue = integer1.compareTo(integer2);

                        return returnValue;

                    } catch (NumberFormatException numberFormatException) {
                        return 0;
                    }
                } else {
                    throw new RuntimeException("Found an invalid chromosome name: " + name1 + " or " + name2);
                }
            }
        });
    }
}
