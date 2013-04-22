package org.broad.igv.feature.genome;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.util.ParsingUtils;

import java.io.*;

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
     * Export a "chrom.sizes" file
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
