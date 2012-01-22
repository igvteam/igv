package org.broad.igv.dev.affective;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeListItem;

import java.io.*;
import java.util.Date;

/**
 * @author Jim Robinson
 * @date 1/11/12
 */
public class AffectiveUtils {

    public static int POINTS_PER_SECOND = 8;

    // Length of day in hours
    public static int DAY_LENGTH_HOURS = 8;

    // Start time in "data points" units (8:00 AM)
    public static int START_TIME_HR = 8;
    public static int START_TIME = START_TIME_HR * 60 * 60 * POINTS_PER_SECOND;
    public static final GenomeListItem GENOME_DESCRIPTOR = new GenomeListItem("Affective", "", "affective", false);
    private static AffectiveGenome genome;


    public static void createCytoband(String[] args) throws IOException {

        File root = new File("/Users/jrobinso/IGV/Miriah/Participant3");

        int start = 0;
        for (File dir : root.listFiles()) {
            final File[] files = dir.listFiles();
            if (dir.isDirectory()) {
                for (File f : files) {
                    if (f.getName().endsWith(".csv")) {
                        int nLines = 0;
                        BufferedReader br = new BufferedReader(new FileReader(f));
                        while (br.readLine() != null) {
                            nLines++;
                        }
                        System.out.println(dir.getName().replace("Day ", "") + "\t" + 0 + "\t" + nLines);
                        break;
                    }
                }
            }
        }
    }

    public static Genome getGenome() {

        genome = new AffectiveGenome();
        genome.addChromosome (new AffectiveChromosome(new Date()));
        return genome;
    }
}
