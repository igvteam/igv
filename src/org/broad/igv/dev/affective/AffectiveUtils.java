package org.broad.igv.dev.affective;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.LineplotRenderer;
import org.broad.igv.renderer.XYPlotRenderer;
import org.broad.igv.track.Track;

import java.awt.*;
import java.io.*;
import java.util.Date;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 1/11/12
 */
public class AffectiveUtils {

    public static int POINTS_PER_SECOND = 8;

    // Length of day in hours
    public static int DAY_LENGTH_HOURS = 9;

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
        //genome.addChromosome (new AffectiveChromosome(new Date()));
        return genome;
    }

    public static void doAffectiveHacks(List<Track> newTracks) {
         for(Track track : newTracks) {
             String name = track.getName();
             if(name.endsWith("-axis")) {
                 track.setDataRange(new DataRange(-1.5f, 0, 1.5f));
             }
             else if(name.equals("Battery")) {
                 track.setDataRange(new DataRange(-1, 1));
             }
             else if(name.endsWith("Celsius")) {
                 track.setDataRange(new DataRange(20, 30));
                 track.setRendererClass(LineplotRenderer.class);
             }
             else if(name.startsWith("EDA")) {
                 track.setDataRange(new DataRange(0, 10));
                 track.setColor(new Color(0, 150, 0));
             }

         }
    }
}
