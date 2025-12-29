package org.igv.variant.util;

import java.io.*;

/**
 * @author jrobinso
 * @date Apr 25, 2011
 */
public class AutismUtils {


    public static void main(String[] args) throws IOException {
        createSampleMap();
    }

    private static void createLoadXML() throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader("/Users/jrobinso/IGV/autism/bam_locations.txt"));

        PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("/Users/jrobinso/IGV/autism/bam.xml")));

        String nextLine;
        while ((nextLine = reader.readLine()) != null) {

            String bamPath = nextLine.trim();
            String url = bamPath.replace("/seq", "http://iwww.broadinstitute.org/igvdata");
            String name = (new File(bamPath)).getName().replace(".bam", "");

            writer.println("   <Resource name=\"" + name + "\"");
            writer.println("             path=\"" + url + "\"/>");

        }

        reader.close();
        writer.close();
    }

    private static void createSampleMap() throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader("/Users/jrobinso/IGV/autism/bam_locations.txt"));

        PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("/Users/jrobinso/IGV/autism/samplemap.txt")));

        String nextLine;
        while ((nextLine = reader.readLine()) != null) {

            String bamPath = nextLine.trim();
            String url = bamPath.replace("/seq", "http://iwww.broadinstitute.org/igvdata");
            String name = (new File(bamPath)).getName().replace(".bam", "");

            writer.println(name + "\t" + url);

        }

        reader.close();
        writer.close();
    }
}
