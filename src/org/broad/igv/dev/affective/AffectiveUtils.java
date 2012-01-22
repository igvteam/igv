package org.broad.igv.dev.affective;

import java.io.*;

/**
 * @author Jim Robinson
 * @date 1/11/12
 */
public class AffectiveUtils {




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
}
