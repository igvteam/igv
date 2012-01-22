package org.broad.igv.dev.affective;

import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

/**
 * Special parser for time-series experiment.
 * <p/>
 * Head of sample file follows
 * Log File Created by Q V1.0 - (c) 2010 Affectiva Inc.
 * File Version: 1.01
 * Firmware Version: 1.50
 * UUID: AFQ441000BU
 * Sampling Rate: 8
 * Start Time: 2011-03-22 08:52:57 Offset:-04
 * Z-axis | Y-axis | X-axis | Battery | °Celsius | EDA(uS)
 * ---------------------------------------------------------
 * -0.050,-0.100,-1.150,-1,27.400,0.001
 *
 * @author Jim Robinson
 * @date 11/25/11
 *
 */
public class TSParser {

    DataConsumer dataConsumer;
    File root;


    public TSParser(String file, DataConsumer dataConsumer) {
        this.root = new File(file);
        this.dataConsumer = dataConsumer;
    }


    public void parse() throws IOException {

       for (File dir : root.listFiles()) {
            final File[] files = dir.listFiles();
            if (dir.isDirectory()) {
                String chrName = dir.getName().replace("Day ", "");
                for (File f : files) {
                    if (f.getName().endsWith(".csv")) {
                        parseFile(f, chrName);
                        break;
                    }
                }
            }
        }

        dataConsumer.parsingComplete();
    }

    private void parseFile(File file, String chrName) throws IOException {
        BufferedReader br = null;

        try {
            br = ParsingUtils.openBufferedReader(file.getAbsolutePath());

            for (int i = 0; i < 6; i++) br.readLine();

            String headerLine = br.readLine();
            String[] headers = headerLine.split("\\|");
            for (int i = 0; i < headers.length; i++) {
                headers[i] = headers[i].trim();
            }
            dataConsumer.setTrackParameters(TrackType.AFFECTIVE, null, headers);

            br.readLine(); // Skip ----------

            String nextLine;
            int step = 0;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("...")) {
                    System.out.println("Skipping " + nextLine);

                } else {
                    try {
                        String[] tokens = nextLine.split(",");
                        if (tokens.length == headers.length) {
                            float[] values = new float[tokens.length];
                            for (int i = 0; i < tokens.length; i++) {
                                values[i] = Float.parseFloat(tokens[i]);
                            }
                            dataConsumer.addData(chrName, step, step + 1, values, null);
                        }
                    } catch (Exception e) {
                        System.out.println(nextLine);
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }
                }
                step++;

            }
        } finally {

            if (br != null) br.close();
        }
    }


}
