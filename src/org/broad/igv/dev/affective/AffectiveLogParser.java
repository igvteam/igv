/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.dev.affective;

import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.tools.parsers.UnsortedException;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.Arrays;
import java.util.Comparator;
import java.util.regex.Pattern;

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
 * Z-axis | Y-axis | X-axis | Battery | Celsius | EDA(uS)
 * ---------------------------------------------------------
 * -0.050,-0.100,-1.150,-1,27.400,0.001
 *
 * @author Jim Robinson
 * @date 11/25/11
 */
public class AffectiveLogParser {


    static int CHR_SAMPLING_RATE = 8; // 8 points per second
    static int CHR_START_TIME = (AffectiveUtils.START_TIME_HR * 60) * 60;

    DataConsumer dataConsumer;
    File root;
    Pattern commaPattern = Pattern.compile(",");

    public AffectiveLogParser(String file, DataConsumer dataConsumer) {
        this.root = new File(file);
        this.dataConsumer = dataConsumer;
    }


    public void parse() throws IOException {

        if (root.isDirectory()) {
            final File[] files = root.listFiles();
            Arrays.sort(files, new Comparator<File>() {
                public int compare(File file, File file1) {
                    return file.getName().compareTo(file1.getName());
                }
            });
            for (File f : files) {
                if (f.getName().endsWith(".csv")) {
                    System.out.println(f);
                    parseFile(f);
                }
            }
        } else {
            parseFile(root);
        }

        dataConsumer.parsingComplete();
    }

    private void parseFile(File file) throws IOException {
        BufferedReader br = null;

        try {
            br = new BufferedReader(new InputStreamReader(
                    ParsingUtils.openInputStream(file.getAbsolutePath()), Charset.forName("UTF-8")));

            Header header = parseHeader(br);
            String chrName = header.getStartDate();
            final String[] trackNames = header.getTrackNames();
            dataConsumer.setTrackParameters(TrackType.AFFECTIVE, null, trackNames, false);

            // Set file/data/chromosome specific attributes
            String prefix = "ATTR:" + chrName + ":";
            dataConsumer.setAttribute(prefix + "uuid", header.getUuid());
            dataConsumer.setAttribute(prefix + "samplingRate", String.valueOf(header.getSamplingRate()));
            dataConsumer.setAttribute(prefix + "startTime", String.valueOf(CHR_START_TIME)); //    header.getStartTime()));

            int stepSize = (int) (Math.round((double) CHR_SAMPLING_RATE) / header.getSamplingRate());
            int startOffset = (header.getStartTime() - CHR_START_TIME) * CHR_SAMPLING_RATE;

            String nextLine;

            int startTime = startOffset;
            while ((nextLine = br.readLine()) != null) {

                try {
                    String[] tokens = commaPattern.split(nextLine);
                    if (tokens.length == trackNames.length) {
                        float[] values = new float[tokens.length];
                        for (int i = 0; i < tokens.length; i++) {
                            float v;
                            try {
                                v = Float.parseFloat(tokens[i]);
                                if (i < 3) v = Math.abs(v);
                            } catch (NumberFormatException e) {
                                v = Float.NaN;
                            }
                            values[i] = v;
                        }

                        try {
                            dataConsumer.addData(chrName, startTime, startTime + stepSize, values, null);
                        } catch (UnsortedException e) {
                            // Ignore
                        }
                    } else {
                        System.out.println("Skipping line: " + nextLine);
                    }
                } catch (Exception e) {
                    System.out.println(nextLine);
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }

                startTime += stepSize;

            }
        } finally {

            if (br != null) br.close();
        }
    }


    /**
     * Parse the header fields, and leave the reader positioned at the start of the data section.
     * <p/>
     * Assumption: reader is positioned at the beginning of the file
     * <p/>
     * This is not private so it can be unit tested.
     *
     * @param reader
     * @return
     */
    Header parseHeader(BufferedReader reader) throws IOException {

        String nextLine;
        Header header = new Header();
        while ((nextLine = reader.readLine()) != null) {
            if (nextLine.startsWith("UUID:")) {
                String[] tokens = nextLine.split("\\s+");
                if (tokens.length > 1) {
                    header.uuid = tokens[1];
                }
            } else if (nextLine.startsWith("Sampling Rate:")) {
                String[] tokens = nextLine.split("\\s+");
                if (tokens.length > 2) {
                    header.samplingRate = Integer.parseInt(tokens[2]);
                }

            } else if (nextLine.startsWith("Start Time:")) {
                String[] tokens = nextLine.split("\\s+");
                if (tokens.length > 2) {
                    header.startDate = tokens[2];  //(new SimpleDateFormat("yyyy-MM-dd"))
                    String timeString = tokens[3];
                    String[] hms = timeString.split(":");

                    int hour = Integer.parseInt(hms[0]);
                    int minute = Integer.parseInt(hms[1]);
                    int second = Integer.parseInt(hms[2]);
                    header.startTime = (hour * 60 + minute) * 60 + second;

                    String headerLine = reader.readLine();
                    String[] trackNames = headerLine.split("\\|");
                    for (int i = 0; i < trackNames.length; i++) {
                        trackNames[i] = trackNames[i].trim();
                    }
                    header.trackNames = trackNames;

                }

            } else if (nextLine.startsWith("------")) {
                break;
            }
        }
        return header;
    }


    /**
     * Class representing an Affectiva log file header.  Its assumed these headers are fixed length.
     * <p/>
     * Example:
     * Log File Created by Q V1.0 - (c) 2010 Affectiva Inc.
     * File Version: 1.01
     * Firmware Version: 1.50
     * UUID: AFQ441000BU
     * Sampling Rate: 8
     * Start Time: 2011-04-06 08:57:35 Offset:-04
     * Z-axis | Y-axis | X-axis | BAttery | Celsius | EDA(uS)
     * ---------------------------------------------------------
     */
    static class Header {
        private String uuid;
        private int samplingRate;
        private String startDate;
        private int startTime;
        private String[] trackNames;

        public String getUuid() {
            return uuid;
        }

        public String getStartDate() {
            return startDate;
        }

        public int getStartTime() {
            return startTime;
        }

        public String[] getTrackNames() {
            return trackNames;
        }

        public int getSamplingRate() {
            return samplingRate;
        }
    }

}
