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

import org.broad.igv.data.DataSource;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.PointsRenderer;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.DataSourceTrack;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TrackType;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
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
    public static final GenomeListItem GENOME_DESCRIPTOR = new GenomeListItem("Affective", "", "affective");
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


    public static void loadTDFFile(ResourceLocator locator, List<Track> newTracks, Genome genome, TDFReader reader) {

        TrackType type = reader.getTrackType();
        TrackProperties props = null;
        String trackLine = reader.getTrackLine();
        if (trackLine != null && trackLine.length() > 0) {
            props = new TrackProperties();
            ParsingUtils.parseTrackLine(trackLine, props);
        }

        // In case of conflict between the resource locator display name and the track properties name,
        // use the resource locator
        String name = locator.getName();
        if (name != null && props != null) {
            props.setName(name);
        }

        if (name == null) {
            name = props == null ? locator.getTrackName() : props.getName();
        }

        int trackNumber = 0;
        String path = locator.getPath();
        boolean multiTrack = reader.getTrackNames().length > 1;

        for (String heading : reader.getTrackNames()) {

            String trackId = multiTrack ? path + "_" + heading : path;
            String trackName = multiTrack ? heading : name;
            final DataSource dataSource = new AffectiveDataSource(reader, trackNumber, heading, genome);
            DataSourceTrack track = new DataSourceTrack(locator, trackId, trackName, dataSource);

            String displayName = (name == null || multiTrack) ? heading : name;
            track.setName(displayName);
            track.setTrackType(type);
            if (props != null) {
                track.setProperties(props);
            }

            if (trackName.endsWith("-axis")) {
                track.setDataRange(new DataRange(0, 1.5f));
            } else if (trackName.equals("Battery")) {
                track.setDataRange(new DataRange(-1, 0, 1));
                track.setRendererClass(PointsRenderer.class);
            } else if (trackName.endsWith("Celsius")) {
                track.setDataRange(new DataRange(20, 30));
                track.setRendererClass(PointsRenderer.class);
            } else if (trackName.startsWith("EDA")) {
                track.setDataRange(new DataRange(0, 10));
                track.setColor(new Color(0, 150, 0));
            }


            newTracks.add(track);
            trackNumber++;
        }


    }
}
