package org.broad.igv.cursor;

import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.AsciiFeatureCodec;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 1/14/14
 *         Time: 4:35 PM
 */
public class CursorUtils {


    public static List<CursorRegion> createRegions(CursorTrack track) {

        List<CursorRegion> frames = new ArrayList();
        Map<String, List<BasicFeature>> featureMap = track.getFeatureMap();

        for (List<BasicFeature> flist : featureMap.values()) {

            for (BasicFeature feature : flist) {

                int locaction = (feature.getStart() + feature.getEnd()) / 2;
                frames.add(new CursorRegion(feature.getChr(), locaction));
            }
        }

        return frames;

    }


    public static CursorTrack loadTrack(String peakFile) throws IOException {

        AsciiFeatureCodec codec = (AsciiFeatureCodec) CodecFactory.getCodec(new ResourceLocator(peakFile), null);
//        if(codec == null) {
//            // TODO -- inform user
//            System.out.println("Skipping " + peakFile);
//            return null;
//        }

//        EncodePeakCodec codec = new EncodePeakCodec();

        Map<String, List<BasicFeature>> featureMap = new HashMap<String, List<BasicFeature>>();
        BufferedReader br = null;
        TrackProperties props = null;
        try {
            br = ParsingUtils.openBufferedReader(peakFile);

            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                if (nextLine.startsWith("track")) {
                    props = new TrackProperties();
                    ParsingUtils.parseTrackLine(nextLine, props);
                } else {
                    BasicFeature f = (BasicFeature) codec.decode(nextLine);
                    if (f != null) {
                        String chr = f.getChr();
                        List<BasicFeature> featureList = featureMap.get(chr);
                        if (featureList == null) {
                            featureList = new ArrayList<BasicFeature>();
                            featureMap.put(chr, featureList);
                        }
                        featureList.add(f);
                    }
                }
            }
        } finally {
            if (br != null) br.close();
        }


        for (List<BasicFeature> featureList : featureMap.values()) {
            Collections.sort(featureList, new Comparator<BasicFeature>() {
                @Override
                public int compare(BasicFeature o1, BasicFeature o2) {
                    return o1.getStart() - o2.getStart();
                }
            });

        }

        CursorTrack track = new CursorTrack(featureMap, codec.getFeatureType());
        String trackName = (new File(peakFile)).getName();
        Color trackColor = null;
        if (props != null) {
            trackColor = props.getColor();
            trackName = props.getName();
        }
        if (trackColor == null && trackName != null) {
            if (trackName.toLowerCase().contains("k4me3")) {
                trackColor = new Color(0, 150, 0);
            } else if (trackName.toLowerCase().contains("k27me3")) {
                trackColor = new Color(150, 0, 0);
            }
        }

        if(trackName != null) track.setName(trackName);
        if(trackColor != null) track.setColor(trackColor);

        return track;
    }


}
