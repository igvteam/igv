package org.broad.igv.feature.bedpe;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ParsingUtils;

import java.awt.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by jrobinso on 6/29/18.
 */
public class BedPEParser {

    private static Logger log = Logger.getLogger(BedPEParser.class);

    public static List<BedPEFeature> parse(String file, boolean isClusters, Genome genome) throws IOException {

        int colorColumn = -1;
        int thicknessColumn = -1;
        boolean tenx = false;
        boolean parsedHeader = true;
        String[] columns;
        boolean col7isNumeric = true;   // Until proven otherwise

        Map<String, Color> colorCache = new HashMap<>();
        List<BedPEFeature> features = new ArrayList<>();
        BufferedReader br = null;
        br = ParsingUtils.openBufferedReader(file);
        String nextLine;
        while ((nextLine = br.readLine()) != null) {

            if (nextLine.startsWith("#columns")) {
                // An IGV hack, not sure anyone is using this
                try {
                    String[] t1 = ParsingUtils.WHITESPACE_PATTERN.split(nextLine);
                    if (t1.length == 2) {
                        String[] t2 = ParsingUtils.SEMI_COLON_PATTERN.split(t1[1]);
                        for (String keyValue : t2) {
                            String[] t = keyValue.split("=");
                            if (t[0].equals("color")) {
                                colorColumn = Integer.parseInt(t[1]) - 1;
                            } else if (t[0].equals("thickness")) {
                                thicknessColumn = Integer.parseInt(t[1]) - 1;
                            }
                        }
                    }
                } catch (NumberFormatException e) {
                    log.error("Error parsing #column line.", e);
                }
            } else if (nextLine.startsWith("#")) {
                columns = Globals.tabPattern.split(nextLine);
                if (nextLine.trim().equals("#chrom1\tstart1\tstop1\tchrom2\tstart2\tstop2\tname\tqual\tstrand1\tstrand2\tfilters\tinfo")) {
                    tenx = true;
                } else {
                    for (int i = 6; i < columns.length; i++) {
                        if (columns[i].equalsIgnoreCase("color")) {
                            colorColumn = i;
                        } else if (columns[i].toLowerCase().equalsIgnoreCase("thickness")) {
                            thicknessColumn = i;
                        }
                    }
                }

            } else if (nextLine.startsWith("track") || nextLine.startsWith("##track")) {
                TrackProperties trackProperties = new TrackProperties();
                ParsingUtils.parseTrackLine(nextLine, trackProperties);
            } else {
                String[] tokens = Globals.tabPattern.split(nextLine);

                if (tokens.length < 6) {
                    log.info("Skipping line: " + nextLine);
                    continue;
                }

                String chr1 = genome == null ? tokens[0] : genome.getCanonicalChrName(tokens[0]);
                String chr2 = genome == null ? tokens[3] : genome.getCanonicalChrName(tokens[3]);
                int start1 = Integer.parseInt(tokens[1]);
                int end1 = Integer.parseInt(tokens[2]);
                int start2 = Integer.parseInt(tokens[4]);
                int end2 = Integer.parseInt(tokens[5]);

                BedPEFeature feature = new BedPEFeature(chr1, start1, end1, chr2, start2, end2);

                if (tokens.length > 6) {
                    if (isClusters) {
                        feature.score = Double.parseDouble(tokens[6]);
                    } else {
                        feature.name = tokens[6];
                        col7isNumeric = col7isNumeric && isNumeric(tokens[6]);
                    }
                } else {
                    col7isNumeric = false;
                }

                if (tokens.length > 7) {
                    feature.score = Double.parseDouble(tokens[7]);
                }

                if (tenx) {
                    Map<String, String> attributes = new HashMap<>();
                    if (!tokens[8].equals(".")) {
                        attributes.put("filters", tokens[8]);
                    }
                    String[] kvPairs = Globals.semicolonPattern.split(tokens[11]);
                    for (String kvPair : kvPairs) {
                        String[] kv = Globals.equalPattern.split(kvPair);
                        attributes.put(kv[0], kv[1]);
                    }
                    feature.attributes = attributes;
                    feature.type = attributes.get("TYPE");
                } else {
                    if (colorColumn > 0) {
                        String colorString = tokens[colorColumn];
                        Color c = colorCache.get(colorString);
                        if (c == null) {
                            c = ColorUtilities.stringToColor(colorString);
                            colorCache.put(colorString, c);
                        }
                        feature.color = c;
                    }

                    if (thicknessColumn > 0) {
                        feature.thickness = Integer.parseInt(tokens[thicknessColumn]);
                    }
                }
                // Skipping remaining fields for now

                features.add(feature);
            }
        }


        // A hack to detect "interaction" bedpe files, which are not spec compliant.  Interaction score is column 7
        if (col7isNumeric) {
            for (BedPEFeature f : features) {
                f.score = Double.parseDouble(f.name);
                f.name = null;
            }
        }

        return features;


    }


    public static boolean isNumeric(String strNum) {
        return strNum.matches("-?\\d+(\\.\\d+)?");
    }

}
