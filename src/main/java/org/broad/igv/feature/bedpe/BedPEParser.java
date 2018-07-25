package org.broad.igv.feature.bedpe;

import org.apache.log4j.Logger;
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

    public static List<BedPEFeature> parse(String file) throws IOException {

        int colorColumn = -1;
        int thicknessColumn = -1;

        Map<String, Color> colorCache = new HashMap<>();

        List<BedPEFeature> features = new ArrayList<>();

        BufferedReader br = null;

        br = ParsingUtils.openBufferedReader(file);

        String nextLine;
        while ((nextLine = br.readLine()) != null) {

            if (nextLine.startsWith("#")) {
                if (nextLine.startsWith("#columns")) {
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
                }

            } else {
                String[] tokens = ParsingUtils.WHITESPACE_PATTERN.split(nextLine);

                if (tokens.length < 6) {
                    log.info("Skipping line: " + nextLine);
                    continue;
                }

                BedPEFeature feature = new BedPEFeature();
                feature.chr1 = tokens[0];
                feature.start1 = Integer.parseInt(tokens[1]);
                feature.end1 = Integer.parseInt(tokens[2]);
                feature.chr2 = tokens[3];
                feature.start2 = Integer.parseInt(tokens[4]);
                feature.end2 = Integer.parseInt(tokens[5]);

                if (tokens.length > 6) {
                    feature.name = tokens[6];
                }

                if (tokens.length > 7) {
                    feature.score = tokens[7];
                }

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

                // Skipping remaining fields for now

                features.add(feature);
            }

        }

        return features;


    }


}
