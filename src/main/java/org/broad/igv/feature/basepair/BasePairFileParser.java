package org.broad.igv.feature.basepair;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.StringUtils;
import htsjdk.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.util.*;

public class BasePairFileParser {

    static Logger log = Logger.getLogger(BasePairFileParser.class);

    public static void loadData(ResourceLocator locator,
                                Genome genome,
                                BasePairData basePairData,
                                BasePairTrack.RenderOptions renderOptions) {
        AsciiLineReader reader = null;

        String nextLine = null;
        int rowCounter = 0;

        java.util.List<String> colors = new ArrayList();
        java.util.List<String> colorLabels = new ArrayList();

        try {
            reader = ParsingUtils.openAsciiReader(locator);
            // read lines specifying arc colors
            nextLine = reader.readLine();
            rowCounter++;
            if (nextLine.substring(0, 6).equals("color:")) { // hopefully handle the empty-file case
                while (nextLine.substring(0, 6).equals("color:")) {
                    String[] tokens = Globals.whitespacePattern.split(nextLine, -1);
                    int r = Integer.parseInt(tokens[1]);
                    int g = Integer.parseInt(tokens[2]);
                    int b = Integer.parseInt(tokens[3]);
                    String label = "";
                    try {
                        label = StringUtils.join(Arrays.copyOfRange(tokens, 4, tokens.length), " ");
                    } catch (IndexOutOfBoundsException e) {
                    }
                    Color color = new Color(r, g, b, 255);
                    colors.add(ColorUtilities.colorToString(color));
                    colorLabels.add(label);
                    nextLine = reader.readLine();
                    rowCounter++;
                }

                while (nextLine != null) {

                    String[] tokens = Globals.whitespacePattern.split(nextLine, -1);

                    String chr = (genome == null ? tokens[0] : genome.getCanonicalChrName(tokens[0]));   // TODO Future use

                    int startLeftNuc = Integer.parseInt(tokens[1]) - 1; // stick to IGV's 0-based coordinate convention
                    int startRightNuc = Integer.parseInt(tokens[2]) - 1;
                    int endLeftNuc = Integer.parseInt(tokens[3]) - 1;
                    int endRightNuc = Integer.parseInt(tokens[4]) - 1;
                    int colorIndex = Integer.parseInt(tokens[5]);

                    BasePairFeature feature;
                    if (startLeftNuc <= endRightNuc) {
                        feature = new BasePairFeature(chr,
                                Math.min(startLeftNuc, startRightNuc),
                                Math.max(startLeftNuc, startRightNuc),
                                Math.min(endLeftNuc, endRightNuc),
                                Math.max(endLeftNuc, endRightNuc),
                                colorIndex);
                    } else {
                        feature = new BasePairFeature(chr,
                                Math.min(endLeftNuc, endRightNuc),
                                Math.max(endLeftNuc, endRightNuc),
                                Math.min(startLeftNuc, startRightNuc),
                                Math.max(startLeftNuc, startRightNuc),
                                colorIndex);
                    }

                    basePairData.addFeature(feature);

                    nextLine = reader.readLine();
                    rowCounter++;

                }
            }
            renderOptions.setColors(colors);
            renderOptions.setColorLabels(colorLabels);

            basePairData.finish();   // ensure features are sorted by start position, important for rendering optimization

        } catch (Exception e) {
            log.error("Error parsing base pair file", e);
            if (nextLine != null && rowCounter != 0) {
                throw new ParserException(e.getMessage(), e, rowCounter, nextLine);
            } else {
                throw new RuntimeException(e);
            }
        } finally {
            if (reader != null) {
                reader.close();
            }
        }
    }
}
