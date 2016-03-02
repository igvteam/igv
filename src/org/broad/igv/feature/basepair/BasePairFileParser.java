package org.broad.igv.feature.basepair;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class BasePairFileParser {

    static Logger log = Logger.getLogger(BasePairFileParser.class);

    BasePairTrack track;

    public BasePairTrack loadTrack(ResourceLocator locator, Genome genome) {
        AsciiLineReader reader = null;

        List<Color> colors = new ArrayList();
        HashMap<Color, List<Object>> rowsByColor = new HashMap<Color, List<Object>>();
        BasePairData basePairData = new BasePairData();

        String nextLine = null;
        int rowCounter = 0;

        try {
            reader = ParsingUtils.openAsciiReader(locator);
            // read lines specifying arc colors
            nextLine = reader.readLine();
            rowCounter++;
            while (nextLine.substring(0,6).equals("color:")){
                String[] tokens = Globals.tabPattern.split(nextLine, -1);
                int r = Integer.parseInt(tokens[1]);
                int g = Integer.parseInt(tokens[2]);
                int b = Integer.parseInt(tokens[3]);
                Color color = new Color(r,g,b,255);
                colors.add(color);
                nextLine = reader.readLine();
                rowCounter++;
            }

            for (Color color : colors){
                rowsByColor.put(color, new ArrayList());
            }

            while (nextLine != null) {

                String[] tokens = Globals.tabPattern.split(nextLine, -1);

                int nTokens = tokens.length;

                String chr = (genome == null ? tokens[0] : genome.getCanonicalChrName(tokens[0]));   // TODO Future use

                int startLeftNuc = Integer.parseInt(tokens[1])-1; // stick to IGV's 0-based coordinate convention
                int startRightNuc = Integer.parseInt(tokens[2])-1;
                int endLeftNuc = Integer.parseInt(tokens[3])-1;
                int endRightNuc = Integer.parseInt(tokens[4])-1;
                Color color = colors.get(Integer.parseInt(tokens[5]));

                BasePairFeature feature;
                if (startLeftNuc <= endRightNuc){
                    feature = new BasePairFeature(chr,
                            Math.min(startLeftNuc, startRightNuc),
                            Math.max(startLeftNuc, startRightNuc),
                            Math.min(endLeftNuc, endRightNuc),
                            Math.max(endLeftNuc, endRightNuc),
                            color);
                } else {
                    feature = new BasePairFeature(chr,
                            Math.min(endLeftNuc, endRightNuc),
                            Math.max(endLeftNuc, endRightNuc),
                            Math.min(startLeftNuc, startRightNuc),
                            Math.max(startLeftNuc, startRightNuc),
                            color);
                }

                basePairData.addFeature(feature);

                nextLine = reader.readLine();
                rowCounter++;

            }

            basePairData.finish();   // Insure features are sorted by start position, important for rendering optimization

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

        //for (Object row : rows){
        //    System.out.println(row);
        //}

        track = new BasePairTrack(basePairData, locator.getTrackName());
        return track;
    }

}