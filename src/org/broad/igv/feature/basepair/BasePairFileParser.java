package org.broad.igv.feature.basepair;

import org.apache.batik.dom.svg12.Global;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.ParserException;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.basepair.BasePairData;
import org.broad.igv.track.BasePairTrack;
import org.broad.igv.renderer.BasePairRenderer;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import htsjdk.tribble.readers.AsciiLineReader;

import java.awt.*;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BasePairFileParser {

    static Logger log = Logger.getLogger(BasePairFileParser.class);

    BasePairData basePairData;
    BasePairTrack track;

    public BasePairTrack loadTrack(ResourceLocator locator, Genome genome) {
        AsciiLineReader reader = null;

        List<Object> rows = new ArrayList();
        List<Color> colors = new ArrayList();
        HashMap<Color, List<Object>> rowsByColor = new HashMap<Color, List<Object>>();

        int parseColumn = -1;
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

                int startLeftNuc = Integer.parseInt(tokens[0])-1; // stick to IGV's 0-based coordinate convention
                int startRightNuc = Integer.parseInt(tokens[1])-1;
                int endLeftNuc = Integer.parseInt(tokens[2])-1;
                int endRightNuc = Integer.parseInt(tokens[3])-1;
                Color color = colors.get(Integer.parseInt(tokens[4]));

                List<Object> columns = new ArrayList();
                if (startLeftNuc <= endRightNuc){
                    columns.add(Math.min(startLeftNuc, startRightNuc));
                    columns.add(Math.max(startLeftNuc, startRightNuc));
                    columns.add(Math.min(endLeftNuc, endRightNuc));
                    columns.add(Math.max(endLeftNuc, endRightNuc));
                } else {
                    // maintain left-to-right basepair order even if swapped in file
                    columns.add(Math.min(endLeftNuc, endRightNuc));
                    columns.add(Math.max(endLeftNuc, endRightNuc));
                    columns.add(Math.min(startLeftNuc, startRightNuc));
                    columns.add(Math.max(startLeftNuc, startRightNuc));
                }
                rowsByColor.get(color).add(columns);
                nextLine = reader.readLine();
                rowCounter++;
            }

            basePairData = new BasePairData(colors, rowsByColor);

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