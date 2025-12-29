package org.broad.igv.feature;

//~--- non-JDK imports --------------------------------------------------------

import org.broad.igv.logging.*;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Class description
 *
 * @author Enter your name here...
 * @version Enter version here..., 08/10/16
 */
public class CytoBandFileParser {

    private static Logger logger = LogManager.getLogger(CytoBandFileParser.class);


    /**
     * Validate cytoband file
     *
     * @param reader
     * @param filename
     * @return
     */
    public static boolean isValid(BufferedReader reader, String filename) {

        if (reader == null) {
            return false;
        }

        try {
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                String[] tokens = nextLine.split("\t");
                String chr = tokens[0].trim();
                Cytoband cytoData = new Cytoband(chr);
                parseData(tokens, cytoData);
            }
            return true;

        } catch (Exception e) {
            return false;
        }
    }

    /**
     * Method description
     *
     * @param reader
     * @return
     */
    public static LinkedHashMap<String, List<Cytoband>> loadData(BufferedReader reader) {

        LinkedHashMap<String, List<Cytoband>> dataMap = new LinkedHashMap<String, List<Cytoband>>();
        try {

            String nextLine;
            while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
                String[] data = nextLine.split("\t");
                String chr = data[0].trim();
                List<Cytoband> cytobands = dataMap.get(chr);
                if (cytobands == null) {
                    cytobands = new ArrayList<Cytoband>();
                    dataMap.put(chr, cytobands);
                }
                Cytoband cytoData = new Cytoband(chr);
                parseData(data, cytoData);
                cytobands.add(cytoData);
            }

            reader.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
        return dataMap;

    }

    private static void parseData(String[] tokens, Cytoband cytoData) {

        cytoData.setStart(Integer.parseInt(tokens[1].trim()));
        cytoData.setEnd(Integer.parseInt(tokens[2].trim()));
        if (tokens.length > 3) {
            cytoData.setName(tokens[3]);
        }
        if (tokens.length > 4) {
            if (tokens[4].equals("acen")) {
                cytoData.setType('c');
            } else {
                cytoData.setType(tokens[4].charAt(1));
                if (cytoData.getType() == 'p') {
                    String stainString = tokens[4].substring(4).trim();
                    short stain = stainString.length() == 0 ? 100 : Short.parseShort(stainString);
                    cytoData.setStain(stain);
                }
            }
        }

    }

    private static String parseChromosome(String[] data) {
        String chr = data[0].substring(3);
        int underscore = chr.indexOf('_');
        if (underscore > 0) {
            chr = chr.substring(0, underscore);
        }
        return chr;
    }
}
