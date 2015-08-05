/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.synteny;

import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

public class SyntenyUtils {
    static final Pattern SPACE = Pattern.compile(" ");
    static final Pattern TAB = Pattern.compile("\t");
    static String usageString = "USAGE: java -jar synteny.jar <mapping> <inputFile.seg> <outputFile.seg>";
    static String mappingString = "Recognized mappings:  canFam2tohg18  mm9tohg18";

    public static void main(String[] args) throws IOException {
        if (args.length < 3) {
            System.out.println(usageString);
            System.exit(-1);
        }

        String mapping = args[0];
        String inputFile = args[1];
        String outputFile = args[2];

        if (!inputFile.endsWith(".seg")) {
            System.out.println("Input file type not supported (" + inputFile + ").  Currently only segmented files (.seg) are supported.");
            System.exit(-1);
        }
        if (!outputFile.endsWith(".seg")) {
            System.out.println("Output file type not supported (" + inputFile + ").  Currently only segmented files (.seg) are supported.");
            System.exit(-1);
        }

        String mappingFile = null;
        if (mapping.equals("canFam2tohg18")) {
            mappingFile = "/resources/canFam2/mapWithUn_ph_MA4_ML10K.map";
        } else if (mapping.equals("mm9tohg18")) {
            mappingFile = "/resources/mm9/mapWithUn_ph_MA4_ML10K.clean.map";
        } else {
            System.out.println("Unsupported mapping: " + mapping);
            System.out.println(mappingString);
            System.exit(-1);
        }

        //InputStream is = SyntenyUtils.class.getResourceAsStream(mappingFile);
        mapSegments(mapping, inputFile, outputFile);
        //is.close();
    }

    public static void mapCNFile(String path, String cnFile, String oFile) {
        Map mappings = loadMappings(path, true);
        System.out.println("Mappings loaded");
        int chrCol = 1;
        int posCol = 2;
        try {
            BufferedReader br = new BufferedReader(new FileReader(cnFile));
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(oFile)));
            String nextLine = br.readLine();
            pw.println(nextLine);
            String lastChr = "";
            while ((nextLine = br.readLine()) != null) {
                String[] tokens = TAB.split(nextLine);
                String chr = tokens[chrCol];
                int position = Integer.parseInt(tokens[posCol]);

                if (!lastChr.equals(chr)) {
                    lastChr = chr;
                }

                List mList = (List) mappings.get(chr);

                if (mList != null) {
                    Region region = (Region) getMappingContaining(mList, position);
                    if (region != null) {
                        int toPos = (int) region.mapPosition(position);
                        if (toPos > 0) {
                            tokens[chrCol] = region.getToChr();
                            tokens[posCol] = String.valueOf(toPos);
                            int nTokens = tokens.length;
                            pw.print(tokens[0]);
                            for (int i = 1; i < nTokens; i++) {
                                pw.print("\t");
                                pw.print(tokens[i]);
                            }
                            pw.println();
                        }
                    }
                }
            }

            br.close();
            pw.close();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
        }
    }


    /**
     * Map segments from a ".seg" file
     *
     * @param mappingFile
     * @param iFile
     * @param oFile
     */
    private static void mapSegments(String  mappingFile, String iFile, String oFile) {
        Map mappings = loadMappings(mappingFile, true);

        int sampleColumn = 0;
        int chrColumn = 1;
        int startColumn = 2;
        int endColumn = 3;

        BufferedReader reader = null;
        PrintWriter pw = null;
        String nextLine = null;
        try {
            reader = new BufferedReader(new FileReader(iFile));
            pw = new PrintWriter(new BufferedWriter(new FileWriter(oFile)));


            nextLine = reader.readLine();
            while ((nextLine.startsWith("#")) || (nextLine.trim().length() == 0)) {
                pw.println(nextLine);
                nextLine = reader.readLine();
            }
            pw.println(nextLine);

            String sample;
            String chr;
            int start;
            int end;
            while (((nextLine = reader.readLine()) != null) && (nextLine.trim().length() > 0)) {

                String[] tokens = TAB.split(nextLine);
                int nTokens = tokens.length;
                if (nextLine.startsWith("#") || nTokens <= 4) {
                    System.out.println(nextLine);
                    continue;
                } else {
                    sample = tokens[sampleColumn];
                    chr = tokens[chrColumn];
                    start = Integer.parseInt(tokens[startColumn].trim());
                    end = Integer.parseInt(tokens[endColumn].trim());

                    String mappingChr = chr.startsWith("chr") ? chr : "chr" + chr;
                    List tmp = (List) mappings.get(mappingChr);
                    if (tmp == null) {
                        System.out.println("No mappings for chr: " + chr);
                        continue;
                    }

                    // get all rows overlapping the start-end interval.
                    // region R:chr18:chr7:D9 chr18 40411303 43748821 + chr7 46638925 49478088 -
                    // anchor A:chr2:39250:chr10:44583 chr2 63071936 63074348 + chr10 65929165 65931581 + 2469.0
                    List<Mapping> overlappingMappings = getMappingsOverlapping(tmp, start, end);
                    if (overlappingMappings == null) {
                        System.out.println("No mapping for: " + chr + ":" + start + "-" + end);
                    } else
                        for (Mapping mapping : overlappingMappings) {

                            int adjustedStart = Math.max(start, mapping.getFromStart());
                            int adjustedEnd = Math.min(end, mapping.getFromEnd() - 1);
                            int p1 = (int) mapping.mapPosition(adjustedStart);
                            int p2 = (int) mapping.mapPosition(adjustedEnd);
                            if ((p1 < 0) || (p2 < 0)) {
                                System.out.println("Unmapped position: " + chr + ":" + start + "-" + end + " -> (" + p1 + "  " + p2 + ")");
                            }

                            if (p1 == p2) {
                                System.out.println("Warning: start == end in mapped position " + mapping.toString());
                            }

                            int toStart = Math.min(p1, p2);
                            int toEnd = Math.max(p1, p2);
                            pw.print(sample + "\t" + mapping.getToChr() + "\t" + toStart + "\t" + toEnd);
                            for (int i = endColumn + 1; i < tokens.length; i++) {
                                pw.print("\t" + tokens[i]);
                            }
                            pw.println();
                        }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                if (reader == null) {
                    reader.close();
                }
                if (pw != null)
                    pw.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public static void toBed(String mappingFile, String ofile) {
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new BufferedWriter(new FileWriter(ofile)));

            Map<String, List<Mapping>> mappings = loadMappings(mappingFile, false);
            for (String chr : mappings.keySet())
                for (Mapping m : mappings.get(chr))
                    pw.println(((Region) m).toBed());
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            pw.close();
        }
    }

    public static Map<String, List<Mapping>> loadMappings(String path, boolean reverse) {
        BufferedReader reader = null;
        Map<String, List<Mapping>> mappings = new HashMap();
        Region currentRegion = null;
        try {
            InputStream stream = ParsingUtils.openInputStream(path);
            reader = new BufferedReader(new InputStreamReader(stream));

            Pattern.compile("\t");
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                if ((!nextLine.startsWith("region")) && (!nextLine.startsWith("anchor"))) {
                    continue;
                }
                String[] tokens = SPACE.split(nextLine);
                String type = tokens[0];

                Class c = type.equals("region") ? Region.class : Anchor.class;

                String name = tokens[1];
                String fromChr = tokens[2];
                int fromStart = Integer.parseInt(tokens[3]);
                int fromEnd = Integer.parseInt(tokens[4]);
                String fromStrand = tokens[5];

                if (!fromStrand.equals("+")) {
                    System.out.println("Negative from strand");  // <= unexpected!
                }

                String toChr = tokens[6];
                int toStart = Integer.parseInt(tokens[7]);
                int toEnd = Integer.parseInt(tokens[8]);
                String toStrand = tokens[9];

                AbstractMapping syntenyMapping = (AbstractMapping) c.newInstance();
                if (reverse) {
                    syntenyMapping.setParameters(name, toChr, toStart, toEnd, toStrand, fromChr, fromStart, fromEnd, fromStrand);
                } else {
                    syntenyMapping.setParameters(name, fromChr, fromStart, fromEnd, fromStrand, toChr, toStart, toEnd, toStrand);
                }

                if (c == Region.class) {
                    currentRegion = (Region) syntenyMapping;
                    List syntenyMappingList = (List) mappings.get(syntenyMapping.getFromChr());
                    if (syntenyMappingList == null) {
                        syntenyMappingList = new ArrayList(1000);
                        mappings.put(syntenyMapping.getFromChr(), syntenyMappingList);
                    }
                    syntenyMappingList.add(syntenyMapping);
                } else {
                    currentRegion.addAnchor((Anchor) syntenyMapping);
                }

            }

            for (List<Mapping> syntenyMappingList : mappings.values()) {
                sortMappingList(syntenyMappingList);
            }

            return mappings;
        } catch (Exception exception) {
            exception.printStackTrace();
            return null;
        } finally {
            if (reader != null)
                try {
                    reader.close();
                } catch (IOException iOException) {
                }
        }
    }

    public static Mapping getMappingContaining(List<Mapping> mappings, int fromPosition) {
        int idx = getIndexBefore(mappings, fromPosition);
        for (int i = idx; i < mappings.size(); i++) {
            Mapping mapping = mappings.get(i);
            if (mapping.containsFromPosition(fromPosition)) {
                return mapping;
            }
        }
        return null;
    }

    public static List<Mapping> getMappingsOverlapping(List<Mapping> mappings, int fromStart, int fromEnd) {
        ArrayList overlaps = new ArrayList();
        int idx = getIndexBefore(mappings, fromStart);

        for (int i = idx; i < mappings.size(); i++) {
            AbstractMapping m = (AbstractMapping) mappings.get(i);

            if ((m.getFromEnd() >= fromStart) && (m.getFromStart() <= fromEnd))
                overlaps.add(m);
            else {
                if (m.getFromStart() > fromEnd) {
                    break;
                }
            }
        }
        return overlaps;
    }

    private static void sortMappingList(List<Mapping> features) {
        Collections.sort(features, new Comparator<Mapping>() {
            public int compare(Mapping o1, Mapping o2) {
                return o1.getFromStart() - o2.getFromStart();
            }
        });
    }

    private static int getIndexBefore(List<Mapping> values, int x) {
        return getIndexBefore(values, x, 0, values.size());
    }

    private static int getIndexBefore(List<Mapping> values, int x, int leftBound, int rightBound) {
        int idx = (leftBound + rightBound) / 2;

        if ((idx == 0) || (idx == values.size() - 1)) {
            return idx;
        }
        if (((AbstractMapping) values.get(idx)).getFromStart() == x) {
            return idx;
        }

        if (((AbstractMapping) values.get(idx)).getFromStart() < x) {
            if (((AbstractMapping) values.get(idx + 1)).getFromStart() >= x) {
                return idx;
            }
            leftBound = idx;
            return getIndexBefore(values, x, leftBound, rightBound);
        }

        if (((AbstractMapping) values.get(idx - 1)).getFromStart() <= x) {
            return idx - 1;
        }
        rightBound = idx;
        return getIndexBefore(values, x, leftBound, rightBound);
    }
}