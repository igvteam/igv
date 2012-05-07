package org.broad.igv.tools.sort;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.tribble.readers.AsciiLineReader;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 4/12/12
 */
public class GWASSorter extends Sorter {

    static Logger log = Logger.getLogger(GWASSorter.class);

    private int chrCol = -1;
    private int startCol = -1;
    List<String> headerLines = new ArrayList<String>();


    public GWASSorter(File inputFile, File outputFile) {
        super(inputFile, outputFile);
        findColumns(inputFile);
    }

    private void findColumns(File inputFile) {
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(inputFile));
            String nextLine;
            while ((nextLine = br.readLine()) != null) {
                headerLines.add(nextLine);
                if (!nextLine.startsWith("#")) {
                    String[] tokens = Globals.singleTabMultiSpacePattern.split(nextLine);
                    for (int i = 0; i < tokens.length; i++) {
                        if (tokens[i].equalsIgnoreCase("chr")) {
                            chrCol = i;
                        } else if (tokens[i].equalsIgnoreCase("bp")) {
                            startCol = i;
                        }
                    }
                    break;
                }
            }
            if (chrCol < 0) {
                throw new RuntimeException("Could not find chromosome column");
            }
            if (startCol < 0) {
                throw new RuntimeException("Could not find start column");
            }
        } catch (IOException e) {
            log.error("Error reading GWAS file", e);
            throw new RuntimeException("Error reading GWAS file" + e.getMessage(), e);
        } finally {
            if (br != null) try {
                br.close();
            } catch (IOException e) {

            }
        }
    }


    @Override
    Parser getParser() {
        return new Parser(chrCol, startCol, true);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException {

        String nextLine;
        while((nextLine = reader.readLine()) != null) {
            writer.println(nextLine);
            if (!nextLine.startsWith("#")) {
                break;
            }
        }
        return null;

    }
}
