package org.broad.igv.tools.sort;

import htsjdk.tribble.readers.AsciiLineReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * User: jrobinso
 * Date: Jun 28, 2010
 * Time: 2:37:24 PM
 */
public class GFFSorter extends AsciiSorter {
    private final int chrCol = 0;
    private final int startCol = 3;


    public GFFSorter(File inputFile, File outputFile) {
        super(inputFile, outputFile);
    }

    @Override
    Parser getParser() {
        return new Parser(chrCol, startCol);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) {
        String nextLine = null;

        try {
            nextLine = reader.readLine();
            while (nextLine.startsWith("#")) {
                writer.println(nextLine);
                nextLine = reader.readLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        return nextLine;
    }
}
