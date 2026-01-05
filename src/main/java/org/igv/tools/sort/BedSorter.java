package org.igv.tools.sort;

import htsjdk.tribble.readers.AsciiLineReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author jrobinso
 */
public class BedSorter extends AsciiSorter {
    private int chrCol;
    private int startCol;


    public BedSorter(File inputFile, File outputFile) {
        super(inputFile, outputFile);
        if (inputFile.getName().endsWith(".psl") || inputFile.getName().endsWith(".pslx")) {
            chrCol = 13;
            startCol = 15;
        } else {
            chrCol = 0;
            startCol = 1;
        }
    }

    @Override
    Parser getParser() {
        return new Parser(chrCol, startCol);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException {

        String nextLine = reader.readLine();

        if (nextLine.startsWith("psLayout")) {
            do {
                writer.println(nextLine);
                nextLine = reader.readLine();
            } while (!nextLine.startsWith("-"));
            nextLine = reader.readLine();
        }
        // TODO -- check "browser" line syntax,  is it a multi-line directive?
        while (nextLine.startsWith("#") ||
                nextLine.startsWith("browser") ||
                nextLine.startsWith("track") ||
                nextLine.trim().length() == 0) {
            writer.println(nextLine);
            nextLine = reader.readLine();
        }

        return nextLine;


    }
}
