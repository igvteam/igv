package org.igv.tools.sort;

import htsjdk.tribble.readers.AsciiLineReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Created with IntelliJ IDEA.
 * User: jrobinso
 * Date: 4/26/13
 * Time: 10:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class EQTLSorter extends AsciiSorter {

    int chrCol;
    int startCol;


    public EQTLSorter(File inputFile, File outputFile) {

        super(inputFile, outputFile);
        chrCol = 1;
        startCol =2;
    }


    @Override
    Parser getParser() throws IOException {

        return new Parser(chrCol, startCol);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException {

        // EQTL files have a single header line
        String headerLine = reader.readLine();
        writer.println(headerLine);  // Header

        // Contract is to return the first data line!
        return reader.readLine();

    }
}
