/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.tools.sort;

import htsjdk.tribble.readers.AsciiLineReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author jrobinso
 */
public class InteractionSorter extends AsciiSorter {
    private int chrCol;
    private int startCol;


    public InteractionSorter(File inputFile, File outputFile) {
        super(inputFile, outputFile);
        chrCol = 1;
        startCol = 3;

    }

    @Override
    Parser getParser() {
        return new Parser(chrCol, startCol);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException {

        String nextLine = reader.readLine();
        writer.println(nextLine);
        return nextLine;


    }
}
