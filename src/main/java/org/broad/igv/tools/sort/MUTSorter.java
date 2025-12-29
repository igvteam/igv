package org.broad.igv.tools.sort;

import org.broad.igv.feature.tribble.MUTCodec;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author jrobinso
 *         Date: 4/8/13
 *         Time: 4:23 PM
 */
public class MUTSorter extends AsciiSorter {

    int chrCol;
    int startCol;


    public MUTSorter(File inputFile, File outputFile) {

        super(inputFile, outputFile);
        MUTCodec codec = new MUTCodec(inputFile.getAbsolutePath(), null);
        chrCol = codec.getChrColumn();
        startCol = codec.getStartColumn();
    }


    @Override
    Parser getParser() throws IOException {

        return new Parser(chrCol, startCol);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException {

        String nextLine;
        // TODO -- check "browser" line syntax,  is it a multi-line directive?
        while ((nextLine = reader.readLine()).startsWith("#")) {
            writer.println(nextLine);   // Comments, directives
        }
        writer.println(nextLine);  // Header

        // Contract is to return the first data line!
        nextLine = reader.readLine();

        return nextLine;


    }
}
