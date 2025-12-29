package org.broad.igv.tools.sort;

import org.broad.igv.logging.*;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jun 28, 2010
 * Time: 2:46:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class VCFSorter extends AsciiSorter {

    private static Logger log = LogManager.getLogger(VCFSorter.class);


    public VCFSorter(File inputFile, File outputFile) {
        super(inputFile, outputFile);
    }

    @Override
    Parser getParser() {
        return new Parser(0, 1);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) {
        try {
            String nextLine = reader.readLine();
            // TODO -- check "browser" line syntax,  is it a multi-line directive?
            while (nextLine.startsWith("#")) {
                writer.println(nextLine);
                nextLine = reader.readLine();
            }

            return nextLine;
        } catch (IOException e) {
            log.error("Error writing header", e);
            return null;
        }


    }
}
