/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.igv.tools.sort;

import htsjdk.tribble.readers.AsciiLineReader;
import org.igv.Globals;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Comparator;

/**
 * @author jrobinso
 */
public class SAMSorter extends AsciiSorter {

    public SAMSorter(File inputFile, File outputFile) {
        super(inputFile, outputFile);
    }

    @Override
    Parser getParser() {
        return new Parser(2, 3);
    }

    @Override
    String writeHeader(AsciiLineReader reader, PrintWriter writer) throws IOException {
        String nextLine = reader.readLine();
        while (nextLine != null && nextLine.startsWith("@")) {
            writer.println(nextLine);
            nextLine = reader.readLine();
        }

        // First alignment row
        return nextLine;
    }


    public static Comparator<SortableRecord> ReadNameComparator = new Comparator<SortableRecord>() {

        public int compare(SortableRecord o1, SortableRecord o2) {
            String[] t1 = Globals.tabPattern.split(o1.getText());
            String[] t2 = Globals.tabPattern.split(o2.getText());
            return t1[0].compareTo(t2[0]);

        }
    };
}
