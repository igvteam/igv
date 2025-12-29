package org.igv.tools.sort;

import org.igv.Globals;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.IOException;

/**
 * @author: nazaire
 */
public class Parser {

    private int chrCol;
    private int startCol;
    boolean splitOnWhiteSpace;
    private String commentPrefix;

    public Parser(int chrCol, int startCol) {
        this(chrCol, startCol, false, "#");
    }

    public Parser(int chrCol, int startCol, boolean splitOnWhiteSpace) {
        this(chrCol, startCol, splitOnWhiteSpace, "#");
    }

    public Parser(int chrCol, int startCol, boolean splitOnWhitespace, String commentPrefix) {
        this.chrCol = chrCol;
        this.startCol = startCol;
        this.splitOnWhiteSpace = splitOnWhitespace;
        this.commentPrefix = commentPrefix;
    }


    public SortableRecord readNextRecord(AsciiLineReader reader) {
        String nextLine = null;
        try {
            nextLine = reader.readLine();
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        if (nextLine == null) {
            return null;
        } else if (nextLine.trim().length() == 0 || nextLine.startsWith(commentPrefix)) {
            return readNextRecord(reader);
        }

        try {
            return createRecord(nextLine);
        } catch (ArrayIndexOutOfBoundsException e) {
            e.printStackTrace();
            System.out.println(nextLine);
            throw e;
        }
    }

    public SortableRecord createRecord(String nextLine) {
        String[] fields = splitOnWhiteSpace ?
                Globals.singleTabMultiSpacePattern.split(nextLine) :
                Globals.tabPattern.split(nextLine, -1);

        String chr = fields[chrCol];

        int start;
        try {
            start = Integer.parseInt(fields[startCol].trim());
        } catch (NumberFormatException e) {
            start = Integer.MAX_VALUE;
        }
        String text = nextLine;

        return new SortableRecord(chr, start, text);
    }
}
