package org.igv.tools.sort;

import org.igv.Globals;
import org.igv.logging.*;
import org.igv.gwas.GWASParser;
import htsjdk.tribble.readers.AsciiLineReader;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * @author Jim Robinson
 * @date 4/12/12
 */
public class GWASSorter extends AsciiSorter {

    static Logger log = LogManager.getLogger(GWASSorter.class);

    List<String> headerLines = new ArrayList<String>();
    private GWASParser.GWASColumns columns;


    public GWASSorter(File inputFile, File outputFile) {
        super(inputFile, outputFile);
        this.columns = new GWASParser.GWASColumns();
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
                    Pattern delimiter = nextLine.indexOf('\t') > 0 ? Globals.tabPattern : Globals.whitespacePattern;
                    this.columns.parseHeader(nextLine, delimiter);
                    break;
                }
            }
            if (this.columns.chrCol < 0) {
                throw new RuntimeException("Could not find chromosome column");
            }
            if (this.columns.locationCol < 0) {
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
        return new Parser(this.columns.chrCol, this.columns.locationCol, true);
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
