package org.broad.igv.feature.tribble;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.DataLoadException;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.util.ParsingUtils;

import java.io.IOException;

/**
 * Created by jrobinso on 10/31/17.
 * <p>
 * reference: https://github.com/lh3/miniasm/blob/master/PAF.md#paf-a-pairwise-mapping-format
 */
public class PAFCodec extends AsciiFeatureCodec<PAFFeature> {

    private static Logger log = Logger.getLogger(MUTCodec.class);

    private static int chrColumn = 5;
    private static int startColumn = 7;
    private static int endColumn = 8;
    private static int scoreColumn = 11;

    private String path;  // for error messages
    private Genome genome;
    private int errorCount = 0;


    public PAFCodec(String path, Genome genome) {
        super(PAFFeature.class);
        this.path = path;
        this.genome = genome;
        try {
            LineIterator reader = new LineIteratorImpl(new AsciiLineReader(ParsingUtils.openInputStream(path)));
            readActualHeader(reader);
        } catch (IOException e) {
            log.error(e.getMessage(), e);
        }
    }

    public Object readActualHeader(LineIterator reader) {
        return null;
    }


    //   public PAFFeature(String chr, int start, int end, Strand strand, float score, String name, MultiMap<String, String> attributes) {


    @Override
    public PAFFeature decode(String line) {

        String[] tokens = Globals.tabPattern.split(line);

        StringBuffer description = new StringBuffer();

        String chr = genome == null ? tokens[chrColumn].trim() : genome.getCanonicalChrName(tokens[chrColumn].trim());

        int start;
        try {
            start = Integer.parseInt(tokens[startColumn].trim());
        } catch (NumberFormatException e) {
            throw new DataLoadException("Column " + (startColumn + 1) + " must be a numeric value.", path);
        }

        int end;
        try {
            end = Integer.parseInt(tokens[endColumn].trim());
        } catch (NumberFormatException e) {
            throw new DataLoadException("Column " + (endColumn + 1) + " must be a numeric value.", path);
        }

        Strand strand;
        if (tokens[4].equals("+")) strand = Strand.POSITIVE;
        else if (tokens[4].equals("-")) strand = Strand.NEGATIVE;
        else strand = Strand.NONE;

        float score;
        try {
            score = (float) Double.parseDouble(tokens[scoreColumn]);
        } catch (NumberFormatException e) {
            throw new DataLoadException("Column " + (scoreColumn + 1) + " must be a numeric value.", path);
        }

        String name = tokens[0] + ":" + tokens[2] + "-" + tokens[3];


        description.append("Query sequence name: " + tokens[0] + "<br>");
        description.append("Query sequence length: " + tokens[1] + "<br>");
        description.append("Query start (0-based): " + tokens[2] + "<br>");
        description.append("Query end (0-based):" + tokens[3] + "<br>");
        description.append("Relative strand: " + tokens[4] + "<br>");
        description.append("Target sequence name: " + tokens[5] + "<br>");
        description.append("Target sequence length: " + tokens[6] + "<br>");
        description.append("Target start (0-based): " + tokens[7] + "<br>");
        description.append("Target end (0-based): " + tokens[8] + "<br>");
        description.append("Number of residual matches: " + tokens[9] + "<br>");
        description.append("Alignment block length: " + tokens[10] + "<br>");
        description.append("Mapping quality (0-255): " + tokens[11] + "<br>");

        for (int i = 12; i < tokens.length; i++) {
            String tag = tokens[i];
            if (!tag.startsWith("cg")) {    // Filter cigar string
                description.append(tag + "<br>");
            }
        }


        return new PAFFeature(chr, start, end, strand, score, name, description.toString());

    }


    @Override
    public boolean canDecode(String path) {
        String fn = path.toLowerCase();
        if (fn.endsWith(".gz")) fn = fn.substring(0, fn.length() - 3);
        return fn.endsWith(".paf");
    }
}
