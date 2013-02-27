package org.broad.igv.maf;

import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.index.Interval;
import org.broad.igv.util.index.IntervalTree;
import org.broad.tribble.util.SeekableStream;
import org.broad.tribble.util.SeekableStreamFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 2/18/13
 *         Time: 9:59 AM
 */
public class MAFParser implements MAFReader {

    String path;
    MAFIndex index;

    public MAFParser(String path) {
        this.path = path;

        String indexPath = path + ".index";
        try {
            if (ParsingUtils.pathExists(indexPath)) {
                index = MAFIndex.loadIndex(indexPath);
            } else {
                index = MAFIndex.createIndex(path);
                MAFIndex.writeIndex(index, indexPath);
            }
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }


    @Override
    public List<MultipleAlignmentBlock> loadAlignments(String chr, int start, int end) throws IOException {

        IntervalTree ivTree = index.getIntervalTree(chr);
        if (ivTree == null) return null;

        List<Interval> intervals = ivTree.findOverlapping(start, end);
        if (intervals.isEmpty()) {
            return null;
        }


        // Find the starting (left most) interval.  Alignment blocks do not overlap, so we can start at the
        // minimum file offset and just proceed until the end of the interval.
        long startPosition = Long.MAX_VALUE;
        for(Interval iv : intervals) {
            startPosition = Math.min(startPosition, iv.getValue());
        }


        SeekableStream is = null;


        is = SeekableStreamFactory.getStreamFor(path);
        is.seek(startPosition);

        BufferedReader reader = new BufferedReader(new InputStreamReader(is));

        List<MultipleAlignmentBlock> alignments = new ArrayList<MultipleAlignmentBlock>();

        String line;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("a ")) {
                // TODO -- parse score (optional)
                MultipleAlignmentBlock block = parseBlock(reader);
                if (block.getStart() > end || !block.getChr().equals(chr)) {
                    break;
                } else {
                    alignments.add(block);
                }
            }
        }
        return alignments;
    }


    @Override
    public Collection<String> getChrNames() {
        return index.getChromosomes();
    }

    @Override
    public Collection<String> getSpecies() {
        return index.getSpecies();  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getSpeciesName(String speciesId) {
        return speciesId;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getRefId() {
        return index.getRefId();  //To change body of implemented methods use File | Settings | File Templates.
    }


    void parseHeader() throws IOException {

        BufferedReader reader = null;

        try {
            InputStream is = SeekableStreamFactory.getStreamFor(path);
            reader = new BufferedReader(new InputStreamReader(is));

            List<MultipleAlignmentBlock> alignments = new ArrayList<MultipleAlignmentBlock>();

            String line;

            while ((line = reader.readLine()) != null) {

                if (line.startsWith("a ")) {
                    return;  // Done with header
                } else if (line.startsWith("track")) {
                    parseTrackLine(line);
                }
            }
        } finally {
            if (reader != null) reader.close();
        }
    }



    private void parseTrackLine(String line) {
        //To change body of created methods use File | Settings | File Templates.
    }


    /**
     * Parse an alignment block.  The reader has been advanced to the first sequencce
     *
     * @param reader
     */
    private MultipleAlignmentBlock parseBlock(BufferedReader reader) throws IOException {

        String line;
        MultipleAlignmentBlock ma = new MultipleAlignmentBlock();

        while ((line = reader.readLine()) != null) {
            if (line.trim().length() == 0) {
                return ma;
            }
            if (line.startsWith("s ")) {

                String[] tokens = Globals.whitespacePattern.split(line);

                String src = tokens[1];
                String species = src;
                String chr = src;
                if (src.contains(".")) {
                    String[] srcTokens = ParsingUtils.PERIOD_PATTERN.split(src);
                    species = srcTokens[0];
                    chr = srcTokens[1];
                }
                int start = Integer.parseInt(tokens[2]);
                int size = Integer.parseInt(tokens[3]);
                char strand = tokens[4].charAt(0);
                int srcSize = Integer.parseInt(tokens[5]);
                String text = tokens[6];

                ma.addSequence(new MultipleAlignmentBlock.Sequence(species, chr, start, size, strand, srcSize, text));


            }


        }
        return ma;
    }
}
