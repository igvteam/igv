package org.broad.igv.maf;

import org.broad.igv.Globals;
import org.broad.igv.util.ParsingUtils;
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


    private void parseTrackLine(String line) {
        //To change body of created methods use File | Settings | File Templates.
    }

    @Override
    public List<MultipleAlignmentBlock> loadAligments(String chr, int start, int end, List<String> species) throws IOException {
        return parse();
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


    public List<MultipleAlignmentBlock> parse() throws IOException {

        BufferedReader reader = null;


        InputStream is = SeekableStreamFactory.getStreamFor(path);
        reader = new BufferedReader(new InputStreamReader(is));

        List<MultipleAlignmentBlock> alignments = new ArrayList<MultipleAlignmentBlock>();

        String line = reader.readLine();
        if (line.startsWith("track")) {
            parseTrackLine(line);
            line = reader.readLine();
        }

        do {
            if (line.startsWith("#")) continue;

            if (line.startsWith("a ")) {
                // TODO -- parse score (optional)
                alignments.add(parseBlock(reader));
            }
        } while ((line = reader.readLine()) != null);

        return alignments;
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
