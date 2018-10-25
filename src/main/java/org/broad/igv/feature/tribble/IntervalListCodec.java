package org.broad.igv.feature.tribble;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;
import org.broad.igv.Globals;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;

import java.io.IOException;

/**
 * Created by jrobinso on 8/13/18.
 */
public class IntervalListCodec extends AsciiFeatureCodec<BasicFeature> {

    Genome genome;

    protected IntervalListCodec(Genome genome) {
        super(BasicFeature.class);
        this.genome = genome;
    }

    @Override
    public BasicFeature decode(String line) {

        if (line.startsWith("@")) {
            return null;
        } else {
            //chr11	37781917	40267815	+	name1
            
            final String[] tokens = Globals.whitespacePattern.split(line);

            String chr = tokens[0];
            if(genome != null) {
                chr = genome.getCanonicalChrName(chr);
            }

            final int start = Integer.parseInt(tokens[1]) - 1;
            final int end = Integer.parseInt(tokens[2]);

            final BasicFeature feature = new BasicFeature(chr, start, end);

            if (tokens.length > 3) {
                final String strandString = tokens[3];
                final Strand strand = Strand.fromString(strandString);
                feature.setStrand(strand);
            }

            if (tokens.length > 4) {
                feature.setName(tokens[4]);
            }

            return feature;
        }
    }

    @Override
    public Object readActualHeader(LineIterator lineIterator) {
        return null;
    }

    @Override
    public boolean canDecode(String s) {
        return s.endsWith("interval_list");
    }
}
