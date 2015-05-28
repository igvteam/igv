package org.broad.igv.feature.tribble;

import org.broad.igv.Globals;
import org.broad.igv.feature.UCSCSnpFeature;
import org.broad.igv.feature.genome.Genome;

/**
 * Created by jrobinso on 5/26/15.
 */
public class UCSCSnpCodec extends UCSCCodec<UCSCSnpFeature> {

    private Genome genome;

    public UCSCSnpCodec(Genome genome) {
        super(UCSCSnpFeature.class);
        this.genome = genome;
    }

    @Override
    public UCSCSnpFeature decode(String s) {

        String[] tokens = Globals.tabPattern.split(s);

        if (tokens.length < 25) return null;

        String chr = tokens[1];
        if (genome != null) {
            chr = genome.getChromosomeAlias(chr);
        }
        int start = Integer.parseInt(tokens[2]);
        int end = Integer.parseInt(tokens[3]);
        return new UCSCSnpFeature(chr, start, end, tokens);

    }
}
