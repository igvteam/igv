package org.igv.feature.tribble;

import org.igv.Globals;
import org.igv.feature.IGVFeature;
import org.igv.feature.UCSCSnpFeature;
import org.igv.feature.genome.Genome;

/**
 * Created by jrobinso on 5/26/15.
 */
public class UCSCSnpCodec extends UCSCCodec<IGVFeature> {

    private Genome genome;

    public UCSCSnpCodec(Genome genome) {
        super(UCSCSnpFeature.class);
        this.genome = genome;
    }

    @Override
    public UCSCSnpFeature decode(String [] tokens) {

        if (tokens.length < 25) return null;

        String chr = tokens[1];
        if (genome != null) {
            chr = genome.getCanonicalChrName(chr);
        }
        int start = Integer.parseInt(tokens[2]);
        int end = Integer.parseInt(tokens[3]);
        return new UCSCSnpFeature(chr, start, end, tokens);

    }

    @Override
    public boolean canDecode(String path) {
        String fn = path.toLowerCase();
        if(fn.endsWith(".gz")) fn = fn.substring(0, fn.length()-3);
        return fn.toLowerCase().endsWith(".snp");
    }
}
