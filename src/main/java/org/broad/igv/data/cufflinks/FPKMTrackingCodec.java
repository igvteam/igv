package org.broad.igv.data.cufflinks;

import org.broad.igv.logging.*;
import org.broad.igv.feature.Locus;
import org.broad.igv.util.ParsingUtils;

/**
 * Codec for Cufflinks FPKM files, extension fpkm_tracking
 *
* @author jacob
* @date 2013-Apr-18
*/
public class FPKMTrackingCodec extends CufflinksCodec<FPKMValue>{

    private static Logger log = LogManager.getLogger(FPKMTrackingCodec.class);

    int geneColumn = 4;
    int locusColumn  = 6;

    static final int startfpkmCol = 9;
    static final int colsPerSample = 4;

    private int numSamples = 1;

    public FPKMTrackingCodec(String path){
        super(FPKMValue.class, path);
    }

    @Override
    public Object readHeader(String[] tokens) {
        for(int i=0; i<tokens.length; i++) {
            final String tk = tokens[i];
            if(tk.equals("locus")) locusColumn = i;
            else if(tk.equals("gene_short_name")) geneColumn = i;
        }
        numSamples = (tokens.length - startfpkmCol) / colsPerSample;
        return tokens;
    }

    @Override
    public FPKMValue decode(String line) {
        return decode(ParsingUtils.TAB_PATTERN.split(line));
    }

    //@Override
    public FPKMValue decode(String[] tokens) {
        //Skip header line
        if (tokens[0].equalsIgnoreCase("tracking_id") || tokens[geneColumn].equalsIgnoreCase("gene_short_name")) {
            return null;
        }
        if (tokens.length >= (startfpkmCol + numSamples*colsPerSample)) {
            String locusString = tokens[locusColumn];
            if (locusString == null) return null;

            Locus locus = Locus.fromString(locusString);
            if(locus == null || locus.getChr() == null) return null;

            String gene = tokens[geneColumn];
            float[] fpkm = new float[numSamples];
            float[] confLo = new float[numSamples];
            float[] confHi = new float[numSamples];

            for(int sampNum = 0; sampNum < numSamples; sampNum++){
                int startCol = startfpkmCol + sampNum*colsPerSample;
                fpkm[sampNum] = Float.parseFloat(tokens[startCol]);
                confLo[sampNum] = Float.parseFloat(tokens[startCol+1]);
                confHi[sampNum] = Float.parseFloat(tokens[startCol+2]);
            }

            return new FPKMValue(locus.getChr(), locus.getStart() - 1, locus.getEnd(), gene,
                    fpkm, confLo, confHi);
        } else {
            log.error("Unexpected # of columns.  Expected at least 12,  found " + tokens.length);
            return null;
        }
    }

    public int getNumSamples() {
        return numSamples;
    }

    @Override
    public boolean canDecode(String path) {

        String fn = path.toLowerCase();
        if(fn.endsWith(".gz")) fn = fn.substring(0, fn.length()-3);
        return fn.endsWith("fpkm_tracking");
    }
}
