package org.broad.igv.feature.tribble;

import org.broad.igv.Globals;
import org.broad.igv.feature.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.tribble.Feature;
import org.broad.tribble.util.ParsingUtils;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Codec for UCSC / ENCDOE "broad and narrow peak" files (http://genome.ucsc.edu/FAQ/FAQformat.html#format13)
 * <p/>
 * This format is used to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. It is a BED6+4 format.
 * <p/>
 * chrom - Name of the chromosome (or contig, scaffold, etc.).
 * chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
 * chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
 * name - Name given to a region (preferably unique). Use '.' if no name is assigned.
 * score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
 * strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
 * signalValue - Measurement of overall (usually, average) enrichment for the region.
 * pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
 * qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
 * peak (Narrow peak only) - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.
 *
 * @author jrobinso
 *         Date: 10/16/12
 *         Time: 11:35 AM
 */
public class EncodePeakCodec extends UCSCCodec {

    Genome genome;

    public EncodePeakCodec() {
        this(null);
    }

    public EncodePeakCodec(Genome genome) {
        super(BasicFeature.class);
        this.genome = genome;
    }

    @Override
    public Feature decode(String nextLine) {


        if (nextLine.trim().length() == 0) {
            return null;
        }

        if (nextLine.startsWith("#") || nextLine.startsWith("track") || nextLine.startsWith("browser")) {
            return null;
        }

        String[] tokens = Globals.tabPattern.split(nextLine);

        int tokenCount = tokens.length;

        // The first 3 columns are non optional for BED.  We will relax this
        // and only require 2.

        if (tokenCount < 9) {
            return null;
        }

        String c = tokens[0];
        String chr = genome == null ? c : genome.getChromosomeAlias(c);

        //BED format, and IGV, use starting element as 0.
        int start = Integer.parseInt(tokens[1]);
        int end = Integer.parseInt(tokens[2]);
        BasicFeature feature = new BasicFeature(chr, start, end);

        feature.setName(tokens[3]);
        feature.setScore(Float.parseFloat(tokens[4]));


        Strand strand;
        String strandString = tokens[5].trim();
        char strandChar = (strandString.length() == 0) ? ' ' : strandString.charAt(0);

        if (strandChar == '-') {
            feature.setStrand(Strand.NEGATIVE);
        } else if (strandChar == '+') {
            feature.setStrand(Strand.POSITIVE);
        } else {
            feature.setStrand(Strand.NONE);
        }


        // Store the remaining features in description string */
        StringBuffer desc = new StringBuffer();
        desc.append("Signal value: " + tokens[6]);
        desc.append("<br>P value: " + tokens[7]);
        desc.append("<br>Q value: " + tokens[8]);
        if (tokens.length > 9) {
            int peak = Integer.parseInt(tokens[9]);
            if (peak >= 0) {
                desc.append("<br>Peak: " + (feature.getStart() + peak + 1));
            }
        }
        feature.setDescription(desc.toString());

        return feature;
    }

    @Override
    public boolean canDecode(String s) {
        return s.toLowerCase().contains("narrowpeak");
    }
}
