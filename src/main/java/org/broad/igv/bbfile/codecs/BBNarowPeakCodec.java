package org.broad.igv.bbfile.codecs;

import org.broad.igv.Globals;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.feature.BasicFeature;
import org.broad.igv.feature.EncodePeakFeature;
import org.broad.igv.feature.Exon;
import org.broad.igv.feature.Strand;

import java.util.List;


/*

table bigNarrowPeak
"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
    string name;	 "Name given to a region (preferably unique). Use . if no name is assigned"
    uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
    char[1]  strand;     "+ or - or . for unknown"
    float  signalValue;  "Measurement of average enrichment for the region"
    float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
    float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
    int   peak;         "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."
)
 */

public class BBNarowPeakCodec extends BBBedCodec {

    public BBNarowPeakCodec(int standardFieldCount) {
        super(standardFieldCount);
    }

    @Override
    public EncodePeakFeature decode(BedFeature feat) {

        EncodePeakFeature feature = new EncodePeakFeature(feat.getChromosome(), feat.getStartBase(), feat.getEndBase());

        String[] restOfFields = feat.getRestOfFields();
        feature.setName(restOfFields[3 - 3]);
        feature.setScore(Float.parseFloat(restOfFields[4 - 3]));


        Strand strand;
        String strandString = restOfFields[5-3].trim();
        char strandChar = (strandString.length() == 0) ? ' ' : strandString.charAt(0);

        if (strandChar == '-') {
            strand = Strand.NEGATIVE;
        } else if (strandChar == '+') {
            strand = Strand.POSITIVE;
        } else {
            strand = Strand.NONE;
        }
        feature.setStrand(strand);


        // Store the remaining features in description string */
        feature.setSignal((float) Double.parseDouble(restOfFields[3]));
        feature.setPValue((float) Double.parseDouble(restOfFields[4]));
        feature.setQValue((float) Double.parseDouble(restOfFields[5]));
        feature.setPeakPosition(Integer.parseInt(restOfFields[6]));

        return feature;
    }
}
