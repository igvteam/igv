package org.broad.igv.sam;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.WindowFunction;

import java.awt.*;

/**
 * Created by jrobinso on 6/17/14.
 * <p/>
 * id
 */
public class Ga4ghAlignment extends SAMAlignment {

    protected int alignmentStart;
    protected int alignmentEnd;
    int inferredInsertSize;
    int mappingQuality = 255;  // 255 by default
    String readName;
    protected String cigarString;
    protected String readSequence;

    public Ga4ghAlignment(JsonObject json) {

        String refName = json.get("referenceSequenceName").getAsString();

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        this.chr = genome == null ? refName : genome.getChromosomeAlias(refName);

        // SAMRecord is 1 based inclusive.  IGV is 0 based exclusive.
        this.readName = json.get("name").getAsString();
        this.flags = json.get("flags").getAsInt();
        this.alignmentStart = json.get("position").getAsInt() - 1;
        this.start = this.alignmentStart;   // might be modified later for soft clipping
        this.mappingQuality = json.get("mappingQuality").getAsInt();
        this.inferredInsertSize = json.get("templateLength").getAsInt();

        this.readSequence = json.get("originalBases").getAsString();


        this.cigarString = json.get("cigar").getAsString();
        this.alignmentEnd = this.alignmentStart + getReferenceLength(cigarString);
        this.end = alignmentEnd;

        if (isPaired()) {
            String mateReferenceName = json.get("mateReferenceSequenceName").getAsString();
            String mateChr = genome == null ? mateReferenceName : genome.getChromosomeAlias(mateReferenceName);
            this.setMate(new ReadMate(mateChr,
                    json.get("matePosition").getAsInt() - 1,
                    (this.flags & MATE_STRAND_FLAG) != 0,
                    (this.flags & MATE_UNMAPPED_FLAG) != 0));
        }




//        Object colorTag = record.getAttribute("YC");
//        if (colorTag != null) {
//            try {
//                color = ColorUtilities.stringToColor(colorTag.toString());
//            } catch (Exception e) {
//                log.error("Error interpreting color tag: " + colorTag, e);
//            }
//        }

        String baseQualities = json.get("baseQuality").toString();
        setPairOrientation();
        setPairStrands();
        createAlignmentBlocks(this.cigarString, this.readSequence.getBytes(), baseQualities.getBytes(),
                null, null, -1);

    }

    /**
     * @return the unclippedStart
     */
    public int getAlignmentStart() {
        return alignmentStart;
    }


    public int getAlignmentEnd() {

        return alignmentEnd;
    }


    public String getReadName() {
        return readName;
    }

    public int getMappingQuality() {
        return mappingQuality;
    }

    public int getInferredInsertSize() {
        return inferredInsertSize;
    }


    public String getCigarString() {
        return cigarString;
    }
    public int getReadLength() {
        return readSequence.length();
    }

    public String getReadSequence() {
        return readSequence;
    }

    // TODO -- implement
    @Override
    protected String getAttributeString(boolean truncate) {
        return "";
    }


    @Override
    public Object getAttribute(String key) {
        return null;
    }


    public int getReferenceLength(String cigarString) {
               // Use htsjdk class for now
        TextCigarCodec codec = new TextCigarCodec();
        Cigar cigar = codec.decode(cigarString);
        return cigar.getReferenceLength();
    }

    /*
    "id": "",
   "name": "SRR062634.3899514",
   "readsetId": "",
   "flags": 121,
   "referenceSequenceName": "1",
   "position": 70230915,
   "mappingQuality": 37,
   "cigar": "100M",
   "mateReferenceSequenceName": "1",
   "matePosition": 70230915,
   "templateLength": 0,
   "originalBases": "AATTGACTTTCCTTAGGCCACATGGCATGTTAGAGTCCTATCCAAAATTATTACCAAAGACTATGTATTAGGAAAGGAGAGGTTATGCTGCAGTAAAAAC",
   "alignedBases": "AATTGACTTTCCTTAGGCCACATGGCATGTTAGAGTCCTATCCAAAATTATTACCAAAGACTATGTATTAGGAAAGGAGAGGTTATGCTGCAGTAAAAAC",
   "baseQuality": "@ABD?DB?AEDHGILIEFKKDCFDE;\u003eKGFGKHKIEICBHHDJHQGPHFHH?CHHKGIFAFDFIFPJQDRRPPPFRQRQRQPQCPRDIQFCGHHID@@A/",
   "tags": {
    "AM": [
     "0"
    ],
     */
}
