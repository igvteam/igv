/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.ga4gh;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.sam.ReadMate;
import org.broad.igv.sam.SAMAlignment;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

/**
 * Created by jrobinso on 6/17/14.
 * <p/>
 * id
 */
public class Ga4ghAlignment extends SAMAlignment {

    private final Map<String, String> tags;
    protected int alignmentStart;
    protected int alignmentEnd;
    int inferredInsertSize;
    int mappingQuality = 255;  // 255 by default
    String readName;
    protected String readSequence;
    private boolean negativeStrand;
    private int readNumber;
    private boolean duplicateFragment;
    private int numberReads;
    private boolean properPlacement;
    private boolean supplementaryAlignment;
    private boolean failedVendorQualityChecks;
    private boolean secondaryAlignment;
    private String cigarString;
    private boolean mapped;

    public Ga4ghAlignment(JsonObject json) {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();

        this.readName = json.get("fragmentName").getAsString();
        this.properPlacement = json.has("properPlacement") ? json.get("properPlacement").getAsBoolean() : true;
        this.duplicateFragment = json.has("duplicateFragment") ? json.get("duplicateFragment").getAsBoolean() : false;
        this.numberReads = json.has("numberReads") ? json.get("numberReads").getAsInt() : 1;
        this.inferredInsertSize = json.has("fragmentLength") ? json.get("fragmentLength").getAsInt() : 0;
        this.readNumber = json.has("readNumber") ? json.get("readNumber").getAsInt() : 0;
        this.failedVendorQualityChecks = json.has("failedVendorQualityChecks") ? json.get("failedVendorQualityChecks").getAsBoolean() : false;

        JsonObject alignmentObject = json.getAsJsonObject("alignment");
        if (alignmentObject == null) {
            this.mapped = false;
        } else {
            this.mapped = true;

            JsonObject positionObject = alignmentObject.getAsJsonObject("position");
            String refName = positionObject.get("referenceName").getAsString();
            this.setChr(genome == null ? refName : genome.getChromosomeAlias(refName));


            this.alignmentStart = positionObject.get("position").getAsInt();
            this.mappingQuality = alignmentObject.has("mappingQuality") ? alignmentObject.get("mappingQuality").getAsInt() : 256;
            this.negativeStrand = positionObject.get("reverseStrand").getAsBoolean();

            this.cigarString = generateCigarString(alignmentObject.getAsJsonArray("cigar"));
            this.start = this.alignmentStart;   // might be modified later for soft clipping
            this.alignmentEnd = this.alignmentStart + getReferenceLength(cigarString);
            this.end = alignmentEnd;

        }

        this.secondaryAlignment = json.has("secondaryAlignment") ? json.get("secondaryAlignment").getAsBoolean() : false;
        this.supplementaryAlignment = json.has("supplementaryAlignment") ? json.get("supplementaryAlignment").getAsBoolean() : false;
        this.readSequence = json.has("alignedSequence") ? json.get("alignedSequence").getAsString() : null;
        byte[] baseQualities = json.has("alignedQuality") ? generateBaseQualities(json.getAsJsonArray("alignedQuality")) : null;

        JsonObject mateObject = json.getAsJsonObject("nextMatePosition");
        if (mateObject == null) {
            this.setMate(new ReadMate("*", 0, false, true));
        } else {
            String mateReferenceName = mateObject.get("referenceName").getAsString();
            String mateChr = genome == null ? mateReferenceName : genome.getChromosomeAlias(mateReferenceName);
            int matePosition = Integer.parseInt(mateObject.get("position").getAsString());
            boolean mateNegStrand = mateObject.get("reverseStrand").getAsBoolean();
            this.setMate(new ReadMate(mateChr,
                    matePosition,
                    mateNegStrand,
                    false));           // Assuming mate is mapped

        }

        JsonObject infoObject = json.getAsJsonObject("info");

        this.tags = generateTags(infoObject);
        setPairOrientation();
        setPairStrands();
        createAlignmentBlocks(this.cigarString, this.readSequence.getBytes(), baseQualities, null, null, -1);
    }

    private Map<String, String> generateTags(JsonObject infoObject) {

        Map<String, String> tags = new HashMap<String, String>();
        if (infoObject != null) {
            for (Map.Entry<String, JsonElement> entry : infoObject.entrySet()) {
                String key = entry.getKey();

                JsonArray valueArray = entry.getValue().getAsJsonArray();
                String value = valueArray.get(0).getAsString();

                for (int i = 1; i < valueArray.size(); i++) {
                    value += "," + valueArray.get(i).getAsString();
                }
                tags.put(key, value);
            }
        }

        return tags;

    }

    private byte[] generateBaseQualities(JsonArray alignedQuality) {

        byte[] baseQualities = new byte[alignedQuality.size()];
        Iterator<JsonElement> iter = alignedQuality.iterator();
        int i = 0;
        while (iter.hasNext()) {
            baseQualities[i++] = iter.next().getAsByte();
        }

        return baseQualities;
    }

    private String generateCigarString(JsonArray cigar) {

        StringBuffer cigarStr = new StringBuffer();

        Iterator<JsonElement> iter = cigar.iterator();

        while (iter.hasNext()) {
            JsonObject op = iter.next().getAsJsonObject();
            cigarStr.append(op.getAsJsonPrimitive("operationLength").getAsString());
            cigarStr.append(CigarMap.get(op.getAsJsonPrimitive("operation").getAsString()));
        }
        return cigarStr.toString();

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

    @Override
    protected String getAttributeString(boolean truncate) {

        StringBuffer buffer = new StringBuffer();
        for (Map.Entry<String, String> entry : tags.entrySet()) {
            buffer.append(entry.getKey() + ": " + entry.getValue());
        }
        return buffer.toString();
    }

    @Override
    public boolean isFirstOfPair() {
        return readNumber == 0;
    }

    @Override
    public boolean isSecondOfPair() {
        return readNumber > 0;
    }

    @Override
    public boolean isDuplicate() {
        return duplicateFragment;
    }

    @Override
    public boolean isMapped() {
        return mapped;
    }

    @Override
    public boolean isPaired() {
        return numberReads > 1;
    }

    @Override
    public boolean isProperPair() {
        return properPlacement;
    }

    @Override
    public boolean isSupplementary() {
        return supplementaryAlignment;
    }

    @Override
    public boolean isVendorFailedRead() {
        return failedVendorQualityChecks;
    }

    @Override
    public boolean isPrimary() {
        return !secondaryAlignment;
    }


    @Override
    public Object getAttribute(String key) {
        return tags.get(key);
    }


    public int getReferenceLength(String cigarString) {
        // Use htsjdk class for now
        TextCigarCodec codec = new TextCigarCodec();
        Cigar cigar = codec.decode(cigarString);
        return cigar.getReferenceLength();
    }

    public boolean isNegativeStrand() {
        return negativeStrand;

    }

    static Map<String, String> CigarMap = new HashMap<String, String>();

    static {
        CigarMap.put("ALIGNMENT_MATCH", "M");
        CigarMap.put("INSERT", "I");
        CigarMap.put("DELETE", "D");
        CigarMap.put("SKIP", "N");
        CigarMap.put("CLIP_SOFT", "S");
        CigarMap.put("CLIP_HARD", "H");
        CigarMap.put("PAD", "P");
        CigarMap.put("SEQUENCE_MATCH", "=");
        CigarMap.put("SEQUENCE_MISMATCH", "X");
    }

    ;
}
