/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.ga4gh;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.apache.log4j.Logger;
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

    private static Logger log = Logger.getLogger(Ga4ghAlignment.class);

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
        this.properPlacement = hasNonNullValue(json, "properPlacement") ? json.get("properPlacement").getAsBoolean() : true;
        this.duplicateFragment = hasNonNullValue(json, "duplicateFragment") ? json.get("duplicateFragment").getAsBoolean() : false;
        this.numberReads = hasNonNullValue(json, "numberReads") ? json.get("numberReads").getAsInt() : 1;
        this.inferredInsertSize = hasNonNullValue(json, "fragmentLength") ? json.get("fragmentLength").getAsInt() : 0;
        this.readNumber = hasNonNullValue(json, "readNumber") ? json.get("readNumber").getAsInt() : 0;
        this.failedVendorQualityChecks = hasNonNullValue(json, "failedVendorQualityChecks") ? json.get("failedVendorQualityChecks").getAsBoolean() : false;

        JsonObject alignmentObject = json.getAsJsonObject("alignment");
        if (alignmentObject == null) {
            this.mapped = false;
        } else {
            this.mapped = true;

            JsonObject positionObject = alignmentObject.getAsJsonObject("position");
            String refName = positionObject.get("referenceName").getAsString();
            this.setChr(genome == null ? refName : genome.getChromosomeAlias(refName));

            this.alignmentStart = positionObject.get("position").getAsInt();
            this.mappingQuality = hasNonNullValue(alignmentObject, "mappingQuality") ? alignmentObject.get("mappingQuality").getAsInt() : 256;
            this.negativeStrand = hasNonNullValue(positionObject, "reverseStrand") && positionObject.get("reverseStrand").getAsBoolean();

            this.cigarString = generateCigarString(alignmentObject.getAsJsonArray("cigar"));
            this.start = this.alignmentStart;   // might be modified later for soft clipping
            this.alignmentEnd = this.alignmentStart + getReferenceLength(cigarString);
            this.end = alignmentEnd;

        }

        this.secondaryAlignment = hasNonNullValue(json, "secondaryAlignment") ? json.get("secondaryAlignment").getAsBoolean() : false;
        this.supplementaryAlignment = hasNonNullValue(json, "supplementaryAlignment") ? json.get("supplementaryAlignment").getAsBoolean() : false;
        this.readSequence = hasNonNullValue(json, "alignedSequence") ? json.get("alignedSequence").getAsString() : null;
        byte[] baseQualities = hasNonNullValue(json, "alignedQuality") ? generateBaseQualities(json.getAsJsonArray("alignedQuality")) : null;

        JsonObject mateObject = json.getAsJsonObject("nextMatePosition");
        if (mateObject == null) {
            this.setMate(new ReadMate("*", 0, false, true));
        } else {
            String mateReferenceName = mateObject.get("referenceName").getAsString();
            String mateChr = genome == null ? mateReferenceName : genome.getChromosomeAlias(mateReferenceName);
            int matePosition = Integer.parseInt(mateObject.get("position").getAsString());
            boolean mateNegStrand = hasNonNullValue(mateObject, "reverseStrand") && mateObject.get("reverseStrand").getAsBoolean();
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

    public static boolean hasNonNullValue(JsonObject json, String name) {
        return json.has(name) && !json.get(name).isJsonNull();
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
