package org.broad.igv.sam.sbx;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.util.Locatable;
import org.broad.igv.feature.Strand;
import org.broad.igv.sam.*;
import org.broad.igv.sam.mods.BaseModificationSet;
import org.broad.igv.sam.mods.BaseModificationUtils;
import org.broad.igv.sam.smrt.SMRTKinetics;
import org.broad.igv.track.WindowFunction;

import java.awt.*;
import java.util.List;
import java.util.Map;

public class NullAlignment implements Alignment {

    static NullAlignment instance = new NullAlignment();

    public static NullAlignment getInstance() {
        return instance;
    }

    private NullAlignment() {
    }

    @Override
    public String getReadName() {
        return Alignment.super.getReadName();
    }

    @Override
    public String getReadSequence() {
        return Alignment.super.getReadSequence();
    }

    @Override
    public String getChr() {
        return Alignment.super.getChr();
    }

    @Override
    public int getAlignmentStart() {
        return 0;
    }

    @Override
    public int getAlignmentEnd() {
        return 0;
    }

    @Override
    public boolean contains(double location) {
        return false;
    }

    @Override
    public AlignmentBlock[] getAlignmentBlocks() {
        return new AlignmentBlock[0];
    }

    @Override
    public AlignmentBlock[] getInsertions() {
        return new AlignmentBlock[0];
    }

    @Override
    public Alignment trimSimplexTails() {
        return this;
    }

    @Override
    public String getCigarString() {
        return Alignment.super.getCigarString();
    }

    @Override
    public Cigar getCigar() {
        return Alignment.super.getCigar();
    }

    @Override
    public List<Gap> getGaps() {
        return List.of();
    }

    @Override
    public int getInferredInsertSize() {
        return 0;
    }

    @Override
    public int getLeadingHardClipLength() {
        return Alignment.super.getLeadingHardClipLength();
    }

    @Override
    public int getMappingQuality() {
        return 0;
    }

    @Override
    public ReadMate getMate() {
        return Alignment.super.getMate();
    }

    @Override
    public Strand getReadStrand() {
        return null;
    }

    @Override
    public boolean isProperPair() {
        return false;
    }

    @Override
    public boolean isMapped() {
        return false;
    }

    @Override
    public boolean isPaired() {
        return false;
    }

    @Override
    public boolean isFirstOfPair() {
        return false;
    }

    @Override
    public boolean isSecondOfPair() {
        return false;
    }

    @Override
    public boolean isNegativeStrand() {
        return false;
    }

    @Override
    public boolean isDuplicate() {
        return false;
    }

    @Override
    public boolean isPrimary() {
        return false;
    }

    @Override
    public boolean isSupplementary() {
        return false;
    }

    @Override
    public byte getBase(double position) {
        return 0;
    }

    @Override
    public byte getPhred(double position) {
        return 0;
    }

    @Override
    public Object getAttribute(String key) {
        return Alignment.super.getAttribute(key);
    }

    @Override
    public void setMateSequence(String sequence) {

    }

    @Override
    public String getPairOrientation() {
        return "";
    }

    @Override
    public Strand getFirstOfPairStrand() {
        return null;
    }

    @Override
    public Strand getSecondOfPairStrand() {
        return null;
    }

    @Override
    public boolean isVendorFailedRead() {
        return false;
    }

    @Override
    public Color getYcColor() {
        return null;
    }

    @Override
    public String getSample() {
        return Alignment.super.getSample();
    }

    @Override
    public String getReadGroup() {
        return Alignment.super.getReadGroup();
    }

    @Override
    public String getLibrary() {
        return Alignment.super.getLibrary();
    }

    @Override
    public String getClipboardString(double location, int mouseX) {
        return Alignment.super.getClipboardString(location, mouseX);
    }

    @Override
    public void finish() {
        Alignment.super.finish();
    }

    @Override
    public AlignmentBlock getInsertionAt(int position) {
        return Alignment.super.getInsertionAt(position);
    }

    @Override
    public Gap getDeletionAt(int position) {
        return Alignment.super.getDeletionAt(position);
    }

    @Override
    public ClippingCounts getClippingCounts() {
        return Alignment.super.getClippingCounts();
    }

    @Override
    public void setClusterName(String hap) {
        Alignment.super.setClusterName(hap);
    }

    @Override
    public String getClusterName() {
        return Alignment.super.getClusterName();
    }

    @Override
    public void setHapDistance(int dist) {
        Alignment.super.setHapDistance(dist);
    }

    @Override
    public int getClusterDistance() {
        return Alignment.super.getClusterDistance();
    }

    @Override
    public Map<Integer, BaseModificationUtils> getBaseModificationMap() {
        return Alignment.super.getBaseModificationMap();
    }

    @Override
    public List<BaseModificationSet> getBaseModificationSets() {
        return Alignment.super.getBaseModificationSets();
    }

    @Override
    public SMRTKinetics getSmrtKinetics() {
        return Alignment.super.getSmrtKinetics();
    }

    @Override
    public String getAlignmentValueString(double position, int mouseX, AlignmentTrack.RenderOptions renderOptions) {
        return Alignment.super.getAlignmentValueString(position, mouseX, renderOptions);
    }

    @Override
    public Alignment getSpecificAlignment(double location) {
        return Alignment.super.getSpecificAlignment(location);
    }

    @Override
    public float getScore() {
        return 0;
    }

    @Override
    public String getValueString(double position, int mouseX, WindowFunction windowFunction) {
        return Alignment.super.getValueString(position, mouseX, windowFunction);
    }

    @Override
    public String getContig() {
        return "";
    }

    @Override
    public int getStart() {
        return 0;
    }

    @Override
    public int getEnd() {
        return 0;
    }

    @Override
    public int getLengthOnReference() {
        return Alignment.super.getLengthOnReference();
    }

    @Override
    public boolean overlaps(Locatable other) {
        return Alignment.super.overlaps(other);
    }

    @Override
    public boolean withinDistanceOf(Locatable other, int distance) {
        return Alignment.super.withinDistanceOf(other, distance);
    }

    @Override
    public boolean contains(Locatable other) {
        return Alignment.super.contains(other);
    }

    @Override
    public boolean contigsMatch(Locatable other) {
        return Alignment.super.contigsMatch(other);
    }
}
