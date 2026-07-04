package org.igv.alignment;

import org.igv.feature.Range;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.prefs.Constants;
import org.igv.prefs.IGVPreferences;
import org.igv.prefs.PreferencesManager;
import org.igv.alignment.mods.BaseModficationFilter;
import org.igv.track.Track;
import org.igv.util.collections.CollUtils;
import org.json.JSONObject;
import org.w3c.dom.Element;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

import static org.igv.prefs.Constants.*;

public class RenderOptions implements Cloneable {

    public static final String NAME = "RenderOptions";
    private static final Logger log = LogManager.getLogger(RenderOptions.class);

    private AlignmentTrack track;
    private Boolean shadeBasesOption;
    private Boolean shadeCenters;
    private Boolean showCenterline;
    private Boolean flagUnmappedPairs;
    private Boolean showAllBases;
    private Integer minInsertSize;
    private Integer maxInsertSize;
    private AlignmentTrack.ColorOption colorOption;
    private SortOption sortOption;
    private AlignmentTrack.GroupOption groupByOption;
    private AlignmentTrack.ShadeAlignmentsOption shadeAlignmentsOption;
    private AlignmentTrack.DuplicatesOption duplicatesOption;
    private Integer mappingQualityLow;
    private Integer mappingQualityHigh;
    private boolean viewPairs = false;
    private String colorByTag;
    private String groupByTag;
    private String sortByTag;
    private String linkByTag;
    private Boolean linkedReads;
    private Boolean quickConsensusMode;
    private Boolean showMismatches;
    private Boolean indelQualColoring;
    private Boolean indelQualUsesMin;
    private Boolean indelQualSbx;
    private Boolean tailQualSbx;
    private Boolean hideTailSbx;
    private Boolean insertQualColoring;
    Boolean computeIsizes;
    private Double minInsertSizePercentile;
    private Double maxInsertSizePercentile;
    private Boolean pairedArcView;
    private Boolean flagZeroQualityAlignments;
    private Range groupByPos;
    private Boolean invertSorting;
    private boolean invertGroupSorting;
    private Boolean hideSmallIndels;
    private Integer smallIndelThreshold;
    private BaseModficationFilter basemodFilter;
    private Float basemodThreshold;
    private Boolean basemodDistinguishStrands;
    private int baseQualityMin;
    private int baseQualityMax;

    private Integer minJunctionCoverage;


    AlignmentTrack.BisulfiteContext bisulfiteContext = AlignmentTrack.BisulfiteContext.CG;
    Map<String, PEStats> peStats;

    RenderOptions(AlignmentTrack track) {
        this.track = track;
        peStats = new HashMap<>();

        // Set some constants -- for efficiency
        this.baseQualityMin = track == null ? 5 : track.getPreferences().getAsInt(SAM_BASE_QUALITY_MIN);
        this.baseQualityMax = track == null ? 20 : track.getPreferences().getAsInt(SAM_BASE_QUALITY_MAX);
    }

    IGVPreferences getPreferences() {
        return this.track != null ? this.track.getPreferences() : AlignmentTrack.getPreferences(AlignmentTrack.ExperimentType.OTHER);
    }

    public Boolean getShowCenterline() {
        return showCenterline;
    }

    public void setShowCenterline(Boolean showCenterline) {
        this.showCenterline = showCenterline;
    }

    public int getMinJunctionCoverage() {
        return minJunctionCoverage != null ? minJunctionCoverage : PreferencesManager.getPreferences(Constants.RNA).getAsInt(SAM_JUNCTION_MIN_COVERAGE);
    }

    public void setMinJunctionCoverage(int minJunctionCoverage) {
        this.minJunctionCoverage = minJunctionCoverage;
    }

    public int getBaseQualityMin() {
        return baseQualityMin;
    }

    public int getBaseQualityMax() {
        return baseQualityMax;
    }

    public HashMap<String, Color> getSelectedReadNames() {
        return this.track.getSelectedReadNames();
    }

    public Track.DisplayMode getDisplayMode() {
        return this.track.getDisplayMode();
    }

    void setShowAllBases(boolean showAllBases) {
        this.showAllBases = showAllBases;
    }

    void setShowMismatches(boolean showMismatches) {
        this.showMismatches = showMismatches;
    }

    void setMinInsertSize(int minInsertSize) {
        this.minInsertSize = minInsertSize;
        //updateColorScale();
    }

    public void setViewPairs(boolean viewPairs) {
        this.viewPairs = viewPairs;
    }

    void setIndelQualColoring(boolean indelQualColoring) {
        this.indelQualColoring = indelQualColoring;
    }

    void setIndelQualUsesMin(boolean indelQualUsesMin) {
        this.indelQualUsesMin = indelQualUsesMin;
    }

    void setComputeIsizes(boolean computeIsizes) {
        this.computeIsizes = computeIsizes;
    }


    void setMaxInsertSizePercentile(double maxInsertSizePercentile) {
        this.maxInsertSizePercentile = maxInsertSizePercentile;
    }

    void setMaxInsertSize(int maxInsertSize) {
        this.maxInsertSize = maxInsertSize;
    }

    void setMinInsertSizePercentile(double minInsertSizePercentile) {
        this.minInsertSizePercentile = minInsertSizePercentile;
    }

    void setColorByTag(String colorByTag) {
        this.colorByTag = colorByTag;
    }

    void setColorOption(AlignmentTrack.ColorOption colorOption) {
        this.colorOption = colorOption;
    }

    void setSortOption(SortOption sortOption) {
        this.sortOption = sortOption;
    }

    void setSortByTag(String sortByTag) {
        this.sortByTag = sortByTag;
    }

    void setGroupByTag(String groupByTag) {
        this.groupByTag = groupByTag;
    }

    void setGroupByPos(Range groupByPos) {
        this.groupByPos = groupByPos;
    }

    void setInvertSorting(boolean invertSorting) {
        this.invertSorting = invertSorting;
    }

    void setInvertGroupSorting(boolean invertGroupSorting) {
        this.invertGroupSorting = invertGroupSorting;
    }

    void setLinkByTag(String linkByTag) {
        this.linkByTag = linkByTag;
    }

    void setQuickConsensusMode(boolean quickConsensusMode) {
        this.quickConsensusMode = quickConsensusMode;
    }

    public void setGroupByOption(AlignmentTrack.GroupOption groupByOption) {
        this.groupByOption = (groupByOption == null) ? AlignmentTrack.GroupOption.NONE : groupByOption;
    }

    void setShadeAlignmentsOption(AlignmentTrack.ShadeAlignmentsOption shadeAlignmentsOption) {
        this.shadeAlignmentsOption = shadeAlignmentsOption;
    }

    void setShadeBasesOption(boolean shadeBasesOption) {
        this.shadeBasesOption = shadeBasesOption;
    }

    void setLinkedReads(boolean linkedReads) {
        this.linkedReads = linkedReads;
    }

    public void setHideSmallIndels(boolean hideSmallIndels) {
        this.hideSmallIndels = hideSmallIndels;
    }

    public void setSmallIndelThreshold(int smallIndelThreshold) {
        this.smallIndelThreshold = smallIndelThreshold;
    }

    // getters
    public int getMinInsertSize() {
        return minInsertSize == null ? getPreferences().getAsInt(SAM_MIN_INSERT_SIZE_THRESHOLD) : minInsertSize;
    }

    public int getMaxInsertSize() {
        return maxInsertSize == null ? getPreferences().getAsInt(SAM_MAX_INSERT_SIZE_THRESHOLD) : maxInsertSize;
    }

    public boolean isFlagUnmappedPairs() {
        return flagUnmappedPairs == null ? getPreferences().getAsBoolean(SAM_FLAG_UNMAPPED_PAIR) : flagUnmappedPairs;
    }

    public boolean getShadeBasesOption() {
        return shadeBasesOption == null ? getPreferences().getAsBoolean(SAM_SHADE_BASES) : shadeBasesOption;
    }

    public boolean isShowMismatches() {
        return showMismatches == null ? getPreferences().getAsBoolean(SAM_SHOW_MISMATCHES) : showMismatches;
    }

    public boolean isShowAllBases() {
        return showAllBases == null ? getPreferences().getAsBoolean(SAM_SHOW_ALL_BASES) : showAllBases;
    }

    public boolean isShadeCenters() {
        return shadeCenters == null ? getPreferences().getAsBoolean(SAM_SHADE_CENTER) : shadeCenters;
    }

    public boolean isFlagZeroQualityAlignments() {
        return flagZeroQualityAlignments == null ? getPreferences().getAsBoolean(SAM_FLAG_ZERO_QUALITY) : flagZeroQualityAlignments;
    }

    public boolean isViewPairs() {
        return viewPairs;
    }

    public boolean isIndelQualColoring() {
        return indelQualColoring == null ? getPreferences().getAsBoolean(SAM_INDEL_QUAL_COLORING) : indelQualColoring;
    }

    public boolean isIndelQualUsesMin() {
        return indelQualUsesMin == null ? getPreferences().getAsBoolean(SAM_INDEL_QUAL_USES_MIN) : indelQualUsesMin;
    }


    // SBX Options
    public boolean isIndelQualSbx() {
        return AlignmentTrack.ExperimentType.SBX == track.experimentType && (indelQualSbx == null ? getPreferences().getAsBoolean(SAM_INDEL_QUAL_SBX) : indelQualSbx);
    }

    public void setTailQualSbx(Boolean tailQualSbx) {
        this.tailQualSbx = tailQualSbx;
    }

    public boolean isTailQualSbx() {
        return AlignmentTrack.ExperimentType.SBX == track.experimentType && (tailQualSbx == null ? getPreferences().getAsBoolean(SAM_TAIL_QUAL_SBX) : tailQualSbx);
    }

    public void setHideTailSbx(Boolean hideTailSbx) {
        this.hideTailSbx = hideTailSbx;
    }

    public boolean isHideTailSbx() {
        return AlignmentTrack.ExperimentType.SBX == track.experimentType && (hideTailSbx == null ? getPreferences().getAsBoolean(SAM_HIDE_TAIL_SBX) : hideTailSbx);
    }

    public void setIndelQualSbx(Boolean indelQualSbx) {
        this.indelQualSbx = indelQualSbx;
    }
    // End SBX options

    public boolean isComputeIsizes() {
        return computeIsizes == null ? getPreferences().getAsBoolean(SAM_COMPUTE_ISIZES) : computeIsizes;
    }

    public double getMinInsertSizePercentile() {
        return minInsertSizePercentile == null ? getPreferences().getAsFloat(SAM_MIN_INSERT_SIZE_PERCENTILE) : minInsertSizePercentile;
    }

    public double getMaxInsertSizePercentile() {
        return maxInsertSizePercentile == null ? getPreferences().getAsFloat(SAM_MAX_INSERT_SIZE_PERCENTILE) : maxInsertSizePercentile;
    }

    public AlignmentTrack.ColorOption getColorOption() {
        return colorOption == null ?
                CollUtils.valueOf(AlignmentTrack.ColorOption.class, getPreferences().get(SAM_COLOR_BY), AlignmentTrack.ColorOption.NONE) :
                colorOption;
    }

    public String getColorByTag() {
        return colorByTag == null ? getPreferences().get(SAM_COLOR_BY_TAG) : colorByTag;
    }

    public AlignmentTrack.ShadeAlignmentsOption getShadeAlignmentsOption() {
        if (shadeAlignmentsOption != null) {
            return shadeAlignmentsOption;
        } else {
            try {
                return AlignmentTrack.ShadeAlignmentsOption.valueOf(getPreferences().get(SAM_SHADE_ALIGNMENT_BY));
            } catch (IllegalArgumentException e) {
                log.error("Error parsing alignment shade option: " + AlignmentTrack.ShadeAlignmentsOption.valueOf(getPreferences().get(SAM_SHADE_ALIGNMENT_BY)));
                return AlignmentTrack.ShadeAlignmentsOption.NONE;
            }
        }
    }

    public AlignmentTrack.DuplicatesOption getDuplicatesOption() {
        final IGVPreferences prefs = getPreferences();
        if (duplicatesOption != null) {
            return duplicatesOption;
        } else {
            duplicatesOption = prefs.getAsBoolean(SAM_FILTER_DUPLICATES)
                    ? AlignmentTrack.DuplicatesOption.FILTER
                    : AlignmentTrack.DuplicatesOption.SHOW;
        }
        return duplicatesOption;
    }

    public void setDuplicatesOption(final AlignmentTrack.DuplicatesOption duplicatesOption) {
        this.duplicatesOption = duplicatesOption;
    }

    public int getMappingQualityLow() {
        return mappingQualityLow == null ? getPreferences().getAsInt(SAM_SHADE_QUALITY_LOW) : mappingQualityLow;
    }

    public int getMappingQualityHigh() {
        return mappingQualityHigh == null ? getPreferences().getAsInt(SAM_SHADE_QUALITY_HIGH) : mappingQualityHigh;
    }

    SortOption getSortOption() {
        return sortOption == null ? SortOption.fromString(getPreferences().get(SAM_SORT_OPTION)) : sortOption;
    }

    String getSortByTag() {
        return sortByTag == null ? getPreferences().get(SAM_SORT_BY_TAG) : sortByTag;
    }

    public String getGroupByTag() {
        return groupByTag == null ? getPreferences().get(SAM_GROUP_BY_TAG) : groupByTag;
    }

    public Range getGroupByPos() {
        if (groupByPos == null) {
            String pos = getPreferences().get(SAM_GROUP_BY_POS);
            if (pos != null) {
                String[] posParts = pos.split(" ");
                if (posParts.length != 2) {
                    groupByPos = null;
                } else {
                    int posChromStart = Integer.parseInt(posParts[1]);
                    groupByPos = new Range(posParts[0], posChromStart, posChromStart + 1);
                }
            }
        }
        return groupByPos;
    }

    public boolean isInvertSorting() {
        return invertSorting == null ? getPreferences().getAsBoolean(SAM_INVERT_SORT) : invertSorting;
    }

    public boolean isInvertGroupSorting() {
        return invertGroupSorting;
    }

    public String getLinkByTag() {
        return linkByTag == null ? getPreferences().get(SAM_LINK_TAG) : linkByTag;
    }

    public AlignmentTrack.GroupOption getGroupByOption() {
        AlignmentTrack.GroupOption gbo = groupByOption;
        // Interpret null as the default option.
        gbo = (gbo == null) ?
                CollUtils.valueOf(AlignmentTrack.GroupOption.class, getPreferences().get(SAM_GROUP_OPTION), AlignmentTrack.GroupOption.NONE) :
                gbo;
        // Add a second check for null in case defaultValues.groupByOption == null
        gbo = (gbo == null) ? AlignmentTrack.GroupOption.NONE : gbo;

        return gbo;
    }

    public boolean isLinkedReads() {
        return linkedReads == null ? getPreferences().getAsBoolean(SAM_LINKED_READS) : linkedReads;
    }

    public boolean isQuickConsensusMode() {
        return quickConsensusMode == null ? getPreferences().getAsBoolean(SAM_QUICK_CONSENSUS_MODE) : quickConsensusMode;
    }

    public boolean isHideSmallIndels() {
        return hideSmallIndels == null ? getPreferences().getAsBoolean(SAM_HIDE_SMALL_INDEL) : hideSmallIndels;
    }

    public int getSmallIndelThreshold() {
        return smallIndelThreshold == null ? getPreferences().getAsInt(SAM_SMALL_INDEL_BP_THRESHOLD) : smallIndelThreshold;
    }

    public BaseModficationFilter getBasemodFilter() {
        return basemodFilter;
    }

    public void setBasemodFilter(BaseModficationFilter basemodFilter) {

        this.basemodFilter = basemodFilter;
    }

    public float getBasemodThreshold() {
        return basemodThreshold == null ? getPreferences().getAsFloat(BASEMOD_THRESHOLD) : basemodThreshold.floatValue();
    }

    public void setBasemodThreshold(float basemodThreshold) {
        this.basemodThreshold = basemodThreshold;
    }

    public boolean getBasemodDistinguishStrands() {
        return basemodDistinguishStrands == null ? getPreferences().getAsBoolean(BASEMOD_DISTINGUISH_STRANDS) : basemodDistinguishStrands;
    }

    public void setBasemodDistinguishStrands(boolean basemodDistinguishStrands) {
        this.basemodDistinguishStrands = basemodDistinguishStrands;
    }

    public void unmarshalXML(Element element, Integer version) {
        if (element.hasAttribute("shadeBasesOption")) {
            String v = element.getAttribute("shadeBasesOption");
            if (v != null) {
                shadeBasesOption = v.equalsIgnoreCase("quality") || v.equalsIgnoreCase("true");
            }
        }
        if (element.hasAttribute("shadeCenters")) {
            shadeCenters = Boolean.parseBoolean(element.getAttribute("shadeCenters"));
        }
        if (element.hasAttribute("showAllBases")) {
            showAllBases = Boolean.parseBoolean(element.getAttribute("showAllBases"));
        }
        if (element.hasAttribute("flagUnmappedPairs")) {
            flagUnmappedPairs = Boolean.parseBoolean(element.getAttribute("flagUnmappedPairs"));
        }

        if (element.hasAttribute("minInsertSize")) {
            minInsertSize = Integer.parseInt(element.getAttribute("minInsertSize"));
        }
        if (element.hasAttribute("maxInsertSize")) {
            maxInsertSize = Integer.parseInt(element.getAttribute("maxInsertSize"));
        }
        if (element.hasAttribute("colorOption")) {
            // Convert deprecated options
            final String attributeValue = element.getAttribute("colorOption");
            if ("BASE_MODIFICATION_6MA".equals(attributeValue)) {
                colorOption = AlignmentTrack.ColorOption.BASE_MODIFICATION;
                basemodFilter = new BaseModficationFilter("a");
            } else if ("BASE_MODIFICATION_5MC".equals(attributeValue)) {
                colorOption = AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR;
                // basemodFilter = new BaseModficationFilter(null, 'C');
            } else if ("BASE_MODIFICATION_C".equals(attributeValue)) {
                colorOption = AlignmentTrack.ColorOption.BASE_MODIFICATION;
                // basemodFilter = new BaseModficationFilter(null, 'C');
            } else {
                colorOption = AlignmentTrack.ColorOption.valueOf(attributeValue);
            }
        }
        if (element.hasAttribute("sortOption")) {
            sortOption = SortOption.fromString((element.getAttribute("sortOption")));
        }
        if (element.hasAttribute("groupByOption")) {
            String value = element.getAttribute("groupByOption");
            if (value.equals("HAPLOTYPE")) {
                value = "CLUSTER";  // Backward compatibility
            }
            groupByOption = AlignmentTrack.GroupOption.valueOf(value);
        }
        if (element.hasAttribute("shadeAlignmentsByOption")) {
            shadeAlignmentsOption = AlignmentTrack.ShadeAlignmentsOption.valueOf(element.getAttribute("shadeAlignmentsByOption"));
        }
        if (element.hasAttribute("duplicatesOption")) {
            duplicatesOption = CollUtils.valueOf(AlignmentTrack.DuplicatesOption.class, element.getAttribute("duplicatesOption"), null);
        }
        if (element.hasAttribute("mappingQualityLow")) {
            mappingQualityLow = Integer.parseInt(element.getAttribute("mappingQualityLow"));
        }
        if (element.hasAttribute("mappingQualityHigh")) {
            mappingQualityHigh = Integer.parseInt(element.getAttribute("mappingQualityHigh"));
        }
        if (element.hasAttribute("viewPairs")) {
            viewPairs = Boolean.parseBoolean(element.getAttribute("viewPairs"));
        }
        if (element.hasAttribute("colorByTag")) {
            colorByTag = element.getAttribute("colorByTag");
        }
        if (element.hasAttribute("groupByTag")) {
            groupByTag = element.getAttribute("groupByTag");
        }
        if (element.hasAttribute("sortByTag")) {
            sortByTag = element.getAttribute("sortByTag");
        }
        if (element.hasAttribute("linkByTag")) {
            linkByTag = element.getAttribute("linkByTag");
        }
        if (element.hasAttribute("linkedReads")) {
            linkedReads = Boolean.parseBoolean(element.getAttribute("linkedReads"));
        }
        if (element.hasAttribute("quickConsensusMode")) {
            quickConsensusMode = Boolean.parseBoolean(element.getAttribute("quickConsensusMode"));
        }
        if (element.hasAttribute("showMismatches")) {
            showMismatches = Boolean.parseBoolean(element.getAttribute("showMismatches"));
        }
        if (element.hasAttribute("computeIsizes")) {
            computeIsizes = Boolean.parseBoolean(element.getAttribute("computeIsizes"));
        }
        if (element.hasAttribute("minInsertSizePercentile")) {
            minInsertSizePercentile = Double.parseDouble(element.getAttribute("minInsertSizePercentile"));
        }
        if (element.hasAttribute("maxInsertSizePercentile")) {
            maxInsertSizePercentile = Double.parseDouble(element.getAttribute("maxInsertSizePercentile"));
        }
        if (element.hasAttribute("pairedArcView")) {
            pairedArcView = Boolean.parseBoolean(element.getAttribute("pairedArcView"));
        }
        if (element.hasAttribute("flagZeroQualityAlignments")) {
            flagZeroQualityAlignments = Boolean.parseBoolean(element.getAttribute("flagZeroQualityAlignments"));
        }
        if (element.hasAttribute("groupByPos")) {
            groupByPos = Range.fromString(element.getAttribute("groupByPos"));
        }
        if (element.hasAttribute("invertSorting")) {
            invertSorting = Boolean.parseBoolean(element.getAttribute("invertSorting"));
        }
        if (element.hasAttribute("invertGroupSorting")) {
            invertGroupSorting = Boolean.parseBoolean(element.getAttribute("invertGroupSorting"));
        }
        if (element.hasAttribute("hideSmallIndels")) {
            hideSmallIndels = Boolean.parseBoolean(element.getAttribute("hideSmallIndels"));
        }
        if (element.hasAttribute("smallIndelThreshold")) {
            smallIndelThreshold = Integer.parseInt(element.getAttribute("smallIndelThreshold"));
        }
        if (element.hasAttribute("showInsertionMarkers")) {
            // TODO -- something with this
            // showInsertionMarkers = Boolean.parseBoolean(element.getAttribute("showInsertionMarkers"));
        }
        if (element.hasAttribute("basemodFilter")) {
            basemodFilter = BaseModficationFilter.fromString(element.getAttribute("basemodFilter"));
        }
        if (element.hasAttribute("basemodThreshold")) {
            basemodThreshold = Float.parseFloat(element.getAttribute("basemodThreshold"));
        }
        if (element.hasAttribute("minJunctionCoverage")) {
            minJunctionCoverage = Integer.parseInt(element.getAttribute("minJunctionCoverage"));
        }
    }

    public void marshalJSON(JSONObject jsonObject) {
        if (shadeBasesOption != null) {
            jsonObject.put("shadeBasesOption", shadeBasesOption);
        }
        if (shadeCenters != null) {
            jsonObject.put("shadeCenters", shadeCenters);
        }
        if (flagUnmappedPairs != null) {
            jsonObject.put("flagUnmappedPairs", flagUnmappedPairs);
        }
        if (showAllBases != null) {
            jsonObject.put("showAllBases", showAllBases);
        }
        if (minInsertSize != null) {
            jsonObject.put("minTLEN", minInsertSize);
        }
        if (maxInsertSize != null) {
            jsonObject.put("maxTLEN", maxInsertSize);
        }
        if (colorOption != null) {
            jsonObject.put("colorOption", colorOption.toString());
        }
        if (groupByOption != null) {
            jsonObject.put("groupByOption", groupByOption.toString());
        }
        if (shadeAlignmentsOption != null) {
            jsonObject.put("shadeAlignmentsByOption", shadeAlignmentsOption.toString());
        }
        if (duplicatesOption != null) {
            jsonObject.put("duplicatesOption", duplicatesOption.toString());
        }
        if (mappingQualityLow != null) {
            jsonObject.put("mappingQualityLow", mappingQualityLow);
        }
        if (mappingQualityHigh != null) {
            jsonObject.put("mappingQualityHigh", mappingQualityHigh);
        }
        if (viewPairs) {
            jsonObject.put("viewPairs", viewPairs);
        }
        if (colorByTag != null) {
            jsonObject.put("colorByTag", colorByTag);
        }
        if (groupByTag != null) {
            jsonObject.put("groupByTag", groupByTag);
        }
        if (sortByTag != null) {
            jsonObject.put("sortByTag", sortByTag);
        }
        if (linkByTag != null) {
            jsonObject.put("linkByTag", linkByTag);
        }
        if (linkedReads != null) {
            jsonObject.put("linkedReads", linkedReads);
        }
        if (quickConsensusMode != null) {
            jsonObject.put("quickConsensusMode", quickConsensusMode);
        }
        if (showMismatches != null) {
            jsonObject.put("showMismatches", showMismatches);
        }
        if (computeIsizes != null) {
            jsonObject.put("computeIsizes", computeIsizes);
        }
        if (minInsertSizePercentile != null) {
            jsonObject.put("minTLENPercentile", minInsertSizePercentile);
        }
        if (maxInsertSizePercentile != null) {
            jsonObject.put("maxTLENPercentile", maxInsertSizePercentile);
        }
        if (pairedArcView != null) {
            jsonObject.put("pairedArcView", pairedArcView);
        }
        if (flagZeroQualityAlignments != null) {
            jsonObject.put("flagZeroQualityAlignments", flagZeroQualityAlignments);
        }
        if (groupByPos != null) {
            jsonObject.put("groupByPos", groupByPos.toString());
        }
        if (invertSorting != null) {
            jsonObject.put("invertSorting", invertSorting);
        }
        if (sortOption != null) {
            jsonObject.put("sortOption", sortOption.toString());
        }
        if (invertGroupSorting) {
            jsonObject.put("invertGroupSorting", invertGroupSorting);
        }
        if (hideSmallIndels != null) {
            jsonObject.put("hideSmallIndels", hideSmallIndels);
        }
        if (smallIndelThreshold != null) {
            jsonObject.put("indexlSizeThreshold", smallIndelThreshold);
        }
        if (basemodFilter != null) {
            jsonObject.put("basemodFilter", basemodFilter.toString());
        }
        if (basemodThreshold != null) {
            jsonObject.put("baseModificationThreshold", basemodThreshold);
        }
        if (minJunctionCoverage != null) {
            jsonObject.put("minJunctionCoverage", minJunctionCoverage);
        }
    }

    public void unmarshalJSON(JSONObject json) {
        if (json.has("shadeBasesOption")) {
            String v = json.getString("shadeBasesOption");
            if (v != null) {
                shadeBasesOption = v.equalsIgnoreCase("quality") || v.equalsIgnoreCase("true");
            }
        }
        if (json.has("shadeCenters")) {
            shadeCenters = json.getBoolean("shadeCenters");
        }
        if (json.has("showAllBases")) {
            showAllBases = json.getBoolean("showAllBases");
        }
        if (json.has("flagUnmappedPairs")) {
            flagUnmappedPairs = json.getBoolean("flagUnmappedPairs");
        }

        if (json.has("minTLEN")) {
            minInsertSize = json.getInt("minTLEN");
        }
        if (json.has("maxTLEN")) {
            maxInsertSize = json.getInt("maxTLEN");
        }
        if (json.has("colorOption")) {
            // Convert deprecated options
            final String attributeValue = json.getString("colorOption");
            if ("BASE_MODIFICATION_6MA".equals(attributeValue)) {
                colorOption = AlignmentTrack.ColorOption.BASE_MODIFICATION;
                basemodFilter = new BaseModficationFilter("a");
            } else if ("BASE_MODIFICATION_5MC".equals(attributeValue)) {
                colorOption = AlignmentTrack.ColorOption.BASE_MODIFICATION_2COLOR;
                // basemodFilter = new BaseModficationFilter(null, 'C');
            } else if ("BASE_MODIFICATION_C".equals(attributeValue)) {
                colorOption = AlignmentTrack.ColorOption.BASE_MODIFICATION;
                // basemodFilter = new BaseModficationFilter(null, 'C');
            } else {
                colorOption = AlignmentTrack.ColorOption.valueOf(attributeValue);
            }
        }
        if (json.has("sortOption")) {
            sortOption = SortOption.fromString((json.getString("sortOption")));
        }
        if (json.has("groupByOption")) {
            String value = json.getString("groupByOption");
            if (value.equals("HAPLOTYPE")) {
                value = "CLUSTER";  // Backward compatibility
            }
            groupByOption = AlignmentTrack.GroupOption.valueOf(value);
        }
        if (json.has("shadeAlignmentsByOption")) {
            shadeAlignmentsOption = AlignmentTrack.ShadeAlignmentsOption.valueOf(json.getString("shadeAlignmentsByOption"));
        }
        if (json.has("duplicatesOption")) {
            duplicatesOption = CollUtils.valueOf(AlignmentTrack.DuplicatesOption.class, json.getString("duplicatesOption"), null);
        }
        if (json.has("mappingQualityLow")) {
            mappingQualityLow = json.getInt("mappingQualityLow");
        }
        if (json.has("mappingQualityHigh")) {
            mappingQualityHigh = json.getInt("mappingQualityHigh");
        }
        if (json.has("viewPairs")) {
            viewPairs = json.getBoolean("viewPairs");
        }
        if (json.has("colorByTag")) {
            colorByTag = json.getString("colorByTag");
        }
        if (json.has("groupByTag")) {
            groupByTag = json.getString("groupByTag");
        }
        if (json.has("sortByTag")) {
            sortByTag = json.getString("sortByTag");
        }
        if (json.has("linkByTag")) {
            linkByTag = json.getString("linkByTag");
        }
        if (json.has("linkedReads")) {
            linkedReads = json.getBoolean("linkedReads");
        }
        if (json.has("quickConsensusMode")) {
            quickConsensusMode = json.getBoolean("quickConsensusMode");
        }
        if (json.has("showMismatches")) {
            showMismatches = json.getBoolean("showMismatches");
        }
        if (json.has("computeIsizes")) {
            computeIsizes = json.getBoolean("computeIsizes");
        }
        if (json.has("minTLENPercentile")) {
            minInsertSizePercentile = Double.parseDouble(json.getString("minTLENPercentile"));
        }
        if (json.has("maxTLENPercentile")) {
            maxInsertSizePercentile = Double.parseDouble(json.getString("maxTLENPercentile"));
        }
        if (json.has("pairedArcView")) {
            pairedArcView = json.getBoolean("pairedArcView");
        }
        if (json.has("flagZeroQualityAlignments")) {
            flagZeroQualityAlignments = json.getBoolean("flagZeroQualityAlignments");
        }
        if (json.has("groupByPos")) {
            groupByPos = Range.fromString(json.getString("groupByPos"));
        }
        if (json.has("invertSorting")) {
            invertSorting = json.getBoolean("invertSorting");
        }
        if (json.has("invertGroupSorting")) {
            invertGroupSorting = json.getBoolean("invertGroupSorting");
        }
        if (json.has("hideSmallIndels")) {
            hideSmallIndels = json.getBoolean("hideSmallIndels");
        }
        if (json.has("indexlSizeThreshold")) {
            smallIndelThreshold = json.getInt("indexlSizeThreshold");
        }
        if (json.has("showInsertionMarkers")) {
            // TODO -- something with this
            // showInsertionMarkers = json.getBoolean("showInsertionMarkers"));
        }
        if (json.has("basemodFilter")) {
            basemodFilter = BaseModficationFilter.fromString(json.getString("basemodFilter"));
        }
        if (json.has("baseModificationThreshold")) {
            basemodThreshold = json.getFloat("basemodThreshold");
        }
        if (json.has("minJunctionCoverage")) {
            minJunctionCoverage = json.getInt("minJunctionCoverage");
        }
    }

}
