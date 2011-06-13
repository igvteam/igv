/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */


//chr2:128,565,093-128,565,156

package org.broad.igv.vcf;

import org.apache.log4j.Logger;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.renderer.*;
import org.broad.igv.session.SessionReader;
import org.broad.igv.track.*;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.ResourceLocator;
import org.broad.tribble.Feature;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;

import javax.swing.*;
import java.awt.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;

/**
 * User: Jesse Whitworth
 * Date: Jul 16, 2010
 */
public class VCFTrack extends FeatureTrack {

    private static Logger log = Logger.getLogger(VCFTrack.class);

    private static final Color OFF_WHITE = new Color(170, 170, 170);
    private static final int GROUP_BORDER_WIDTH = 3;
    private static final Color BAND1_COLOR = new Color(245, 245, 245);
    private static final Color BAND2_COLOR = Color.white;

    // Map for organizing samples by family (sample -> family).  This is static (shared) by all vcf tracks
    // We need maps in both directions to (1) look up a group quickly,  and (2) maintain proper order in the group
    private static Map<String, String> REFERENCE_SAMPLE_GROUP_MAP = new HashMap();
    private static LinkedHashMap<String, List<String>> REFERENCE_GROUP_SAMPLE_MAP;

    /**
     * When this flag is true, we have detected that the VCF file contains the FORMAT MR column representing
     * methylation data. This will enable the "Color By/Methylation Rate" menu item.
     */
    private boolean enableMethylationRateSupport;
    private static final Color borderGray = new Color(200, 200, 200);

    public boolean isEnableMethylationRateSupport() {

        return enableMethylationRateSupport;
    }

    public static void addSampleGroups(LinkedHashMap<String, List<String>> map) {
        REFERENCE_GROUP_SAMPLE_MAP = map;
        for (Map.Entry<String, List<String>> entry : map.entrySet()) {
            String family = entry.getKey();
            for (String sample : entry.getValue()) {
                REFERENCE_SAMPLE_GROUP_MAP.put(sample, family);
            }
        }

    }

    private final static int DEFAULT_EXPANDED_GENOTYPE_HEIGHT = 15;
    private final int DEFAULT_SQUISHED_GENOTYPE_HEIGHT = 4;
    private final static int DEFAULT_VARIANT_BAND_HEIGHT = 25;
    private final static int MAX_FILTER_LINES = 15;

    private VCFRenderer renderer = new VCFRenderer(this);

    private int squishedHeight = DEFAULT_SQUISHED_GENOTYPE_HEIGHT;

    // A hack, keeps track of last position drawn.  TODO -- need a proper component "model" for tracks, like a lightweight swing panel
    private int top;
    private int variantBandHeight = DEFAULT_VARIANT_BAND_HEIGHT;

    LinkedHashMap<String, List<String>> samplesByGroups = new LinkedHashMap();
    List<String> allSamples;
    List<String> groupNames;

    int sampleCount;
    private boolean grouped;
    private boolean hasGroups;

    private ColorMode coloring = ColorMode.GENOTYPE;
    private boolean hideAncestral = false;

    private boolean hideFiltered = false;
    private boolean renderID = true;

    private static float dash[] = {4.0f, 1.0f};
    DecimalFormat numFormat = new DecimalFormat("#.###");

    Feature selectedVariant;


    public VCFTrack(ResourceLocator locator, TribbleFeatureSource source) {
        super(locator, source);
        VCFHeader header = (VCFHeader) source.getHeader();

        // Test if the input VCF file contains methylation rate data:

        // This is determined by testing for the presence of two sample format fields: MR and GB, used in the
        // rendering of methylation rate.
        // MR is the methylation rate on a scale of 0 to 100% and GB is the number of bases that pass
        // filter for the position. GB is needed to avoid displaying positions for which limited coverage
        // prevents reliable estimation of methylation rate.
        if (header.getFormatHeaderLine("MR") != null &&
                header.getFormatHeaderLine("GB") != null) {
            enableMethylationRateSupport = true;
            // also set the default color mode to Methylation rate:
            coloring = ColorMode.METHYLATION_RATE;
        }

        allSamples = new ArrayList(header.getGenotypeSamples());
        sampleCount = allSamples.size();

        for (String sample : allSamples) {
            String key = REFERENCE_SAMPLE_GROUP_MAP.get(sample);
            if (key == null) {
                key = "Other";
            }
            List<String> sampleList = samplesByGroups.get(key);
            if (sampleList == null) {
                sampleList = new ArrayList();
                samplesByGroups.put(key, sampleList);
            }
            sampleList.add(sample);
        }

        grouped = samplesByGroups.size() > 1;
        hasGroups = samplesByGroups.size() > 1;
        if (grouped && REFERENCE_GROUP_SAMPLE_MAP != null) {
            sortGroups();
        }

        setDisplayMode(DisplayMode.EXPANDED);
        setRenderID(false);

        // Estimate visibility window.   TODO -- set beta based on available memory
        int cnt = Math.max(1, sampleCount);
        int beta = 100000;
        int visWindow = Math.min(500000, (beta / cnt) * 1000);
        setVisibilityWindow(visWindow);
    }

    private void sortGroups() {

        // First sort the sample groups (keys)
        final Map<String, Integer> groupRank = new HashMap();
        int idx = 0;
        for (String family : REFERENCE_GROUP_SAMPLE_MAP.keySet()) {
            groupRank.put(family, idx++);
        }

        groupNames = new ArrayList(samplesByGroups.keySet());
        Collections.sort(groupNames, new Comparator<String>() {
            public int compare(String s1, String s2) {
                int r1 = groupRank.containsKey(s1) ? groupRank.get(s1) : Integer.MAX_VALUE;
                int r2 = groupRank.containsKey(s2) ? groupRank.get(s2) : Integer.MAX_VALUE;
                return r1 - r2;
            }
        });


        LinkedHashMap<String, List<String>> newSamples = new LinkedHashMap();
        for (String family : groupNames) {
            List<String> sampleList = samplesByGroups.get(family);
            if (sampleList != null && sampleList.size() > 0) {
                final List<String> referenceGroup = REFERENCE_GROUP_SAMPLE_MAP.get(family);
                if (referenceGroup != null) {
                    final Map<String, Integer> rank = new HashMap();
                    for (int i = 0; i < referenceGroup.size(); i++) {
                        rank.put(referenceGroup.get(i), i);
                        Collections.sort(sampleList, new Comparator<String>() {
                            public int compare(String s1, String s2) {
                                int r1 = rank.containsKey(s1) ? rank.get(s1) : Integer.MAX_VALUE;
                                int r2 = rank.containsKey(s2) ? rank.get(s2) : Integer.MAX_VALUE;
                                return r1 - r2;
                            }
                        });
                    }
                }
                newSamples.put(family, sampleList);
            }

        }
        this.samplesByGroups = newSamples;
    }


    public int getGenotypeBandHeight() {
        switch (getDisplayMode()) {
            case SQUISHED:
                return getSquishedHeight();
            case COLLAPSED:
                return 0;
            default:
                return DEFAULT_EXPANDED_GENOTYPE_HEIGHT;

        }
    }

    public int getHeight() {

        final int groupCount = samplesByGroups.size();
        if (getDisplayMode() == Track.DisplayMode.COLLAPSED || groupCount == 0) {
            return variantBandHeight;
        } else {
            int margins = (groupCount - 1) * 3;
            return variantBandHeight + margins + (this.sampleCount * getGenotypeBandHeight());
        }

    }

    @Override
    public int getPreferredHeight() {
        return getHeight();
    }

    public void setHeight(int height) {

        final DisplayMode displayMode = getDisplayMode();

        // If collapsed there's nothing we can do to affect height
        if (displayMode == DisplayMode.COLLAPSED) {
            return;
        }

        // If height is < expanded height try "squishing" track, otherwise expand it
        final int groupCount = samplesByGroups.size();
        final int margins = (groupCount - 1) * 3;
        final int expandedHeight = variantBandHeight + margins + (this.sampleCount * getGenotypeBandHeight());
        if (height < expandedHeight) {
            setDisplayMode(DisplayMode.SQUISHED);
            squishedHeight = Math.max(1, (height - variantBandHeight - margins) / sampleCount);
        } else {
            if (displayMode != DisplayMode.EXPANDED) {
                setDisplayMode(DisplayMode.EXPANDED);
            }
        }


    }

    public void setDisplayMode(DisplayMode mode) {
        super.setDisplayMode(mode);
    }

    @Override
    protected void renderFeatureImpl(RenderContext context, Rectangle trackRectangle, PackedFeatures packedFeatures) {

        Graphics2D g2D = context.getGraphics();

        top = trackRectangle.y;
        final int left = trackRectangle.x;
        final int right = (int) trackRectangle.getMaxX();

        Rectangle visibleRectangle = context.getVisibleRect();

        // A disposable rect -- note this gets modified all over the place, bad practice
        Rectangle rect = new Rectangle(trackRectangle);
        rect.height = getGenotypeBandHeight();
        rect.y = trackRectangle.y + variantBandHeight;
        colorBackground(g2D, rect, visibleRectangle, false);

        if (top > visibleRectangle.y && top < visibleRectangle.getMaxY()) {
            g2D.drawLine(left, top + 1, right, top + 1);
        }

        List<Feature> features = packedFeatures.getFeatures();
        if (features.size() > 0) {

            final double locScale = context.getScale();
            final double origin = context.getOrigin();
            //byte[] reference;
            //int windowStart;
            int lastPX = -1;
            final double pXMin = rect.getMinX();
            final double pXMax = rect.getMaxX();

            for (Feature feature : features) {

                VariantContext variant = (VariantContext) feature;
                //char ref = getReference(variant, windowStart, reference);

                if (hideFiltered && variant.isFiltered()) {
                    continue;
                }

                // 1 -> 0 based coordinates
                int start = variant.getStart() - 1;
                int end = variant.getEnd();
                int pX = (int) ((start - origin) / locScale);
                int dX = (int) Math.max(2, (end - start) / locScale);

                if (pX + dX < pXMin) {
                    continue;
                }
                if (pX > pXMax) {
                    break;
                }
                int w = dX;
                int x = pX;
                if (w < 3) {
                    w = 3;
                    x--;
                }


                if (pX + dX > lastPX) {

                    rect.y = top;
                    rect.height = variantBandHeight;
                    if (rect.intersects(visibleRectangle)) {
                        renderer.renderVariantBand(variant, rect, x, w, context, hideFiltered);
                    }

                    if (getDisplayMode() != Track.DisplayMode.COLLAPSED) {
                        rect.y += rect.height;
                        rect.height = getGenotypeBandHeight();

                        // Loop through groups
                        if (grouped) {
                            for (Map.Entry<String, List<String>> entry : samplesByGroups.entrySet()) {
                                for (String sample : entry.getValue()) {
                                    if (rect.intersects(visibleRectangle)) {
                                        renderer.renderGenotypeBandSNP(variant, context, rect, x, w, sample, coloring,
                                                hideFiltered);
                                    }
                                    rect.y += rect.height;
                                }
                                g2D.setColor(OFF_WHITE);
                                g2D.fillRect(rect.x, rect.y, rect.width, GROUP_BORDER_WIDTH);
                                rect.y += GROUP_BORDER_WIDTH;
                            }
                        } else {
                            for (String sample : allSamples) {
                                if (rect.intersects(visibleRectangle)) {
                                    renderer.renderGenotypeBandSNP(variant, context, rect, x, w, sample, coloring,
                                            hideFiltered);
                                }
                                rect.y += rect.height;
                            }

                        }
                    }


                    boolean isSelected = selectedVariant != null && selectedVariant == variant;
                    if (isSelected) {
                        Graphics2D selectionGraphics = context.getGraphic2DForColor(Color.black);
                        selectionGraphics.drawRect(x, top, w, getHeight());
                    }

                    lastPX = pX + dX;

                }

            }
        } else {
            rect.height = variantBandHeight;
            rect.y = trackRectangle.y;
            g2D.setColor(Color.gray);
            GraphicUtils.drawCenteredText("No Variants Found", trackRectangle, g2D);
        }

        // Variant band border
        if (allSamples.size() > 0) {
            int variantBandY = trackRectangle.y + variantBandHeight;
            if (variantBandY >= visibleRectangle.y && variantBandY <= visibleRectangle.getMaxY()) {

                Graphics2D borderGraphics = context.getGraphic2DForColor(Color.black);
                borderGraphics.drawLine(left, variantBandY, right, variantBandY);
            }
        }

        // Bottom border
        int bottomY = trackRectangle.y + trackRectangle.height;
        if (bottomY >= visibleRectangle.y && bottomY <= visibleRectangle.getMaxY()) {
            g2D.drawLine(left, bottomY, right, bottomY);
        }
    }


    /**
     * Render the names
     * @param g2D
     * @param trackRectangle
     * @param visibleRectangle
     */
    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle visibleRectangle) {


        // A hack
        int top = trackRectangle.y;
        final int left = trackRectangle.x;
        final int right = (int) trackRectangle.getMaxX();

        Rectangle rect = new Rectangle(trackRectangle);
        //g2D.clearRect(rect.x, rect.y, rect.width, rect.height);
        g2D.setFont(FontManager.getFont(fontSize));
        g2D.setColor(BAND2_COLOR);

        if (top > visibleRectangle.y && top < visibleRectangle.getMaxY()) {
            g2D.drawLine(left, top + 1, right, top + 1);
        }


        g2D.setColor(Color.black);
        rect.height = variantBandHeight;
        if (rect.intersects(visibleRectangle)) {
            GraphicUtils.drawWrappedText(getName(), rect, g2D, false);
        }

        if (grouped) {
            g2D.drawLine(rect.x, rect.y + rect.height - 1, rect.x + rect.width, rect.y + rect.height - 1);
        }

        rect.y += rect.height;
        rect.height = getGenotypeBandHeight();
        if (getDisplayMode() != Track.DisplayMode.COLLAPSED) {
            colorBackground(g2D, rect, visibleRectangle, true);
        }

        // Bottom border
        int bottomY = trackRectangle.y + trackRectangle.height;
        if (bottomY >= visibleRectangle.y && bottomY <= visibleRectangle.getMaxY()) {
            g2D.setColor(borderGray);
            g2D.drawLine(left, bottomY, right, bottomY);
        }

        if (allSamples.size() > 0) {
            int variantBandY = trackRectangle.y + variantBandHeight;
            if (variantBandY >= visibleRectangle.y && variantBandY <= visibleRectangle.getMaxY()) {
                g2D.setColor(Color.black);
                g2D.drawLine(left, variantBandY, right, variantBandY);
            }
        }

    }


    private void colorBackground(Graphics2D g2D, Rectangle bandRectangle, Rectangle visibleRectangle, boolean renderNames) {


        if (getDisplayMode() == Track.DisplayMode.SQUISHED && squishedHeight < 4) {
            return;
        }

        boolean coloredLast = true;

        Rectangle textRectangle = new Rectangle(bandRectangle);
        textRectangle.height--;

        int bandFontSize = Math.min(fontSize, (int) bandRectangle.getHeight() - 1);
        Font font = FontManager.getFont(bandFontSize);
        Font oldFont = g2D.getFont();
        g2D.setFont(font);

        if (grouped) {
            for (List<String> sampleList : samplesByGroups.values()) {
                for (String sample : sampleList) {

                    if (coloredLast) {
                        g2D.setColor(BAND1_COLOR);
                        coloredLast = false;
                    } else {
                        g2D.setColor(BAND2_COLOR);
                        coloredLast = true;

                    }

                    if (bandRectangle.intersects(visibleRectangle)) {
                        g2D.fillRect(bandRectangle.x, bandRectangle.y, bandRectangle.width, bandRectangle.height);
                        if (renderNames && bandRectangle.height >= 3) {
                            String printName = sample;
                            textRectangle.y = bandRectangle.y + 1;
                            g2D.setColor(Color.black);
                            GraphicUtils.drawWrappedText(printName, bandRectangle, g2D, false);

                        }
                    }
                    bandRectangle.y += bandRectangle.height;

                }

                g2D.setColor(OFF_WHITE);
                g2D.fillRect(bandRectangle.x, bandRectangle.y, bandRectangle.width, GROUP_BORDER_WIDTH);
                bandRectangle.y += GROUP_BORDER_WIDTH;

            }
        } else {
            for (String sample : allSamples) {

                if (coloredLast) {
                    g2D.setColor(BAND1_COLOR);
                    coloredLast = false;
                } else {
                    g2D.setColor(BAND2_COLOR);
                    coloredLast = true;

                }

                if (bandRectangle.intersects(visibleRectangle)) {
                    g2D.fillRect(bandRectangle.x, bandRectangle.y, bandRectangle.width, bandRectangle.height);
                    if (renderNames && bandRectangle.height >= 3) {
                        String printName = sample;
                        textRectangle.y = bandRectangle.y + 1;
                        g2D.setColor(Color.black);
                        g2D.setFont(font);
                        GraphicUtils.drawWrappedText(printName, bandRectangle, g2D, false);

                    }
                }
                bandRectangle.y += bandRectangle.height;

            }

        }
        g2D.setFont(oldFont);
    }

    public void setRenderID(boolean value) {
        this.renderID = value;
    }

    public boolean getRenderID() {
        return renderID;
    }

    public void setHideAncestral(boolean value) {
        this.hideAncestral = value;
    }

    public boolean getHideFiltered() {
        return hideFiltered;
    }

    public void setHideFiltered(boolean value) {
        this.hideFiltered = value;
    }

    public boolean getHideAncestral() {
        return hideAncestral;
    }

    public ColorMode getColorMode() {
        return coloring;
    }

    public void setColorMode(ColorMode mode) {
        this.coloring = mode;
    }


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        VariantContext variant = (VariantContext) getFeatureAt(chr, position + 1, y, frame); //getVariantAtPosition(chr, (int) position, frame);
        if (variant != null) {

            if (y < top + variantBandHeight) {
                return getVariantToolTip(variant);
            } else {
                String sample = null;
                if (grouped) {
                    // TODO This is a hack for the autism site, fix SOON.  Assumes groups are all trios
                    int groupHeight = 3 * getGenotypeBandHeight() + 3;
                    int groupNumber = Math.min(groupNames.size() - 1, (y - top - variantBandHeight) / groupHeight);
                    String group = groupNames.get(groupNumber);
                    List<String> sampleList = samplesByGroups.get(group);
                    int sampleNumber = (y - top - variantBandHeight - groupNumber * groupHeight) / getGenotypeBandHeight();
                    if (sampleNumber >= 0 || sampleNumber < sampleList.size()) {
                        sample = sampleList.get(sampleNumber);
                    }

                } else {
                    int sampleNumber = (y - top - variantBandHeight) / getGenotypeBandHeight();
                    if (sampleNumber >= 0 && sampleNumber < allSamples.size()) {
                        sample = allSamples.get(sampleNumber);
                    }

                }
                return getSampleToolTip(sample, variant);
            }
        }
        return null;
    }


    private String getVariantToolTip(VariantContext variant) {
        String id = variant.getAttributeAsString(VariantContext.ID_KEY);

        StringBuffer toolTip = new StringBuffer();
        toolTip.append("Chr:" + variant.getChr());
        toolTip.append("<br>Position:" + variant.getStart());
        toolTip.append("<br>ID: " + id);
        toolTip.append("<br>Reference: " + variant.getReference().toString());
        Set alternates = variant.getAlternateAlleles();
        if (alternates.size() > 0) {
            toolTip.append("<br>Alternate: " + alternates.toString());
        }

        toolTip.append("<br>Qual: " + numFormat.format(variant.getPhredScaledQual()));
        toolTip.append("<br>Type: " + variant.getType());
        if (variant.isFiltered()) {
            toolTip.append("<br>Is Filtered Out: Yes</b>");
            toolTip = toolTip.append(getFilterTooltip(variant));
        } else {
            toolTip.append("<br>Is Filtered Out: No</b><br>");
        }
        toolTip.append("<br><b>Alleles:</b>");
        toolTip.append(getAlleleToolTip(getZygosityCounts(variant)));

        double af = getAlleleFreq(variant);
        if (af < 0 && variant.getSampleNames().size() > 0) {
            af = getAlleleFraction(variant);
        }
        toolTip.append("<br>Allele Frequency: " + (af >= 0 ? numFormat.format(af) : "Unknown") + "<br>");

        if (variant.getSampleNames().size() > 0) {
            double afrac = getAlleleFraction(variant);
            toolTip = toolTip.append("<br>Minor Allele Fraction: " + numFormat.format(afrac) + "<br>");
        }

        toolTip.append("<br><b>Genotypes:</b>");
        toolTip.append(getGenotypeToolTip(getZygosityCounts(variant)) + "<br>");
        toolTip.append(getVariantInfo(variant) + "<br>");
        return toolTip.toString();
    }

    protected String getVariantInfo(VariantContext variant) {
        Set<String> keys = variant.getAttributes().keySet();
        if (keys.size() > 0) {
            String toolTip = "<br><b>Variant Attributes</b>";
            int count = 0;

            // Put AF and GMAF and put at the top, if present
            String k = "AF";
            String afValue = variant.getAttributeAsString(k);
            if (afValue != null && afValue.length() > 0 && !afValue.equals("null")) {
                toolTip = toolTip.concat("<br>" + getFullName(k) + ": " + variant.getAttributeAsString(k));
            }
            k = "GMAF";
            afValue = variant.getAttributeAsString(k);
            if (afValue != null && afValue.length() > 0 && !afValue.equals("null")) {
                toolTip = toolTip.concat("<br>" + getFullName(k) + ": " + variant.getAttributeAsString(k));
            }

            for (String key : keys) {
                count++;

                if (key.equals("AF") || key.equals("GMAF")) continue;

                if (count > MAX_FILTER_LINES) {
                    toolTip = toolTip.concat("<br>....");
                    break;
                }
                toolTip = toolTip.concat("<br>" + getFullName(key) + ": " + variant.getAttributeAsString(key));

            }
            return toolTip;
        }
        return " ";
    }

    private String getSampleInfo(Genotype genotype) {
        Set<String> keys = genotype.getAttributes().keySet();
        if (keys.size() > 0) {
            String tooltip = "<br><b>Sample Attributes</b>";
            for (String key : keys) {
                try {
                    tooltip = tooltip.concat("<br>" + getFullName(key) + ": " + genotype.getAttributeAsString(key));
                } catch (IllegalArgumentException iae) {
                    tooltip = tooltip.concat("<br>" + key + ": " + genotype.getAttributeAsString(key));
                }
            }
            return tooltip;
        }
        return null;
    }

    public void clearSelectedVariant() {
        selectedVariant = null;
    }

    public List<String> getAllSamples() {
        return allSamples;
    }

    public void setAllSamples(List<String> samples) {
        this.allSamples = samples;
    }

    public boolean isGrouped() {
        return grouped;
    }

    public void setGrouped(boolean grouped) {
        this.grouped = grouped;
    }

    public boolean isHasGroups() {
        return hasGroups;
    }

    public int getSquishedHeight() {
        return squishedHeight;
    }

    public void setSquishedHeight(int squishedHeight) {
        this.squishedHeight = squishedHeight;
    }

    public static enum ColorMode {
        GENOTYPE, METHYLATION_RATE, ALLELE
    }


    static Map<String, String> fullNames = new HashMap();

    static {
        fullNames.put("AA", "Ancestral Allele");
        fullNames.put("AC", "Allele Count in Genotypes");
        fullNames.put("AN", "Total Alleles in Genotypes");
        fullNames.put("AF", "Allele Frequency");
        fullNames.put("DP", "Depth");
        fullNames.put("MQ", "Mapping Quality");
        fullNames.put("NS", "Number of Samples with Data");
        fullNames.put("BQ", "RMS Base Quality");
        fullNames.put("SB", "Strand Bias");
        fullNames.put("DB", "dbSNP Membership");
        fullNames.put("GQ", "Genotype Quality");
        fullNames.put("GL", "Genotype Likelihoods");  //Hom-ref, het, hom-var break this down into a group
    }

    static String getFullName(String key) {
        return fullNames.containsKey(key) ? fullNames.get(key) : key;
    }


    private String getSampleToolTip(String sample, VariantContext variant) {
        String id = variant.getAttributeAsString(VariantContext.ID_KEY);
        StringBuffer toolTip = new StringBuffer();
        toolTip = toolTip.append("Chr:" + variant.getChr());
        toolTip = toolTip.append("<br>Position:" + variant.getStart());
        toolTip = toolTip.append("<br>ID: " + id + "<br>");
        toolTip = toolTip.append("<br><b>Sample Information</b>");
        toolTip = toolTip.append("<br>Sample: " + sample);
        toolTip = toolTip.append("<br>Position:" + variant.getStart());

        Genotype genotype = variant.getGenotype(sample);
        if (genotype != null) {
            toolTip = toolTip.append("<br>Bases: " + genotype.getGenotypeString());
            toolTip = toolTip.append("<br>Quality: " + numFormat.format(genotype.getPhredScaledQual()));
            toolTip = toolTip.append("<br>Type: " + genotype.getType());
        }
        if (variant.isFiltered()) {
            toolTip = toolTip.append("<br>Is Filtered Out: Yes</b>");
            toolTip = toolTip.append(getFilterTooltip(variant));
        } else {
            toolTip = toolTip.append("<br>Is Filtered Out: No</b><br>");
        }

        if (genotype != null) {
            toolTip = toolTip.append(getSampleInfo(genotype) + "<br>");
        }
        return toolTip.toString();
    }


    private String getFilterTooltip(VariantContext variant) {
        Set filters = variant.getFilters();
        String toolTip = "<br>";
        for (Object filter : filters) {
            toolTip = toolTip.concat("- " + (String) filter + "<br>");
        }

        return toolTip;
    }

    /**
     * Return the allele frequency as annotated with an AF or GMAF attribute.  A value of -1 indicates
     * no annotation (unknown allele frequency).
     */
    public double getAlleleFreq(VariantContext variant) {
        double alleleFreq = Double.parseDouble(variant.getAttributeAsString("AF", "-1"));
        if (alleleFreq < 0) {
            alleleFreq = Double.parseDouble(variant.getAttributeAsString("GMAF", "-1"));
        }
        return alleleFreq;
    }

    /**
     * Return the allele fraction for this variant.  The allele fraction is similiar to allele frequency, but is based
     * on the samples in this VCF as opposed to an AF or GMAF annotation.
     * <p/>
     * A value of -1 indicates unknown
     */
    public double getAlleleFraction(VariantContext variant) {
        ZygosityCount counts = getZygosityCounts(variant);
        int total = counts.getHomVar() + counts.getHet() + counts.getHomRef();
        return total == 0 ? -1 : (((double) counts.getHomVar() + ((double) counts.getHet()) / 2) / total);

    }

    public double getAllelePercent(VariantContext variant) {

        double af = getAlleleFraction(variant);
        return af < 0 ? getAlleleFreq(variant) : af;

    }

    static class ZygosityCount {
        private int homVar = 0;
        private int het = 0;
        private int homRef = 0;
        private int noCall = 0;

        public void incrementCount(Genotype genotype) {
            if (genotype != null) {
                if (genotype.isHomVar()) {
                    homVar++;
                } else if (genotype.isHet()) {
                    het++;
                } else if (genotype.isHomRef()) {
                    homRef++;
                } else {
                    noCall++;
                }
            }
        }

        public int getHomVar() {
            return homVar;
        }

        public int getHet() {
            return het;
        }

        public int getHomRef() {
            return homRef;
        }

        public int getNoCall() {
            return noCall;
        }

        public int getTotalCall() {
            return homVar + homRef + het;
        }

        public int getVarCall() {
            return homVar + het;
        }

        public int getSampleCount() {
            return homVar + homRef + het + noCall;
        }

    }

    public ZygosityCount getZygosityCounts(VariantContext variant) {
        ZygosityCount zc = new ZygosityCount();
        for (String sample : allSamples) {
            Genotype genotype = variant.getGenotype(sample);
            zc.incrementCount(genotype);
        }
        return zc;
    }

    private String getAlleleToolTip(ZygosityCount counts) {
        double noCall = counts.getNoCall() * 2;
        double aNum = (counts.getHet() + counts.getHomRef() + counts.getHomVar()) * 2;
        double aCount = (counts.getHomVar() * 2 + counts.getHet()) * 2;

        String toolTip = "<br>No Call: " + (int) noCall;
        toolTip = toolTip.concat("<br>Allele Num: " + (int) aNum);
        toolTip = toolTip.concat("<br>Allele Count: " + (int) aCount);
        return toolTip;
    }

    private String getGenotypeToolTip(ZygosityCount counts) {
        int noCall = counts.getNoCall();
        int homRef = counts.getHomRef();
        int nonVar = noCall + homRef;
        int het = counts.getHet();
        int homVar = counts.getHomVar();
        int var = het + homVar;

        String toolTip = "<br>Non Variant: " + nonVar;
        toolTip = toolTip.concat("<br> - No Call: " + noCall);
        toolTip = toolTip.concat("<br> - Hom Ref: " + homRef);
        toolTip = toolTip.concat("<br>Variant: " + var);
        toolTip = toolTip.concat("<br> - Het: " + het);
        toolTip = toolTip.concat("<br> - Hom Var: " + homVar);
        return toolTip;
    }


    public JPopupMenu getPopupMenu(final TrackClickEvent te) {
        VariantContext f = null;
        if (te.getFrame() != null && te.getFrame().getName() != null) {
            f = (VariantContext) getFeatureClosest(te.getChromosomePosition(), te.getMouseEvent().getY(), te.getFrame());
            selectedVariant = f;
            IGV.getInstance().doRefresh();
        }

        return new VCFMenu(this, f);
    }

    public static void refresh() {
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                IGV.getInstance().doRefresh();
            }
        });
    }

    public Map<String, String> getPersistentState() {

        Map<String, String> attributes = super.getPersistentState();
        attributes.put(SessionReader.SessionAttribute.RENDER_NAME.getText(), String.valueOf(renderID));

        ColorMode mode = getColorMode();
        if (mode != null) {
            attributes.put(SessionReader.SessionAttribute.COLOR_MODE.getText(), mode.toString());
        }

        if (squishedHeight != DEFAULT_SQUISHED_GENOTYPE_HEIGHT) {
            attributes.put("SQUISHED_ROW_HEIGHT", String.valueOf(squishedHeight));
        }

        return attributes;
    }


    public void restorePersistentState(Map<String, String> attributes) {
        super.restorePersistentState(attributes);

        String rendername = attributes.get(SessionReader.SessionAttribute.RENDER_NAME.getText());
        if (rendername != null) {
            setRenderID(rendername.equalsIgnoreCase("true"));
        }

        String colorModeText = attributes.get(SessionReader.SessionAttribute.COLOR_MODE.getText());
        if (colorModeText != null) {
            try {
                setColorMode(ColorMode.valueOf(colorModeText));
            }
            catch (Exception e) {
                log.error("Error interpreting display mode: " + colorModeText);
            }
        }

        String squishedHeightText = attributes.get("SQUISHED_ROW_HEIGHT");
        if (squishedHeightText != null) {
            try {
                squishedHeight = Integer.parseInt(squishedHeightText);
            }
            catch (Exception e) {
                log.error("Error restoring squished height: " + squishedHeightText);
            }
        }
    }
}
