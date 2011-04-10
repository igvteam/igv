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
import org.broad.igv.ui.UIConstants;
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

    private final int EXPANDED_GENOTYPE_HEIGHT = 15;
    private final int SQUISHED_GENOTYPE_HEIGHT = 4;
    private final int DEFAULT_VARIANT_BAND_HEIGHT = 25;
    private final int MAX_FILTER_LINES = 15;

    private VCFRenderer renderer = new VCFRenderer();

    // A hack, keeps track of last position drawn.  TODO -- need a proper component "model" for tracks, like a lightweight swing panel
    private int top;
    private int variantBandHeight = DEFAULT_VARIANT_BAND_HEIGHT;

    // Samples, organized by pedigree or other grouping
    private LinkedHashMap<String, List<String>> samples = new LinkedHashMap();
    int groupCount;
    int sampleCount;

    private ColorMode coloring = ColorMode.GENOTYPE;
    private boolean hideAncestral = false;


    private boolean hideFiltered = true;
    private boolean renderID = true;

    private static float dash[] = {4.0f, 1.0f};
    DecimalFormat numFormat = new DecimalFormat("#.###");

    Feature selectedVariant;


    public VCFTrack(ResourceLocator locator, TribbleFeatureSource source) {
        super(locator, source);
        VCFHeader header = (VCFHeader) source.getHeader();

        Set<String> allSamples = header.getGenotypeSamples();
        groupCount = 1;
        sampleCount = allSamples.size();
        samples.put("All", new ArrayList<String>(allSamples));


        this.setDisplayMode(DisplayMode.EXPANDED);
        setRenderID(false);

        // Estimate visibility window.   TODO -- set beta based on available memory
        int cnt = Math.max(1, sampleCount);
        int beta = 20000;
        int visWindow = Math.min(500000, (beta / cnt) * 1000);
        setVisibilityWindow(visWindow);
    }


    // TODO -- make this work with grouping

    public List<String> getAllSamples() {
        return samples.get("All");
    }

    public void setAllSamples(List<String> allSamples) {
        samples.put("All", allSamples);
    }


    public int getGenotypeBandHeight() {
        switch (getDisplayMode()) {
            case SQUISHED:
                return SQUISHED_GENOTYPE_HEIGHT;
            case COLLAPSED:
                return 0;
            default:
                return EXPANDED_GENOTYPE_HEIGHT;

        }
    }

    public int getHeight() {
        if (getDisplayMode() == Track.DisplayMode.COLLAPSED) {
            return variantBandHeight;
        } else {
            return variantBandHeight + (groupCount - 1) * 3 + (sampleCount * getGenotypeBandHeight());
        }

    }

    @Override
    public int getPreferredHeight() {
        return getHeight();
    }

    public void setHeight(int height) {
        this.height = Math.max(minimumHeight, height);
    }

    public void setDisplayMode(DisplayMode mode) {
        super.setDisplayMode(mode);
        setHeight(getPreferredHeight());
    }

    @Override
    protected void renderFeatureImpl(RenderContext context, Rectangle trackRectangle, PackedFeatures packedFeatures) {

        Graphics2D g2D = context.getGraphics();

        top = trackRectangle.y;
        Rectangle visibleRectangle = context.getVisibleRect();

        // A disposable rect -- note this gets modified all over the place, bad practice
        Rectangle rect = new Rectangle(trackRectangle);
        rect.height = getGenotypeBandHeight();
        rect.y = trackRectangle.y + variantBandHeight;
        colorBackground(g2D, rect, visibleRectangle, false, false);

        if (top > visibleRectangle.y && top < visibleRectangle.getMaxY()) {
            drawBorderLine(g2D, top + 1, trackRectangle.x, trackRectangle.x + trackRectangle.width);
        }

        List<Feature> features = packedFeatures.getFeatures();
        if (features.size() > 0) {

            double locScale = context.getScale();
            //byte[] reference;
            //int windowStart;
            int lastPX = -1;
            double pXEnd = rect.getMaxX();

            for (Feature feature : features) {

                VariantContext variant = (VariantContext) feature;
                //char ref = getReference(variant, windowStart, reference);

                // 1 -> 0 based coordinates
                int start = variant.getStart() - 1;
                int end = variant.getEnd();
                int pX = (int) ((start - context.getOrigin()) / locScale);
                int dX = (int) Math.max(2, (end - start) / locScale);

                if (pX > pXEnd) {
                    break;
                }
                int w = dX;
                int x = pX;
                if (w < 3) {
                    w = 3;
                    x--;
                }


                if (pX + dX > lastPX) {
                    ZygosityCount zygCounts = getZygosityCounts(variant);

                    rect.y = top;
                    rect.height = variantBandHeight;
                    if (rect.intersects(visibleRectangle)) {
                        if (getAllSamples().size() == 0) {
                            renderer.renderVariant(variant, rect, pX, dX, context);
                        } else {
                            AlleleCount alleleCounts = new AlleleCount(zygCounts);
                            renderer.renderAlleleBand(variant, rect, x, w, context, hideFiltered, alleleCounts);

                        }
                    }

                    if (this.getDisplayMode() != Track.DisplayMode.COLLAPSED) {
                        rect.y += rect.height;
                        rect.height = getGenotypeBandHeight();

                        // Loop through groups
                        for (Map.Entry<String, List<String>> entry : samples.entrySet()) {

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
                    }


                    boolean isSelected = selectedVariant != null && selectedVariant == variant;
                    if(isSelected) {
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

        // Bottom border
        int bottomY = trackRectangle.y + trackRectangle.height;
        if (bottomY >= visibleRectangle.y && bottomY <= visibleRectangle.getMaxY()) {
            final int left = trackRectangle.x;
            final int right = (int) trackRectangle.getMaxX();
            drawBorderLine(g2D, bottomY, left, right);
        }
    }


    private void drawBorderLine(Graphics2D g2D, int bottomY, int left, int right) {
        g2D.setColor(Color.black);
        g2D.drawLine(left, bottomY, right, bottomY);
    }


    private void colorBackground(Graphics2D g2D, Rectangle bandRectangle, Rectangle visibleRectangle, boolean renderNames,
                                 boolean isSelected) {
        boolean coloredLast = true;

        Rectangle textRectangle = new Rectangle(bandRectangle);
        textRectangle.height--;

        Font font = FontManager.getScalableFont((int) bandRectangle.getHeight() - 1);
        Font oldFont = g2D.getFont();
        g2D.setFont(font);
        for (List<String> sampleList : samples.values()) {
            for (String sample : sampleList) {
                if (isSelected) {
                    g2D.setColor(Color.lightGray);
                } else {
                    if (coloredLast) {
                        g2D.setColor(BAND1_COLOR);
                        coloredLast = false;
                    } else {
                        g2D.setColor(BAND2_COLOR);
                        coloredLast = true;
                    }
                }

                if (bandRectangle.intersects(visibleRectangle)) {
                    g2D.fillRect(bandRectangle.x, bandRectangle.y, bandRectangle.width, bandRectangle.height);
                    if (renderNames && bandRectangle.height >= 3) {
                        String printName = sample;
                        if (UIConstants.isSigmaProject()) {
                            if (printName.equals("384") || printName.equals("391")) {
                                printName = "S-" + printName;
                            } else if (printName.equals("467")) {
                                printName = "PA-" + printName;
                            } else if (printName.equals("469")) {
                                printName = "CC-" + printName;
                            } else if (printName.equals("259") || printName.equals("265") || printName.equals("266")) {
                                printName = "AA-" + printName;
                            } else if (printName.equals("701") || printName.equals("564")) {
                                printName = "BIP-" + printName;
                            } else if (printName.equals("352") || printName.equals("414") || printName.equals("375")) {
                                printName = "L-" + printName;
                            } else if (printName.equals("4") || printName.equals("8") || printName.equals("13") || printName.equals("15")) {
                                printName = "LA-" + printName;
                            } else if (printName.equals("563") || printName.equals("566")) {
                                printName = "OK-" + printName;
                            } else if (printName.equals("491") || printName.equals("497")) {
                                printName = "WC-" + printName;
                            }
                        }

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


    @Override
    public void renderName(Graphics2D g2D, Rectangle trackRectangle, Rectangle visibleRectangle) {


        // A hack
        int top = trackRectangle.y;

        Rectangle rect = new Rectangle(trackRectangle);
        g2D.clearRect(rect.x, rect.y, rect.width, rect.height);
        g2D.setFont(FontManager.getScalableFont(10));
        if (isSelected()) {
            g2D.setColor(Color.lightGray);
        } else {
            g2D.setColor(BAND2_COLOR);
        }

        if (top > visibleRectangle.y && top < visibleRectangle.getMaxY()) {
            drawBorderLine(g2D, top + 1, trackRectangle.x, trackRectangle.x + trackRectangle.width);
        }


        g2D.setColor(Color.black);
        rect.height = variantBandHeight;
        if (rect.intersects(visibleRectangle)) {
            GraphicUtils.drawWrappedText(getName(), rect, g2D, false); //getName() + " (" + samples.size() + ")", rect, g2D, false);
        }

        rect.y += rect.height;
        rect.height = getGenotypeBandHeight();
        if (getDisplayMode() != Track.DisplayMode.COLLAPSED) {
            // Sample names printed in colorBackground !!
            colorBackground(g2D, rect, visibleRectangle, true, isSelected());
        }

        // Bottom border
        int bottomY = trackRectangle.y + trackRectangle.height;
        if (bottomY >= visibleRectangle.y && bottomY <= visibleRectangle.getMaxY()) {
            final int left = trackRectangle.x;
            final int right = (int) trackRectangle.getMaxX();
            drawBorderLine(g2D, bottomY, left, right);
        }
    }


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {

        VariantContext variant = (VariantContext) getFeatureAt(chr, position, y, frame); //getVariantAtPosition(chr, (int) position, frame);
        if (variant != null) {

            if (y < top + variantBandHeight) {
                return getVariantToolTip(variant);
            } else {
                int sampleNumber = (y - top - variantBandHeight) / getGenotypeBandHeight();
                String sample = getAllSamples().get(sampleNumber);
                return getSampleToolTip(sample, variant);
            }
        }
        return null;
    }


    protected String getVariantInfo(VariantContext variant) {
        Set<String> keys = variant.getAttributes().keySet();
        if (keys.size() > 0) {
            String toolTip = "<br><b>Variant Attributes</b>";
            int count = 0;
            for (String key : keys) {
                try {
                    count++;
                    if (count > MAX_FILTER_LINES) {
                        toolTip = toolTip.concat("<br>....");
                        break;
                    }
                    toolTip = toolTip.concat("<br>" + InfoFieldName.findEnum(key) + ": " + variant.getAttributeAsString(key));
                } catch (IllegalArgumentException iae) {
                    toolTip = toolTip.concat("<br>" + key + ": " + variant.getAttributeAsString(key));
                }
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
                    tooltip = tooltip.concat("<br>" + InfoFieldName.findEnum(key) + ": " + genotype.getAttributeAsString(key));
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

    public static enum ColorMode {
        GENOTYPE, ALLELE
    }

    ;

    public static enum InfoFieldName {

        AA("Ancestral Allele"),
        AC("Allele Count in Genotypes"),
        AN("Total Alleles in Genotypes"),
        AF("Allele Frequency"),
        DP("Depth"),
        MQ("Mapping Quality"),
        NS("Number of Samples with Data"),
        BQ("RMS Base Quality"),
        SB("Strand Bias"),
        DB("dbSNP Membership"),
        GQ("Genotype Quality"),
        GL("Genotype Likelihoods");  //Hom-ref, het, hom-var break this down into a group

        private String name;

        InfoFieldName(String name) {
            this.name = name;
        }

        public String getText() {
            return name;
        }

        @Override
        public String toString() {
            return getText();
        }

        static public InfoFieldName findEnum(String value) {

            if (value == null) {
                return null;
            } else {
                return InfoFieldName.valueOf(value);
            }
        }
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

    private String getVariantToolTip(VariantContext variant) {
        String id = variant.getAttributeAsString(VariantContext.ID_KEY);
        StringBuffer toolTip = new StringBuffer();
        toolTip = toolTip.append("Chr:" + variant.getChr());
        toolTip = toolTip.append("<br>Position:" + variant.getStart());
        toolTip = toolTip.append("<br>ID: " + id);
        toolTip = toolTip.append("<br>Reference: " + variant.getReference().toString());
        Set alternates = variant.getAlternateAlleles();
        if (alternates.size() > 0) {
            toolTip = toolTip.append("<br>Alternate: " + alternates.toString());
        }

        toolTip = toolTip.append("<br>Qual: " + numFormat.format(variant.getPhredScaledQual()));
        toolTip = toolTip.append("<br>Type: " + variant.getType());
        if (variant.isFiltered()) {
            toolTip = toolTip.append("<br>Is Filtered Out: Yes</b>");
            toolTip = toolTip.append(getFilterTooltip(variant));
        } else {
            toolTip = toolTip.append("<br>Is Filtered Out: No</b><br>");
        }
        toolTip = toolTip.append("<br><b>Alleles:</b>");
        toolTip = toolTip.append(getAlleleToolTip(getZygosityCounts(variant)));
        toolTip = toolTip.append("<br>Allele Frequency: " + numFormat.format(getAlleleFreq(variant)) + "<br>");
        toolTip = toolTip.append("<br><b>Genotypes:</b>");
        toolTip = toolTip.append(getGenotypeToolTip(getZygosityCounts(variant)) + "<br>");
        toolTip = toolTip.append(getVariantInfo(variant) + "<br>");
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

    public double getAlleleFreq(VariantContext variant) {
        double alleleFreq = Double.valueOf(variant.getAttributeAsString("AF", "-1"));
        if (alleleFreq == -1) {
            ZygosityCount counts = getZygosityCounts(variant);
            int total = counts.getHomVar() + counts.getHet() + counts.getHomRef();
            return (((double) counts.getHomVar() + ((double) counts.getHet()) / 2) / total);
        }
        return alleleFreq;
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
        for (String sample : getAllSamples()) {
            Genotype genotype = variant.getGenotype(sample);
            zc.incrementCount(genotype);
        }
        return zc;
    }


    /*
   int noCall = counts[3];
   int homRef = counts[2];
   int nonVar = noCall + homRef;
   int het = counts[1];
   int homVar = counts[0];
   int var = het + homVar;

    */

    class AlleleCount {
        private int totalAlleles;
        private int alleleNum;
        private int alleleCount;


        public AlleleCount(ZygosityCount zygCounts) {
            totalAlleles = getAllSamples().size() * 2;
            alleleNum = (zygCounts.getHomVar() + zygCounts.getHet() + zygCounts.getHomRef()) * 2;
            alleleCount = zygCounts.getHomVar() * 2 + zygCounts.getHet();
        }

        public int getTotalAlleles() {
            return totalAlleles;
        }

        public int getAlleleNum() {
            return alleleNum;
        }

        public int getAlleleCount() {
            return alleleCount;
        }

        public float getAllelePercent() {
            return ((float) alleleCount / alleleNum);
        }
    }

    public AlleleCount getRenderCounts(VariantContext variant) {
        ZygosityCount zygCounts = getZygosityCounts(variant);
        AlleleCount alleleCounts = new AlleleCount(zygCounts);
        return alleleCounts;
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


    // public List<String> getSamples() {
    //     return samples;
    // }

    // protected void setSamples(List<String> samples) {
    //     this.samples = samples;
    // }

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
        return attributes;
    }


    public void restorePersistentState(Map<String, String> attributes) {
        super.restorePersistentState(attributes);

        String rendername = attributes.get(SessionReader.SessionAttribute.RENDER_NAME.getText());
        String colorModeText = attributes.get(SessionReader.SessionAttribute.COLOR_MODE.getText());


        // Set expand
        if (rendername != null) {
            setRenderID(rendername.equalsIgnoreCase("true"));
        }

        if (colorModeText != null) {
            try {
                setColorMode(ColorMode.valueOf(colorModeText));
            }
            catch (Exception e) {
                log.error("Error interpreting display mode: " + colorModeText);
            }
        }

        setHeight(getPreferredHeight());
    }
}
