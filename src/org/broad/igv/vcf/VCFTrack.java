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
import org.broad.igv.track.tribble.TribbleFeatureSource;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.IGVMainFrame;
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

    private final int EXPANDED_GENOTYPE_HEIGHT = 15;
    private final int SQUISHED_GENOTYPE_HEIGHT = 4;

    private final int MAX_GENOTYPE_HEIGHT = 50;
    private final int MIN_GENOTYPE_HEIGHT = 1;
    private final int DEFAULT_ALLELE_HEIGHT = 25;
    private final int MAX_FILTER_LINES = 15;

    private VCFRenderer renderer = new VCFRenderer();

    private VCFMenu menu = new VCFMenu(this);

    int visibleHeight = 0;

    //private int genotypeBandHeight = EXPANDED_GENOTYPE_HEIGHT;
    private int alleleBandHeight = DEFAULT_ALLELE_HEIGHT;
    //private int savedBandHeight = genotypeBandHeight;
    private List<String> samples;
    private ColorMode coloring = ColorMode.ZYGOSITY;
    private boolean hideReference = true;
    private boolean hideAncestral = false;


    private boolean hideFiltered = true;
    private boolean renderID = true;

    private static float dash[] = {4.0f, 1.0f};
    DecimalFormat numFormat = new DecimalFormat("#.###");

    private static Stroke dashedStroke = new BasicStroke(1.0f, BasicStroke.CAP_BUTT,
            BasicStroke.JOIN_MITER, 1.0f, dash, 0.0f);


    public VCFTrack(ResourceLocator locator, TribbleFeatureSource source) {
        super(locator, source);
        VCFHeader header = (VCFHeader) source.getHeader();
        samples = new ArrayList<String>(header.getGenotypeSamples());

        int visWindow = 2500000;
        if (samples.size() > 1) {
            visWindow = Math.max(100000, 2500000 - 100000 * samples.size());
            this.setDisplayMode(DisplayMode.EXPANDED);
        }
        setRenderID(false);
        setVisibilityWindow(visWindow);
    }

    @Override
    public int getPreferredHeight() {
        return getHeight();
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
            return alleleBandHeight;
        } else {
            return alleleBandHeight + (samples.size() * getGenotypeBandHeight());
        }

    }


    public void setHeight(int height) {
        this.height = Math.max(minimumHeight, height);
    }

    public void setDisplayMode(DisplayMode mode) {
        super.setDisplayMode(mode);
        setHeight(getPreferredHeight());
    }

    @Override
    protected void renderFeatureImpl(RenderContext context, Rectangle rect, PackedFeatures packedFeatures) {
        double locScale = context.getScale();
        Rectangle visibleRectangle = context.getVisibleRect();
        Graphics2D g2D = context.getGraphics();

        visibleHeight = visibleRectangle.height;

        int windowStart = (int) context.getOrigin();


        rect.height = getGenotypeBandHeight();
        rect.y = getTop() + alleleBandHeight;
        colorBackground(g2D, rect, visibleRectangle, false, false);


        List<Feature> features = packedFeatures.getFeatures();
        if (features.size() > 0) {
            byte[] reference = null;

            int lastPX = -1;
            double pXEnd = rect.getMaxX();
            for (Feature feature : features) {
                rect.y = getTop();

                VariantContext variant = (VariantContext) feature;
                char ref = getReference(variant, windowStart, reference);

                // 1 -> 0 based coordinates
                int start = variant.getStart() - 1;
                int end = variant.getEnd();
                int pX = (int) ((start - context.getOrigin()) / locScale);
                int dX = (int) Math.max(2, (end - start) / locScale);

                if (pX > pXEnd) {
                    break;
                }

                if (pX > lastPX) {
                    ZygosityCount zygCounts = getZygosityCounts(variant);

                    rect.height = alleleBandHeight;
                    if (rect.intersects(visibleRectangle)) {
                        if (samples.size() == 0) {
                            renderer.renderVariant(variant, rect, pX, dX, context);
                        } else {
                            AlleleCount alleleCounts = new AlleleCount(zygCounts);
                            renderer.renderAlleleBand(variant, rect, pX, dX, context, hideFiltered,
                                    alleleCounts);
                        }
                    }

                    if (this.getDisplayMode() != Track.DisplayMode.COLLAPSED) {
                        rect.y += rect.height;

                        rect.height = getGenotypeBandHeight();
                        for (String sample : samples) {
                            if (rect.intersects(visibleRectangle)) {
                                if (variant.isSNP()) {
                                    renderer.renderGenotypeBandSNP(variant, context, rect, pX, dX, sample, ref, coloring,
                                            hideFiltered);
                                }
                            }
                            rect.y += rect.height;
                        }
                    }

                }
                lastPX = pX;
            }
        } else {
            rect.height = alleleBandHeight;
            rect.y = getTop();
            g2D.setColor(Color.gray);
            GraphicUtils.drawCenteredText("No Variants Found", rect, g2D);

        }
    }

    private void colorBackground(Graphics2D g2D, Rectangle bandRectangle, Rectangle visibleRectangle, boolean renderNames, boolean isSelected) {
        boolean coloredLast = true;
        for (int i = 0; i < samples.size(); i++) {
            if (isSelected) {
                g2D.setColor(Color.lightGray);
            } else {
                if (coloredLast) {
                    g2D.setColor(new Color(240, 240, 245));
                    coloredLast = false;
                } else {
                    g2D.setColor(Color.white);
                    coloredLast = true;
                }
            }

            if (bandRectangle.intersects(visibleRectangle)) {
                g2D.fillRect(bandRectangle.x, bandRectangle.y, bandRectangle.width, bandRectangle.height);
                if (renderNames) {
                    if (bandRectangle.height >= 10) {
                        String sample = samples.get(i);
                        g2D.setColor(Color.black);
                        GraphicUtils.drawWrappedText(sample, bandRectangle, g2D, false);
                    }
                }
            }
            bandRectangle.y += bandRectangle.height;
        }
    }

    public char getReference(VariantContext variant, int start, byte[] reference) {
        if (hideReference) {
            if (reference == null) {
                return ' ';
            }
            int i = Math.max(1, variant.getStart() - start);
            return (char) reference[i - 1];
        } else if (hideAncestral) {
            String ancestralAllele = variant.getAttributeAsString("AA");
            if (ancestralAllele == null) {
                return '.';
            } else {
                return ancestralAllele.charAt(0);
            }
        } else {
            return ' ';
        }
    }

    public void setHideReference(boolean value) {
        this.hideReference = value;
    }

    public boolean getHideReference() {
        return hideReference;
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

        Rectangle rect = new Rectangle(trackRectangle);
        g2D.clearRect(rect.x, rect.y, rect.width, rect.height);
        g2D.setFont(FontManager.getScalableFont(10));
        if (isSelected()) {
            g2D.setColor(Color.lightGray);
        } else {
            g2D.setColor(Color.white);
        }

        g2D.setColor(Color.black);
        rect.height = alleleBandHeight;
        if (rect.intersects(visibleRectangle)) {
            GraphicUtils.drawWrappedText(getName(), rect, g2D, false); //getName() + " (" + samples.size() + ")", rect, g2D, false);
        }

        rect.y += rect.height;
        rect.height = getGenotypeBandHeight();
        if (getDisplayMode() != Track.DisplayMode.COLLAPSED) {
            colorBackground(g2D, rect, visibleRectangle, true, isSelected());
        }
    }


    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        double currentWindow = frame.getScale() * 1000;

        VariantContext variant = (VariantContext) getFeatureAt(chr, position, y, frame); //getVariantAtPosition(chr, (int) position, frame);
        if (variant != null) {
            String poi = getMousePOI(y);
            if (poi.equals("VARIANTBAND")) {
                return getVariantToolTip(variant);
            } else {
                if (currentWindow < source.getFeatureWindowSize()) {
                    return getSampleToolTip(poi, variant);
                }
            }
        }
        return " ";
    }

    private String getMousePOI(int pY) {
        String poi = "VARIANTBAND";
        if (pY > (getTop() + alleleBandHeight)) {
            int sampleNumber = (pY - getTop() - alleleBandHeight) / getGenotypeBandHeight();

            try {
                poi = samples.get(sampleNumber);
            } catch (IndexOutOfBoundsException ioe) {
                log.error(ioe);

            }
        }
        return poi;
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
        return " ";
    }

    public static enum ColorMode {
        ZYGOSITY, VARIANT, ALLELE
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
        String toolTip = "<b>Summary</b>";
        toolTip = toolTip.concat("<br>Chr:" + variant.getChr());
        toolTip = toolTip.concat("<br>Position:" + variant.getStart());
        toolTip = toolTip.concat("<br>ID: " + id + "<br>");
        toolTip = toolTip.concat("<br><b>Sample Information</b>");
        toolTip = toolTip.concat("<br>Sample: " + sample);
        toolTip = toolTip.concat("<br>Position:" + variant.getStart());

        Genotype genotype = variant.getGenotype(sample);
        if (genotype != null) {
            toolTip = toolTip.concat("<br>Bases: " + genotype.getGenotypeString());
            toolTip = toolTip.concat("<br>Quality: " + numFormat.format(genotype.getPhredScaledQual()));
            toolTip = toolTip.concat("<br>Type: " + genotype.getType());
        }
        if (variant.isFiltered()) {
            toolTip = toolTip.concat("<br>Is Filtered Out: Yes</b>");
            toolTip = toolTip.concat(getFilterTooltip(variant));
        } else {
            toolTip = toolTip.concat("<br>Is Filtered Out: No</b><br>");
        }

        if (genotype != null) {
            toolTip = toolTip.concat(getSampleInfo(genotype) + "<br>");
        }

        toolTip = toolTip.concat("<br>" + "<b>Cursor Location</b>");
        return toolTip;
    }

    private String getVariantToolTip(VariantContext variant) {
        String id = variant.getAttributeAsString(VariantContext.ID_KEY);
        String toolTip = "<b>Summary</b>";
        toolTip = toolTip.concat("<br>Chr:" + variant.getChr());
        toolTip = toolTip.concat("<br>Position:" + variant.getStart());
        toolTip = toolTip.concat("<br>ID: " + id);
        toolTip = toolTip.concat("<br>Reference: " + variant.getReference().toString());
        Set alternates = variant.getAlternateAlleles();
        if (alternates.size() > 0) {
            toolTip = toolTip.concat("<br>Alternate: " + alternates.toString());
        }

        toolTip = toolTip.concat("<br>Qual: " + numFormat.format(variant.getPhredScaledQual()));
        toolTip = toolTip.concat("<br>Type: " + variant.getType());
        if (variant.isFiltered()) {
            toolTip = toolTip.concat("<br>Is Filtered Out: Yes</b>");
            toolTip = toolTip.concat(getFilterTooltip(variant));
        } else {
            toolTip = toolTip.concat("<br>Is Filtered Out: No</b><br>");
        }
        toolTip = toolTip.concat("<br><b>Alleles:</b>");
        toolTip = toolTip.concat(getAlleleToolTip(getZygosityCounts(variant)));
        toolTip = toolTip.concat("<br>Allele Frequency: " + numFormat.format(getAlleleFreq(variant)) + "<br>");
        toolTip = toolTip.concat("<br><b>Genotypes:</b>");
        toolTip = toolTip.concat(getGenotypeToolTip(getZygosityCounts(variant)) + "<br>");
        toolTip = toolTip.concat(getVariantInfo(variant) + "<br>");
        toolTip = toolTip.concat("<br>" + "<b>Cursor Location</b>");
        return toolTip;
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

        for (String sample : samples) {
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
            totalAlleles = samples.size() * 2;
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


    public List<String> getSamples() {
        return samples;
    }

    protected void setSamples(List<String> samples) {
        this.samples = samples;
    }

    public JPopupMenu getPopupMenu(final TrackClickEvent te) {
        Feature f = null;
        if (te.getFrame() != null && te.getFrame().getName() != null) {
            f = getFeatureClosest(te.getChromosomePosition(), te.getMouseEvent().getY(), te.getFrame());
        }
        return menu.getDataPanelMenu(te, f);
    }

    public static void refresh() {
        UIUtilities.invokeOnEventThread(new Runnable() {
            public void run() {
                IGVMainFrame.getInstance().doRefresh();
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
