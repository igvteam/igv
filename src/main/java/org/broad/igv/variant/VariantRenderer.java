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

package org.broad.igv.variant;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.*;
import org.broad.igv.logging.*;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.color.ColorPalette;
import org.broad.igv.ui.color.ColorTable;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.variant.vcf.VCFVariant;

import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.broad.igv.prefs.Constants.*;

/**
 * @author Jesse Whitworth
 * @api
 * @since Jul 16, 2010
 */
public class VariantRenderer { //extends FeatureRenderer {

    private static final Logger log = LogManager.getLogger(VariantRenderer.class);

    private static final int BOTTOM_MARGIN = 0;
    private static final int TOP_MARGIN = 3;
    private static final float ALPHA_VALUE = 0.2f;

    static Map<Character, Color> nucleotideColors;

    static {
        nucleotideColors = new HashMap<>();
        nucleotideColors.put('A', Color.GREEN);
        nucleotideColors.put('a', Color.GREEN);
        nucleotideColors.put('C', Color.BLUE);
        nucleotideColors.put('c', Color.BLUE);
        nucleotideColors.put('T', Color.RED);
        nucleotideColors.put('t', Color.RED);
        nucleotideColors.put('G', new Color(242, 182, 65));
        nucleotideColors.put('g', new Color(242, 182, 65));
        nucleotideColors.put('N', Color.gray);
        nucleotideColors.put('n', Color.gray);
        nucleotideColors.put('.', Color.gray);
        nucleotideColors.put(null, Color.BLACK);
    }


    private final VariantTrack track;

    private Color colorAlleleBandVar;
    private Color colorAlleleBandVarAlpha;
    private Color colorAlleleBandRef;
    private Color colorAlleleBandRefAlpha;
    private Color colorHomRef;
    private Color colorHomRefAlpha;
    private Color colorHomVar;
    private Color colorHomVarAlpha;
    private Color colorHet;
    private Color colorHetAlpha;
    private Color colorNoCall;
    private Color colorNoCallAlpha;

    private Color nonRefColor = new Color(200, 200, 215);

    //Variant Type Colors
    private Color snpColor = new Color(34, 179, 246);
    private Color mnpColor = new Color(34, 101, 246);
    private Color insertionColor = new Color(225, 32, 243);
    private Color deletionColor = new Color(255, 170, 63);
    private Color mixedColor = new Color(223, 24, 24);
    private Color symbolicColor = new Color(10, 228, 107);
    private Color missingColor = new Color(97, 97, 97);
    private Color lightMissingColor = new Color( 240, 240 ,240);
    //boolean colors
    private Color trueColor = Color.BLACK;
    private Color falseColor =  Color.red;


    private ColorPalette palette = ColorUtilities.getPalette("Pastel 1");
    private ColorTable defaultTagValueTable = new PaletteColorTable(palette);
    private ColorScale infoTagScale = defaultTagValueTable;
    private ColorScale formatTagScale = defaultTagValueTable;

    public VariantRenderer(VariantTrack track) {
        this.track = track;
        updateColors();
    }

    /**
     * Check colors against user prefs and update if neccessary
     */
    private void updateColors() {

        final IGVPreferences prefMgr = PreferencesManager.getPreferences();

        if (!prefMgr.getAsColor(AF_VAR_COLOR).equals(colorAlleleBandVar)) {
            colorAlleleBandVar = prefMgr.getAsColor(AF_VAR_COLOR);
            colorAlleleBandVarAlpha = ColorUtilities.getCompositeColor(colorAlleleBandVar, ALPHA_VALUE);

        }
        if (!prefMgr.getAsColor(AF_REF_COLOR).equals(colorAlleleBandRef)) {
            colorAlleleBandRef = prefMgr.getAsColor(AF_REF_COLOR);
            colorAlleleBandRefAlpha = ColorUtilities.getCompositeColor(colorAlleleBandRef, ALPHA_VALUE);

        }
        if (!prefMgr.getAsColor(HOMREF_COLOR).equals(colorHomRef)) {
            colorHomRef = prefMgr.getAsColor(HOMREF_COLOR);
            colorHomRefAlpha = ColorUtilities.getCompositeColor(colorHomRef, ALPHA_VALUE);
        }
        if (!prefMgr.getAsColor(HOMVAR_COLOR).equals(colorHomVar)) {
            colorHomVar = prefMgr.getAsColor(HOMVAR_COLOR);
            colorHomVarAlpha = ColorUtilities.getCompositeColor(colorHomVar, ALPHA_VALUE);
        }
        if (!prefMgr.getAsColor(HETVAR_COLOR).equals(colorHet)) {
            colorHet = prefMgr.getAsColor(HETVAR_COLOR);
            colorHetAlpha = ColorUtilities.getCompositeColor(colorHet, ALPHA_VALUE);
        }
        if (!prefMgr.getAsColor(NOCALL_COLOR).equals(colorNoCall)) {
            colorNoCall = prefMgr.getAsColor(NOCALL_COLOR);
            colorNoCallAlpha = ColorUtilities.getCompositeColor(colorNoCall, ALPHA_VALUE);
        }

    }

    /**
     * Render the site track (the top, summary view of the site).
     *
     * @param variant
     * @param bandRectangle
     * @param pixelX        Location of the variant in pixels
     * @param xWidth        Width of the variant in pixels
     * @param context
     * @api
     */
    public void renderSiteBand(Variant variant,
                               Rectangle bandRectangle,
                               int pixelX, int xWidth,
                               RenderContext context) {

        updateColors();

        final boolean useAlpha = variant.isFiltered();
        final Color alleleColor;
        final Color refColor;
        double percent;

        Color colorAlleleRef = variant.isNonRef() ? nonRefColor :   track.getColor();
        //colorAlleleRef = new Color(((int)(Math.random()*255)),((int)(Math.random()*255)),((int)(Math.random()*255)));

        Color colorAlleleRefAlpha = useAlpha ? ColorUtilities.getCompositeColor(colorAlleleRef, ALPHA_VALUE) : colorAlleleRef;
        final List<ColorBand> colorBands;
        final VariantTrack.ColorMode siteColorMode = track.getSiteColorMode();
        if (siteColorMode == VariantTrack.ColorMode.METHYLATION_RATE) {
            alleleColor = this.convertMethylationRateToColor((float) variant.getMethylationRate() / 100);
            percent = variant.getCoveredSampleFraction();
            refColor = useAlpha ? colorAlleleRefAlpha : colorAlleleRef;   // Gray
            colorBands = List.of(new ColorBand(refColor, 1.0 - percent), new ColorBand(alleleColor, percent));
        }
        else if (siteColorMode == VariantTrack.ColorMode.ALLELE_FREQUENCY
                || siteColorMode == VariantTrack.ColorMode.ALLELE_FRACTION ) {
            alleleColor = useAlpha ? colorAlleleBandVarAlpha : colorAlleleBandVar; // Red

            double af = siteColorMode == VariantTrack.ColorMode.ALLELE_FREQUENCY ?
                    variant.getAlternateAlleleFrequency() :
                    variant.getAlleleFraction();
            percent = Math.min(1, af);

            if (percent <= 0) {
                percent = 0;
                refColor = useAlpha ? colorAlleleRefAlpha : colorAlleleRef;
            } else {
                refColor = useAlpha ? colorAlleleBandRefAlpha : colorAlleleBandRef;                      // Blue
            }
            colorBands = List.of(new ColorBand(refColor, 1.0 - percent), new ColorBand(alleleColor, percent));
        }
        else if ( track.getSiteColorMode() == VariantTrack.ColorMode.VARIANT_TYPE){
            Color solidColor = switch (variant.getType()) {
                case "NO_VARIATION" -> nonRefColor;
                case "SNP" -> snpColor;
                case "MNP" -> mnpColor;
                case "INDEL" -> {
                    final int refLength = variant.getReference().length();
                    final List<Allele> alternateAlleles = variant.getAlternateAlleles();
                    boolean seenDeletion = false ;
                    boolean seenInsertion = false;
                    for(Allele allele : alternateAlleles) {
                        final int altLength = allele.getBases().length;
                        if(altLength > refLength){
                            seenInsertion = true;
                        } else if(altLength < refLength) {
                            seenDeletion = true;
                        }
                    }
                    if(seenDeletion && seenInsertion){
                        yield mixedColor;
                    } else if( seenDeletion){
                        yield deletionColor;
                    } else {
                        yield insertionColor;
                    }
                }
                case "SYMBOLIC" -> symbolicColor;
                case "MIXED"  -> mixedColor;
                default -> mixedColor;
            };
            refColor = useAlpha ? ColorUtilities.getCompositeColor(solidColor, ALPHA_VALUE) : solidColor;
            colorBands = List.of(new ColorBand(refColor, 1.0));
        } else if(track.getSiteColorMode() == VariantTrack.ColorMode.INFO_FIELD){
            final SelectVcfFieldDialog.ColorResult infoField = track.getColorByInfoField();
            colorBands = getColorsForInfoField(variant, infoField);
        }
        else {
            colorBands = List.of(new ColorBand(track.getColor(), 1.0));
        }


        drawBars(pixelX, calculateBottomYSiteBand(bandRectangle), xWidth, calculateBarHeightSiteBand(bandRectangle), context, colorBands);
    }

    private void drawBars(int pixelX, int bottomY, int xWidth, int height, RenderContext context, List<ColorBand> colorBands) {

        //keep the sum as a double so we don't end up with repeated rounding errors
        double totalFractionConsumed = 0;
        for(ColorBand band: colorBands){
            double bandHeight = (band.fraction * height);
            totalFractionConsumed  += band.fraction;
            totalFractionConsumed = Math.max(0, Math.min(totalFractionConsumed, 1.0));
            double remainderHeight = (1.0 - totalFractionConsumed) * height;
            Graphics2D g = context.getGraphic2DForColor(band.color);

            if(bandHeight > 0) {
                //use double rect to avoid white spaces due to rounding issues
                Rectangle2D bar = new Rectangle2D.Double(pixelX, bottomY - bandHeight -(remainderHeight), xWidth, bandHeight);
                g.fill(bar);
            }
        }
    }

    /**
     * internal record to pass around colors and how much of the bar they should fill
     * @param color color to use
     * @param fraction what fraction of the total bar this particular band should fill.
     */
    private record ColorBand(Color color, double fraction){}

    private List<ColorBand> getColorsForInfoField(Variant variant, SelectVcfFieldDialog.ColorResult infoField){
        if (infoField != null) {
            if(infoField.colors() != null){
                infoTagScale = infoField.colors();
            }
            final String fieldName = infoField.value();
            if (fieldName != null && !fieldName.isEmpty()) {
                final Object header = track.getHeader();
                if (header instanceof VCFHeader vcfHeader && variant instanceof VCFVariant vcfVariant) {
                    VariantContext vc = vcfVariant.getVariantContext();
                    final VCFCompoundHeaderLine headerLine = vcfHeader.getInfoHeaderLine(fieldName);
                    if (headerLine != null) {
                        //handle flag specially
                        if (headerLine.getType() == VCFHeaderLineType.Flag) {
                            return vc.hasAttribute(fieldName)
                                    ? List.of(new ColorBand(infoTagScale.getColor("true"), 1.0))
                                    : List.of(new ColorBand(infoTagScale.getColor("false"), 1.0));
                        } else if (vc.hasAttribute(fieldName)) {
                            List<String> values = vc.getAttributeAsStringList(fieldName, "");
                            if (values.isEmpty()) {
                                return List.of(new ColorBand(missingColor, 1.0));
                            } else {
                                final double percent = 1.0 / values.size();
                                return values.stream()
                                        .map(infoTagScale::getColor)
                                        .map(color -> new ColorBand(color, percent))
                                        .toList();
                            }
                        }
                    } else {
                        // we don't know anything about this attribute so the best we can do is pick a consistent
                        final String value = variant.getAttributeAsString(fieldName);
                        Color color = (value == null || value.isEmpty() || value.equals(VCFConstants.MISSING_VALUE_v4)) // best guess at missing value here
                                ? missingColor
                                : infoTagScale.getColor(value);
                        return List.of(new ColorBand(color, 1.0));
                    }
                }
            }
        }
        return List.of(new ColorBand(missingColor, 1.0));
    }

    private List<ColorBand> getColorsForFormatField(Genotype genotype, SelectVcfFieldDialog.ColorResult formatField){
        if (formatField != null) {
            if(formatField.colors() != null){
                formatTagScale = formatField.colors();
            }
            final String fieldName = formatField.value();
            if (fieldName != null && !fieldName.isEmpty()) {
                final Object header = track.getHeader();
                if (header instanceof VCFHeader vcfHeader) {
                    final VCFCompoundHeaderLine headerLine = vcfHeader.getFormatHeaderLine(fieldName);
                    if (headerLine != null) {
                      if (genotype.getAttributes().get(fieldName) != null || fieldName.equals("GT")) {
                            List<String> values = fieldName.equals("GT")  // special case for GT TODO use numericGT
                                    ? List.of(genotype.getGenotypeString())
                                    :genotype.getAttributeAsStringList(fieldName);
                            if (values.isEmpty()) {
                                return List.of(new ColorBand(lightMissingColor, 1.0));
                            } else {
                                final double percent = 1.0 / values.size();
                                return values.stream()
                                        .map(formatTagScale::getColor)
                                        .map(color -> new ColorBand(color, percent))
                                        .toList();
                            }
                        }
                    } else {
                        // we don't know anything about this attribute so the best we can do is pick a consistent
                        final Object value = genotype.getAttributes().get(fieldName);
                        Color color = (value == null || value.toString().isEmpty() || value.equals(VCFConstants.MISSING_VALUE_v4)) // best guess at missing value here
                                ? lightMissingColor
                                : formatTagScale.getColor(value.toString());
                        return List.of(new ColorBand(color, 1.0));
                    }
                }
            }
        }
        return List.of(new ColorBand(missingColor, 1.0));
    }

    protected int calculateBottomYSiteBand(Rectangle bandRectangle) {
        return bandRectangle.y + bandRectangle.height - BOTTOM_MARGIN;
    }

    protected int calculateBarHeightSiteBand(Rectangle bandRectangle) {
        return bandRectangle.height - TOP_MARGIN - BOTTOM_MARGIN;
    }


    public void renderGenotypeBandSNP(Variant variant, RenderContext context, Rectangle bandRectangle, int pX0, int dX,
                                      String sampleName, VariantTrack.ColorMode coloring, boolean hideFiltered) {

        updateColors();

        int pY = (int) bandRectangle.getY();
        int dY = (int) bandRectangle.getHeight();

        int tOffset = 6;
        int bOffset = 8;
        Graphics2D g = context.getGraphics2D("GENOTYPE");

        if (dX >= 10) {
            if (dY > 24) {
                Font f = FontManager.getFont(Font.BOLD, Math.min(dX, 12));
                g.setFont(f);
            } else if (dY > 18) {
                Font f = FontManager.getFont(Font.BOLD, Math.min(dX, 8));
                tOffset = 4;
                bOffset = 5;
                g.setFont(f);
            }
        }

        boolean isFiltered = variant.isFiltered() && hideFiltered;

        Genotype genotype = variant.getGenotype(sampleName);

        if (genotype == null) {
            log.error("Genotype not found for sample " + sampleName);
        } else {
//            Color b1Color = Color.gray;
//            Color b2Color = Color.gray;

            //Assign proper coloring
            List<ColorBand> colorBands = switch (coloring) {
                case GENOTYPE -> List.of(new ColorBand(getGenotypeColor(genotype, isFiltered), 1.0));
                case METHYLATION_RATE -> {
                    final double goodBaseCount = genotype.getAttributeAsDouble("GB");
                    Color b1Color = colorNoCall;
                    Color b2Color = b1Color;
                    final double value = genotype.getAttributeAsDouble("MR");
                    if (!Double.isNaN(goodBaseCount) && !Double.isNaN(value)) {
                        if (goodBaseCount < VariantTrack.METHYLATION_MIN_BASE_COUNT || Double.isNaN(value)) {
                            b1Color = colorNoCall;
                            b2Color = b1Color;
                        } else {
                            float mr = (float) value;
                            mr /= 100f;
                            b1Color = convertMethylationRateToColor(mr);
                            b2Color = b1Color;
                        }
                    } else {
                        log.error("GB and MR fields must be defined for all records in a VCF methylation file.");
                    }
                    yield List.of(new ColorBand(b1Color, 0.5), new ColorBand(b2Color, 0.5));
                }
                case FORMAT_FIELD -> getColorsForFormatField(genotype, track.getColorByFormatField());
                default ->  List.of(new ColorBand(missingColor, 1.0));
            };


            int y0 = track.getDisplayMode() == Track.DisplayMode.EXPANDED ? pY + 1 : pY;
            int h = Math.max(1, track.getDisplayMode() == Track.DisplayMode.EXPANDED ? dY - 2 : dY);


            if (coloring == VariantTrack.ColorMode.METHYLATION_RATE || colorBands.size() == 2) {
                // This is just drawBars but horizontal instead of vertical
                // NOTE: I think it used to be for alleles but it was co-opted for methylation and
                // the original allele use was lost
                // Color by allele
                colorBands.getFirst();
                g.setColor(colorBands.getFirst().color());
                g.fillRect(pX0, y0, (dX / 2), h);
                g.setColor(colorBands.getLast().color());
                g.fillRect(pX0 + (dX / 2), y0, (dX / 2), h);
            } else {
                g.setColor(colorBands.getFirst().color());
                g.fillRect(pX0, y0, dX, h);
                 //todo drawBars(pX0, y0 - h, dX, dY, context, colorBands);
            }

            doMysteryThing(pX0, dX, dY, colorBands.getFirst().color(), g, pY, tOffset, bOffset);
        }
    }

    //TODO this seems like it's doing nothing but it's hard to tell
    private static void doMysteryThing(int pX0, int dX, int dY, Color b1Color, Graphics2D g, int pY, int tOffset, int bOffset) {
        if ((dX >= 10) && (dY >= 18)) {
            g.setColor(b1Color == Color.blue ? Color.white : Color.black);
            drawCenteredText(g, new char[]{' '}, pX0, pY - tOffset, dX, dY);
            drawCenteredText(g, new char[]{' '}, pX0, pY + (dY / 2) - bOffset, dX, dY);
        }
    }

    private Color convertMethylationRateToColor(float mr) {
        Color color;
        /*
        if (mr >= .25) {
            return Color.getHSBColor((mr - .25f) * (1f / 0.75f), 1, 1);
        } else {
            // use a light grey between 0 and 0.25 brightness to indicate moderate methylation at the site.
            return new Color(1f - mr, 1f - mr, 1f - mr);
        }
          */
        float v = 1.5f + (mr / 2); //add one to have a floor. getHSBColor removes the floor to yield a fraction betweeen 0 and 1.
        return Color.getHSBColor(v, .75f, 1);
    }

    public char getFirstBase(Allele allele) {
        byte[] bases = allele.getBases();
        if (bases.length > 0) {
            return (char) bases[0];
        } else {
            return '.';
        }
    }

    public Color getGenotypeColor(Genotype genotype, boolean isFiltered) {
        if (genotype.isNoCall()) {
            return isFiltered ? colorNoCallAlpha : colorNoCall;
        } else if (genotype.isHomRef()) {
            return isFiltered ? colorHomRefAlpha : colorHomRef;
        } else if (genotype.isHomVar()) {
            return isFiltered ? colorHomVarAlpha : colorHomVar;
        } else if (genotype.isHet()) {
            return isFiltered ? colorHetAlpha : colorHet;
        }
        return Color.white;
    }

    private static void drawCenteredText(Graphics2D g, char[] chars, int x, int y,
                                  int w, int h) {

        // Get measures needed to center the message
        FontMetrics fm = g.getFontMetrics();

        // How many pixels wide is the string
        int msg_width = fm.charsWidth(chars, 0, 1);

        // How far above the baseline can the font go?
        int ascent = fm.getMaxAscent();

        // How far below the baseline?
        int descent = fm.getMaxDescent();

        // Use the string width to find the starting point
        int msgX = x + w / 2 - msg_width / 2;

        // Use the vertical height of this font to find
        // the vertical starting coordinate
        int msgY = y + h / 2 - descent / 2 + ascent / 2;

        g.drawChars(chars, 0, 1, msgX, msgY);

    }

}
