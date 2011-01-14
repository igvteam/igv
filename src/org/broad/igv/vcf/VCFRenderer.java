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

package org.broad.igv.vcf;

import org.broad.igv.renderer.FeatureRenderer;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.util.ColorUtilities;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.util.*;

/**
 * User: Jesse Whitworth
 * Date: Jul 16, 2010
 */

public class VCFRenderer { //extends FeatureRenderer {


    private static float alphaValue = (float) .1;
    private static Color colorHomRef = Color.lightGray;
    private static Color colorHomRefAlpha = ColorUtilities.getCompositeColor(colorHomRef.getColorComponents(null), alphaValue);
    private static Color colorHomVar = Color.cyan;
    private static Color colorHomVarAlpha = ColorUtilities.getCompositeColor(colorHomVar.getColorComponents(null), alphaValue);
    private static Color colorHet = Color.blue;
    private static Color colorHetAlpha = ColorUtilities.getCompositeColor(colorHet.getColorComponents(null), alphaValue);
    private static Color colorNoCall = Color.gray;
    private static Color colorNoCallAlpha = ColorUtilities.getCompositeColor(colorNoCall.getColorComponents(null), alphaValue);

    private static int variantWidth = 3;
    static Map<Character, Color> nucleotideColors = new HashMap<Character, Color>();
    private static final Color DARK_GREEN = new Color(30, 120, 30);
    private static final Color BLUE = Color.blue.darker();

    static {
        nucleotideColors.put('A', Color.GREEN);
        nucleotideColors.put('a', Color.GREEN);
        nucleotideColors.put('C', Color.BLUE);
        nucleotideColors.put('c', Color.BLUE);
        nucleotideColors.put('T', Color.RED);
        nucleotideColors.put('t', Color.RED);
        nucleotideColors.put('G', new Color(242, 182, 65));
        nucleotideColors.put('g', new Color(242, 182, 65));
        nucleotideColors.put('N', Color.lightGray);
        nucleotideColors.put('n', Color.lightGray);
        nucleotideColors.put('.', Color.lightGray);
        nucleotideColors.put(null, Color.BLACK);
    }

    public VCFRenderer() {
    }

    public void renderVariantBand(VariantContext variant, Rectangle bandRectangle, int pX0, int dX,
                                  RenderContext context, boolean hideFiltered, VCFTrack.ZygosityCount zygCounts,
                                  int idBandHeight) {


        int pY = (int) bandRectangle.getY() + idBandHeight;
        int dY = (int) bandRectangle.getHeight() - idBandHeight;
        int bottomY = bandRectangle.y + bandRectangle.height;

        float sampleHeight = ((float) dY) / zygCounts.getSampleCount();

        //Allele Percentages
        boolean filtered = hideFiltered && variant.isFiltered();

        int totalHeight = (int) (zygCounts.getTotalCall() * sampleHeight);
        Color c = filtered ? colorHomRefAlpha : colorHomRef;
        context.getGraphic2DForColor(c).fillRect(pX0, bottomY - totalHeight, dX, totalHeight);

        int varHeight = (int) (zygCounts.getVarCall() * sampleHeight);
        c = filtered ? colorHomVarAlpha : colorHomVar;
        context.getGraphic2DForColor(c).fillRect(pX0, bottomY - varHeight, dX, varHeight);

        int hetHeight = (int) (zygCounts.getHet() * sampleHeight);
        c = filtered ? colorHetAlpha : colorHet;
        context.getGraphic2DForColor(c).fillRect(pX0, bottomY - hetHeight, dX, hetHeight);

        //Render ID Name

        String id = variant.getAttributeAsString(VariantContext.ID_KEY);
        if (!id.equals(".") && (idBandHeight > 0)) {
            Graphics2D g = (Graphics2D) context.getGraphics().create();
            Font f = FontManager.getScalableFont(Font.BOLD, 10);
            f = f.deriveFont(AffineTransform.getRotateInstance(-Math.PI / 6.0));
            g.setFont(f);
            g.setColor(Color.black);
            g.drawString(variant.getAttributeAsString(VariantContext.ID_KEY), pX0 + (dX / 2), pY - 1);
            g.dispose();
        }

        context.getGraphic2DForColor(Color.black).drawRect(bandRectangle.x, bandRectangle.y, bandRectangle.width,
                bandRectangle.height);


    }

    public void renderAlleleBand(VariantContext variant, Rectangle bandRectangle, int pX0, int dX,
                                 RenderContext context, boolean hideFiltered, VCFTrack.AlleleCount alleleCounts) {


        int bottomY = bandRectangle.y + bandRectangle.height;

        float allelePercent = alleleCounts.getAllelePercent();
        int alleleBarHeight = (int) (allelePercent * bandRectangle.height);

        boolean filtered = hideFiltered && variant.isFiltered();

        Graphics2D g = context.getGraphic2DForColor(Color.red);
        g.fillRect(pX0, bottomY - alleleBarHeight, dX, alleleBarHeight);

        int remainderHeight = bandRectangle.height - alleleBarHeight;
        g = context.getGraphic2DForColor(Color.lightGray);
        g.fillRect(pX0, bandRectangle.y, dX, remainderHeight);


        context.getGraphic2DForColor(Color.black).drawRect(bandRectangle.x, bandRectangle.y, bandRectangle.width,
                bandRectangle.height);

    }

    public void renderVariant(VariantContext variant, Rectangle bandRectangle, int pX0, int dX, RenderContext context) {


        int bottomY = bandRectangle.y + bandRectangle.height;

        int barHeight = (int) (0.8 * bandRectangle.height);

        Graphics2D g = context.getGraphic2DForColor(BLUE);
        g.fillRect(pX0, bottomY - barHeight, dX, barHeight);

        context.getGraphic2DForColor(Color.black).drawRect(bandRectangle.x, bandRectangle.y, bandRectangle.width,
                bandRectangle.height);

    }

    private AlphaComposite makeComposite(float alpha) {
        int type = AlphaComposite.SRC_OVER;
        return (AlphaComposite.getInstance(type, alpha));
    }

    public void renderGenotypeBandSNP(VariantContext variant, RenderContext context, Rectangle bandRectangle, int pX0, int dX,
                                      String sampleName, char ref, VCFTrack.ColorMode coloring, boolean hideFiltered) {

        int pY = (int) bandRectangle.getY();
        int dY = (int) bandRectangle.getHeight();

        int tOffset = 6;
        int bOffset = 8;
        Graphics2D g = (Graphics2D) context.getGraphics().create();

        if (dX >= 10) {
            if (dY > 24) {
                Font f = FontManager.getScalableFont(Font.BOLD, Math.min(dX, 12));
                g.setFont(f);
            } else if (dY > 18) {
                Font f = FontManager.getScalableFont(Font.BOLD, Math.min(dX, 8));
                tOffset = 4;
                bOffset = 5;
                g.setFont(f);
            }
        }

        boolean isFiltered = variant.isFiltered() && hideFiltered;

        Genotype genotype = variant.getGenotype(sampleName);
        if (genotype != null) {
            char b1 = getFirstBase(genotype.getAlleles().get(0));
            char b2 = getFirstBase(genotype.getAlleles().get(1));
            Color b1Color;
            Color b2Color;

            //Assign proper coloring
            if (coloring == VCFTrack.ColorMode.ZYGOSITY) {
                b1Color = getGenotypeColor(genotype, isFiltered);
                b2Color = b1Color;
            } else if (coloring == VCFTrack.ColorMode.VARIANT) {
                b1Color = getGenotypeColor(genotype, isFiltered);
                b2Color = b1Color;
            } else if (coloring == VCFTrack.ColorMode.ALLELE) {
                b1Color = nucleotideColors.get(b1);
                b2Color = nucleotideColors.get(b2);
            } else {
                b1Color = colorNoCall;
                b2Color = b1Color;
            }

            //Suppress Options
            if ((b1 == ref && b2 == ref) || (variant.isFiltered() && hideFiltered)) {
                g.setComposite(makeComposite(alphaValue));
            }


            //Temp remove to see if no-calls are clearer
            if (b1 == '.' || b2 == '.') {
                dX = dX / 2;
                int offset = dX / 2;
                pX0 = pX0 + offset;
            }
            if (coloring == VCFTrack.ColorMode.ZYGOSITY || coloring == VCFTrack.ColorMode.VARIANT) {
                g.setColor(b1Color);
                g.fillRect(pX0, pY, dX, dY);
            } else {
                g.setColor(b1Color);
                g.fillRect(pX0, pY, (dX / 2), Math.max(1, dY));
                g.setColor(b2Color);
                g.fillRect(pX0 + (dX / 2), pY, (dX / 2), Math.max(1, dY));
            }
            if ((dX >= 10) && (dY >= 18)) {
                if (b1Color == Color.blue) {
                    g.setColor(Color.white);
                } else {
                    g.setColor(Color.black);
                }
                drawCenteredText(g, new char[]{b1}, pX0, pY - tOffset, dX, dY);
                drawCenteredText(g, new char[]{b2}, pX0, pY + (dY / 2) - bOffset, dX, dY);
            }
        }
//        g.dispose();
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
            return colorNoCall;
        } else if (genotype.isHomRef()) {
            return isFiltered ? colorHomRefAlpha : colorHomRef;
        } else if (genotype.isHomVar()) {
            return isFiltered ? colorHomVarAlpha : colorHomVar;
        } else if (genotype.isHet()) {
            return isFiltered ? colorHetAlpha : colorHet;
        }
        return Color.white;
    }

    private void drawCenteredText(Graphics2D g, char[] chars, int x, int y,
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
