/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
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

package org.broad.igv.variant;

import org.apache.log4j.Logger;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.FontManager;
import org.broad.igv.util.ColorUtilities;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * User: Jesse Whitworth
 * Date: Jul 16, 2010
 */

public class VariantRenderer { //extends FeatureRenderer {

    private static Logger log = Logger.getLogger(VariantRenderer.class);

    private static float alphaValue = 0.2f;
    public static Color colorHomRef = new Color(235, 235, 235);
    public static Color colorHomRefAlpha = ColorUtilities.getCompositeColor(colorHomRef.getColorComponents(null), alphaValue);
    public static Color colorHomVar = new Color(0, 245, 255);
    public static Color colorHomVarAlpha = ColorUtilities.getCompositeColor(colorHomVar.getColorComponents(null), alphaValue);
    public static Color colorHet = Color.blue.brighter();  //new Color(107, 30, 115); //Color.blue;
    public static Color colorHetAlpha = ColorUtilities.getCompositeColor(colorHet.getColorComponents(null), alphaValue);
    public static Color colorNoCall = Color.white;
    public static Color colorNoCallAlpha = ColorUtilities.getCompositeColor(colorNoCall.getColorComponents(null), alphaValue);
    public static final Color colorAlleleBand = Color.red;
    public static Color colorAlleleBandAlpha = ColorUtilities.getCompositeColor(colorAlleleBand.getColorComponents(null), alphaValue);
    public static final Color colorAlleleRef = Color.gray;
    public static Color colorAlleleRefAlpha = ColorUtilities.getCompositeColor(colorAlleleRef.getColorComponents(null), alphaValue);
    private static final Color blue = new Color(0, 0, 220);
    public static Color blueAlpha = ColorUtilities.getCompositeColor(blue.getColorComponents(null), alphaValue);

    private static int variantWidth = 3;
    static Map<Character, Color> nucleotideColors = new HashMap<Character, Color>();
    private static final Color DARK_GREEN = new Color(30, 120, 30);

    private VariantTrack track;

    static {
        nucleotideColors.put('A', Color.GREEN);
        nucleotideColors.put('a', Color.GREEN);
        nucleotideColors.put('C', Color.BLUE);
        nucleotideColors.put('c', Color.BLUE);
        nucleotideColors.put('T', Color.RED);
        nucleotideColors.put('t', Color.RED);
        nucleotideColors.put('G', new Color(242, 182, 65));
        nucleotideColors.put('g', new Color(242, 182, 65));
        nucleotideColors.put('N', colorAlleleRef);
        nucleotideColors.put('n', colorAlleleRef);
        nucleotideColors.put('.', colorAlleleRef);
        nucleotideColors.put(null, Color.BLACK);
    }

    public VariantRenderer(VariantTrack track) {
        this.track = track;
    }


    public void renderVariantBand(Variant variant,
                                  Rectangle bandRectangle,
                                  int pX0, int dX,
                                  RenderContext context,
                                  boolean hideFiltered) {

        final int bottomMargin = 0;
        final int topMargin = 3;

        final int bottomY = bandRectangle.y + bandRectangle.height - bottomMargin;

        final int barHeight = bandRectangle.height - topMargin - bottomMargin;

        final boolean filtered = variant.isFiltered();
        final Color alleleColor = filtered ? colorAlleleBandAlpha : colorAlleleBand;

        final double allelePercent = Math.min(1, track.getAllelePercent(variant));
        final int alleleBarHeight;
        final int remainderHeight;

        final Color refColor;
        if (allelePercent <= 0) {
            alleleBarHeight = 0;
            remainderHeight = barHeight;
            refColor = filtered ? colorAlleleRefAlpha : colorAlleleRef;
        } else {
            alleleBarHeight = (int) (allelePercent * barHeight);
            remainderHeight = barHeight - alleleBarHeight;

            refColor = filtered ? blueAlpha : blue;
        }

        Graphics2D g = context.getGraphic2DForColor(alleleColor);
        g.fillRect(pX0, bottomY - alleleBarHeight, dX, alleleBarHeight);

        g = context.getGraphic2DForColor(refColor);
        g.fillRect(pX0, bottomY - alleleBarHeight - remainderHeight, dX, remainderHeight);

    }


    public void renderGenotypeBandSNP(Variant variant, RenderContext context, Rectangle bandRectangle, int pX0, int dX,
                                      String sampleName, VariantTrack.ColorMode coloring, boolean hideFiltered) {

        int pY = (int) bandRectangle.getY();
        int dY = (int) bandRectangle.getHeight();

        int tOffset = 6;
        int bOffset = 8;
        Graphics2D g = (Graphics2D) context.getGraphics().create();

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
            log.error("Now what?");
        } else {
            Color b1Color = Color.gray;
            Color b2Color = Color.gray;
            char b1 = ' ';
            char b2 = ' ';
            //Assign proper coloring
            switch (coloring) {
                case GENOTYPE:

                    b1Color = getGenotypeColor(genotype, isFiltered);
                    b2Color = b1Color;
                    break;

                case ALLELE:
                    final List<Allele> alleleList = genotype.getAlleles();
                    if (alleleList.size() > 0) {
                        b1 = getFirstBase(alleleList.get(0));
                        b1Color = nucleotideColors.get(b1);
                    }
                    if (alleleList.size() > 1) {
                        b2 = getFirstBase(alleleList.get(1));
                        b2Color = nucleotideColors.get(b2);
                    }
                    break;
                case METHYLATION_RATE:

                    final Double goodBaseCount = genotype.getAttributeAsDouble("GB");

                    final Double value = genotype.getAttributeAsDouble("MR");
                    if (goodBaseCount < 10 || value == null) {
                        b1Color = colorNoCall;
                        b2Color = b1Color;

                    } else {
                        float mr = (float) value.doubleValue();
                        //   System.out.printf("position %d methylation-rate: %f%n", variant.getStart(), mr);
                        mr /= 100f;
                        b1Color = convertMethylationRateToColor(mr);
                        b2Color = b1Color;
                    }
                    break;

                default:
                    b1Color = colorNoCall;
                    b2Color = b1Color;
            }


            int y0 = track.getDisplayMode() == Track.DisplayMode.EXPANDED ? pY + 1 : pY;
            int h = Math.max(1, track.getDisplayMode() == Track.DisplayMode.EXPANDED ? dY - 2 : dY);

            if (coloring == VariantTrack.ColorMode.GENOTYPE) {
                g.setColor(b1Color);
                g.fillRect(pX0, y0, dX, h);
            } else {
                // Color by allele
                g.setColor(b1Color);
                g.fillRect(pX0, y0, (dX / 2), h);
                g.setColor(b2Color);
                g.fillRect(pX0 + (dX / 2), y0, (dX / 2), h);
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

    private Color convertMethylationRateToColor(float mr) {
        Color color;

        if (mr >= .25) {
            return Color.getHSBColor((mr - .25f) * (1f / 0.75f), 1, 1);
        } else {
            // use a light grey between 0 and 0.25 brightness to indicate moderate methylation at the site.
            return new Color(1f - mr, 1f - mr, 1f - mr);
        }


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
