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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.renderer;

import org.apache.log4j.Logger;
import org.broad.igv.feature.*;
import org.broad.igv.prefs.IGVPreferences;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.LoadedDataInterval;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.SequenceTrack;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.broad.igv.prefs.Constants.*;


/**
 * @author jrobinso
 */
public class SequenceRenderer {


    private static Logger log = Logger.getLogger(SequenceRenderer.class);

    private static final int AMINO_ACID_RESOLUTION = 5;

    public static Map<Character, Color> nucleotideColors;

    public static Map<Character, Color> getNucleotideColors() {
        if (nucleotideColors == null) setNucleotideColors();
        return nucleotideColors;
    }

    private synchronized static void setNucleotideColors() {

        IGVPreferences prefs = PreferencesManager.getPreferences();

        nucleotideColors = new HashMap();

        Color a = ColorUtilities.stringToColor(prefs.get(COLOR_A), new Color(0, 150, 0));
        Color c = ColorUtilities.stringToColor(prefs.get(COLOR_C), Color.blue);
        Color t = ColorUtilities.stringToColor(prefs.get(COLOR_T), Color.red);
        Color g = ColorUtilities.stringToColor(prefs.get(COLOR_G), new Color(209, 113, 5));
        Color n = ColorUtilities.stringToColor(prefs.get(COLOR_N), Color.gray);

        nucleotideColors.put('A', a);
        nucleotideColors.put('a', a);
        nucleotideColors.put('C', c);
        nucleotideColors.put('c', c);
        nucleotideColors.put('T', t);
        nucleotideColors.put('t', t);
        nucleotideColors.put('G', g);
        nucleotideColors.put('g', g);
        nucleotideColors.put('N', n);
        nucleotideColors.put('n', n);
        nucleotideColors.put('-', Color.lightGray);

    }


    private TranslatedSequenceDrawer translatedSequenceDrawer;

    private Strand strand = Strand.POSITIVE;


    public SequenceRenderer() {
        if (nucleotideColors == null) setNucleotideColors();
        translatedSequenceDrawer = new TranslatedSequenceDrawer();
    }


    public void draw(LoadedDataInterval<SequenceTrack.SeqCache> sequenceInterval,
                     RenderContext context,
                     Rectangle trackRectangle,
                     boolean showTranslation,
                     int resolutionThreshold) {


        String chr = context.getChr();
        if (!chr.equals(sequenceInterval.range.chr)) {
            log.error("Chromosome mismatch in sequence track");
            return;
        }

        if (context.getScale() >= resolutionThreshold) {
            // Zoomed out too far to see sequences.  This can happen when in gene list view and one of the frames
            // is zoomed in but others are not
            context.getGraphic2DForColor(UIConstants.LIGHT_GREY).fill(trackRectangle);

        } else {
            double locScale = context.getScale();
            int start = (int) context.getOrigin();
            int end = (int) (start + trackRectangle.width * locScale) + 1;

            SequenceTrack.SeqCache cache = sequenceInterval.getFeatures();
            byte[] seq = cache.seq;
            int sequenceStart = cache.start;
            if (end <= sequenceStart) return;

            //The combined height of sequence and (optionally) colorspace bands
            int untranslatedSequenceHeight = (int) trackRectangle.getHeight();

            if (showTranslation) {
                untranslatedSequenceHeight = (int) (trackRectangle.getHeight() / 4);
                // Draw translated sequence
                Rectangle translatedSequenceRect = new Rectangle(trackRectangle.x, trackRectangle.y + untranslatedSequenceHeight,
                        (int) trackRectangle.getWidth(), (int) trackRectangle.getHeight() - untranslatedSequenceHeight);
                if (context.getScale() < AMINO_ACID_RESOLUTION) {
                    translatedSequenceDrawer.draw(context, sequenceStart, translatedSequenceRect, cache, strand);
                }
            }

            //Rectangle containing the sequence and (optionally) colorspace bands
            Rectangle untranslatedSequenceRect = new Rectangle(trackRectangle.x, trackRectangle.y,
                    (int) trackRectangle.getWidth(), untranslatedSequenceHeight);


            byte[] seqCS = null;

            if (seq != null && seq.length > 0) {
                int yBase = untranslatedSequenceRect.y + 2;
                int yCS = untranslatedSequenceRect.y + 2;
                int dY = untranslatedSequenceRect.height - 4;
                int dX = (int) (1.0 / locScale);
                // Create a graphics to use
                Graphics2D g = context.getGraphics2D("SEQUENCE");

                //dhmay adding check for adequate track height
                int fontSize = Math.min(untranslatedSequenceRect.height, Math.min(dX, 12));
                if (fontSize >= 8) {
                    Font f = FontManager.getFont(Font.BOLD, fontSize);
                    g.setFont(f);
                }

                // Loop through base pair coordinates
                int lastVisibleNucleotideEnd = Math.min(end, seq.length + sequenceStart);
                int lastPx0 = -1;
                int scale = Math.max(1, (int) context.getScale());
                double origin = context.getOrigin();
                for (int loc = start - 1; loc < lastVisibleNucleotideEnd; loc += scale) {
                    int pX0 = (int) ((loc - origin) / locScale);
                    // Skip drawing if we haven't advanced 1 pixel past last nt.  Low zoom
                    if (pX0 > lastPx0) {
                        lastPx0 = pX0;

                        int idx = loc - sequenceStart;
                        if (idx < 0) continue;
                        if (idx >= seq.length) break;

                        char c = (char) seq[idx];
                        if (Strand.NEGATIVE.equals(strand)) c = complementChar(c);

                        Color color = nucleotideColors.get(c);
                        if (fontSize >= 8) {
                            if (color == null) {
                                color = Color.black;
                            }
                            g.setColor(color);
                            drawCenteredText(g, new char[]{c}, pX0, yBase + 2, dX, dY - 2);
                        } else {
                            int bw = Math.max(1, dX - 1);
                            if (color != null) {
                                g.setColor(color);
                                g.fillRect(pX0, yBase, bw, dY);
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Return the complement of a nucleotide or ambiguity code
     *
     * @param inputChar
     * @return
     */
    protected char complementChar(char inputChar) {
        switch (inputChar) {
            case 'A':
                return 'T';
            case 'T':
            case 'U':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
            case 'M':
                return 'K';
            case 'R':
                return 'Y';
            case 'Y':
                return 'R';
            case 'K':
                return 'M';
            case 'V':
                return 'B';
            case 'H':
                return 'D';
            case 'D':
                return 'H';
            case 'B':
                return 'V';

            case 'a':
                return 't';
            case 't':
            case 'u':
                return 'a';
            case 'g':
                return 'c';
            case 'c':
                return 'g';
            case 'm':
                return 'k';
            case 'r':
                return 'y';
            case 'y':
                return 'r';
            case 'k':
                return 'm';
            case 'v':
                return 'b';
            case 'h':
                return 'd';
            case 'd':
                return 'h';
            case 'b':
                return 'v';

            default:
                return inputChar;
        }
    }

    private void drawCenteredText(Graphics2D g, char[] chars, int x, int y, int w, int h) {
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

    public Strand getStrand() {
        return strand;
    }

    public void setStrand(Strand strand) {
        this.strand = strand;
    }

    /**
     * @author Damon May
     * This class draws three amino acid bands representing the 3-frame translation of one strand
     * of the associated SequenceTrack
     */
    public static class TranslatedSequenceDrawer {

        public static final int HEIGHT_PER_BAND = 14;
        public static final int TOTAL_HEIGHT = 3 * HEIGHT_PER_BAND;

        //alternating colors for aminoacids
        public static final Color AA_COLOR_1 = new Color(128, 128, 128);
        public static final Color AA_COLOR_2 = new Color(170, 170, 170);

        public static final Color AA_FONT_COLOR = Color.WHITE;

        //minimum font size for drawing AA characters
        public static final int MIN_FONT_SIZE = 6;
        //minimum vertical buffer around AA characters in band
        public static final int MIN_FONT_VBUFFER = 1;
        //ideal vertical buffer around AA characters in band
        public static final int IDEAL_FONT_VBUFFER = 2;

        protected static final Color STOP_CODON_COLOR = Color.RED;
        protected static final Color METHIONINE_COLOR = Color.GREEN;
        protected static final Color ALT_START_COLOR = new Color(255, 231, 95);

        protected static final Color NUCLEOTIDE_SEPARATOR_COLOR = new Color(150, 150, 150, 120);

        public void draw(RenderContext context, int start, Rectangle trackRectangle, SequenceTrack.SeqCache cache, Strand strand) {

            //each band gets 1/3 of the height, rounded
            int idealHeightPerBand = trackRectangle.height / 3;

            //In this situation, band height is more equal if we tweak things a bit
            if (trackRectangle.height % 3 == 2)
                idealHeightPerBand++;
            int minHeightPerBand = Math.min(idealHeightPerBand, trackRectangle.height - (2 * idealHeightPerBand));

            double locScale = context.getScale();
            double origin = context.getOrigin();

            //Figure out the dimensions of a single box containing an aminoacid, at this zoom
            int oneAcidBoxWidth = getPixelFromChromosomeLocation(context.getChr(), 3, origin, locScale) -
                    getPixelFromChromosomeLocation(context.getChr(), 0, origin, locScale) + 1;
            int oneAcidBoxMinDimension = Math.min(oneAcidBoxWidth, minHeightPerBand);

            //Calculate the font size.  If that's less than MIN_FONT_SIZE, we won't draw amino acids
            int fontSize = 0;
            boolean shouldDrawLetters = false;
            if (oneAcidBoxMinDimension >= 2 * MIN_FONT_VBUFFER + MIN_FONT_SIZE) {
                int idealFontSize = oneAcidBoxMinDimension - 2 * IDEAL_FONT_VBUFFER;

                fontSize = Math.max(idealFontSize, MIN_FONT_SIZE);
                shouldDrawLetters = true;
            }

            boolean shouldDrawNucleotideLines = shouldDrawLetters && oneAcidBoxWidth >= 2.5 * fontSize;

            Rectangle bandRectangle = new Rectangle(trackRectangle.x, 0, trackRectangle.width, 0);
            int heightAlreadyUsed = 0;

            //rf 0
            bandRectangle.y = trackRectangle.y;
            bandRectangle.height = idealHeightPerBand;
            heightAlreadyUsed += bandRectangle.height;

            //This set collects the X positions for nucleotide lines, if we choose to draw them.
            //Technically we could calculate these, but I haven't managed to do that without some wiggle
            Set<Integer> nucleotideLineXPositions = new HashSet<Integer>();

            AminoAcidSequence[] aa = strand == Strand.POSITIVE ? cache.posAA : cache.negAA;

            //only draw nucleotide lines the last time this is called
            drawOneTranslation(context, bandRectangle, 0, shouldDrawLetters, fontSize,
                    nucleotideLineXPositions, aa[0], strand);

            //rf 1
            bandRectangle.y = trackRectangle.y + heightAlreadyUsed;
            bandRectangle.height = idealHeightPerBand;
            heightAlreadyUsed += bandRectangle.height;
            drawOneTranslation(context, bandRectangle, 1, shouldDrawLetters, fontSize,
                    nucleotideLineXPositions, aa[1], strand);

            //rf 2
            bandRectangle.y = trackRectangle.y + heightAlreadyUsed;
            bandRectangle.height = trackRectangle.height - heightAlreadyUsed;
            drawOneTranslation(context, bandRectangle, 2, shouldDrawLetters, fontSize,
                    nucleotideLineXPositions, aa[2], strand);

            if (shouldDrawNucleotideLines) {
                Graphics2D graphicsForNucleotideLines = context.getGraphic2DForColor(NUCLEOTIDE_SEPARATOR_COLOR);
                //use a dashed stroke
                graphicsForNucleotideLines.setStroke(new BasicStroke(1, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL,
                        0, new float[]{1, 2}, 0));

                int topYCoord = trackRectangle.y - 1;
                for (int xVal : nucleotideLineXPositions) {
                    if (xVal >= trackRectangle.x && xVal <= trackRectangle.x + trackRectangle.width)
                        graphicsForNucleotideLines.drawLine(xVal,
                                topYCoord,
                                xVal, topYCoord + trackRectangle.height);
                }
            }
        }

        /**
         * Draw the band representing a translation of the sequence in one reading frame
         *
         * @param context
         * @param bandRectangle
         * @param readingFrame
         * @param shouldDrawLetters
         * @param fontSize
         * @param nucleotideLineXPositions a Set that will accrue all of the x positions that we define here
         */
        protected void drawOneTranslation(RenderContext context,
                                          Rectangle bandRectangle,
                                          int readingFrame,
                                          boolean shouldDrawLetters,
                                          int fontSize,
                                          Set<Integer> nucleotideLineXPositions,
                                          AminoAcidSequence aaSequence,
                                          Strand strand) {

            double locScale = context.getScale();
            double origin = context.getOrigin();

            Graphics2D fontGraphics = context.getGraphics2D("AA_FONT");
            fontGraphics.setColor(AA_FONT_COLOR);

            if (aaSequence != null && aaSequence.hasNonNullSequence()) {
                Graphics2D g = context.getGraphics2D("TRANSLATION");

                //This rectangle holds a single AA glyph. x and width will be updated in the for loop
                Rectangle aaRect = new Rectangle(0, bandRectangle.y, 1, bandRectangle.height);

                //start position for this amino acid. Will increment in for loop below
                int aaSeqStartPosition = aaSequence.getStart(); // + readingFrame;

                //calculated oddness or evenness of first amino acid
                int firstFullAcidIndex = (int) Math.floor((aaSeqStartPosition - readingFrame) / 3);
                boolean odd = (firstFullAcidIndex % 2) == 1;

                if (shouldDrawLetters) {
                    Font f = FontManager.getFont(Font.BOLD, fontSize);
                    g.setFont(f);
                }

                for (CodonAA acid : aaSequence.getSequence()) {
                    if (acid != null) {
                        //calculate x pixel boundaries of this AA rectangle
                        int px = getPixelFromChromosomeLocation(context.getChr(), aaSeqStartPosition, origin, locScale);
                        int px2 = getPixelFromChromosomeLocation(context.getChr(), aaSeqStartPosition + 3,
                                origin, locScale);

                        //if x boundaries of this AA overlap the band rectangle
                        if ((px <= bandRectangle.getMaxX()) && (px2 >= bandRectangle.getX())) {
                            aaRect.x = px;
                            aaRect.width = px2 - px;

                            nucleotideLineXPositions.add(aaRect.x);
                            nucleotideLineXPositions.add(aaRect.x + aaRect.width);

                            char aaSymbol = acid.getAminoAcid().getSymbol();
                            Graphics2D bgGraphics =
                                    context.getGraphic2DForColor(getColorForAminoAcid(aaSymbol, odd, acid.getCodon()));

                            bgGraphics.fill(aaRect);

                            if (shouldDrawLetters) {
                                String acidString = new String(new char[]{aaSymbol});
                                GraphicUtils.drawCenteredText(acidString, aaRect, fontGraphics);
                            }
                        }
                        //need to switch oddness whether we displayed the AA or not,
                        //because oddness is calculated from first AA
                        odd = !odd;

                        aaSeqStartPosition += 3;
                    }
                }

            }
        }

        protected Color getColorForAminoAcid(char acidSymbol, boolean odd, String codon) {

            if (codon.equals("ATG")) {
                return METHIONINE_COLOR;
            } else if (acidSymbol == '*') {
                return STOP_CODON_COLOR;
            } else {
                AminoAcidManager.CodonTable codonTable = AminoAcidManager.getInstance().getCodonTable();
                Set<String> altStartCodons = codonTable.getAltStartCodons();
                if (altStartCodons != null && altStartCodons.contains(codon)) {
                    return ALT_START_COLOR;
                } else {
                    return odd ? AA_COLOR_1 : AA_COLOR_2;
                }
            }
        }

        protected int getPixelFromChromosomeLocation(String chr, int chromosomeLocation, double origin,
                                                     double locationScale) {
            return (int) Math.round((chromosomeLocation - origin) / locationScale);
        }

        public Color getDefaultColor() {
            return Color.BLACK;
        }


    }
}
