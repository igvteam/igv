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
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.AminoAcid;
import org.broad.igv.feature.AminoAcidManager;
import org.broad.igv.feature.AminoAcidSequence;
import org.broad.igv.feature.Strand;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.SOLIDUtils;

import java.awt.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;


/**
 * @author jrobinso
 */
public class SequenceRenderer {


    private static Logger log = Logger.getLogger(SequenceRenderer.class);

    //Maximum scale at which the track is displayed
    //public static final int MAX_SCALE_FOR_RENDER = 1000;
    public static final int AMINO_ACID_RESOLUTION = 5;

    private static Map<Character, Color> nucleotideColors;

    protected TranslatedSequenceDrawer translatedSequenceDrawer;

    //are we rendering positive or negative strand?
    protected Strand strand = Strand.POSITIVE;

    //Have we successfully downloaded sequence info?
    private boolean hasSequence = true;

    public SequenceRenderer() {

        if(nucleotideColors == null) setNucleotideColors();

        translatedSequenceDrawer = new TranslatedSequenceDrawer();
    }


    public static Map<Character, Color> getNucleotideColors() {

        if(nucleotideColors == null) setNucleotideColors();
        return nucleotideColors;

    }

    public static void setNucleotideColors() {

        PreferenceManager prefs = PreferenceManager.getInstance();

        nucleotideColors = new HashMap();

        Color a = ColorUtilities.stringToColor(prefs.get(PreferenceManager.COLOR_A), new Color(0, 150, 0));
        Color c = ColorUtilities.stringToColor(prefs.get(PreferenceManager.COLOR_C), Color.blue);
        Color t = ColorUtilities.stringToColor(prefs.get(PreferenceManager.COLOR_T),  Color.red);
        Color g = ColorUtilities.stringToColor(prefs.get(PreferenceManager.COLOR_G),  Color.gray);
        Color n = ColorUtilities.stringToColor(prefs.get(PreferenceManager.COLOR_N),  Color.gray);

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

    }

    /**
     * @param context
     * @param trackRectangle
     * @param showColorSpace
     * @param showTranslation Should we show the translated amino acids?
     */
    public void draw(RenderContext context, Rectangle trackRectangle,
                     boolean showColorSpace, boolean showTranslation,
                     int resolutionThreshold) {


        if (context.getScale() >= resolutionThreshold) {
            // Zoomed out too far to see sequences.  This can happen when in gene list view and one of the frames
            // is zoomed in but others are not
            context.getGraphic2DForColor(UIConstants.LIGHT_GREY).fill(trackRectangle);

        } else {
            double locScale = context.getScale();
            double origin = context.getOrigin();
            String chr = context.getChr();
            //String genomeId = context.getGenomeId();
            Genome genome = GenomeManager.getInstance().getCurrentGenome();

            //The location of the first base that is loaded, which may include padding around what's visible
            int start = Math.max(0, (int) origin - 1);
            //The location of the last base that is loaded
            int end = (int) (origin + trackRectangle.width * locScale) + 1;

            if (end <= start) return;

            int firstCodonOffset = 0;
            int lastCodonOffset = 0;


            //If we're translating, we need to start with the first bp of the first codon, in frame 3, and
            //end with the last bp of the last codon, in frame 1
            if (showTranslation) {
                if (start > 1) {
                    firstCodonOffset = 2;
                    start -= firstCodonOffset;
                }

                lastCodonOffset = 2;
                end += lastCodonOffset;
            }

            byte[] seq = genome.getSequence(chr, start, end);
            if (seq == null) {
                this.hasSequence = false;
                return;
            } else {
                this.hasSequence = true;
            }

            //The combined height of sequence and (optionally) colorspace bands
            int untranslatedSequenceHeight = (int) trackRectangle.getHeight();


            if (showTranslation) {
                untranslatedSequenceHeight = showColorSpace ? (int) trackRectangle.getHeight() / AMINO_ACID_RESOLUTION * 2 :
                        (int) (trackRectangle.getHeight() / 4);
                // Draw translated sequence
                Rectangle translatedSequenceRect = new Rectangle(trackRectangle.x, trackRectangle.y + untranslatedSequenceHeight,
                        (int) trackRectangle.getWidth(), (int) trackRectangle.getHeight() - untranslatedSequenceHeight);
                if (context.getScale() < AMINO_ACID_RESOLUTION) {
                    translatedSequenceDrawer.draw(context, start, translatedSequenceRect, seq, strand);
                }
            }

            //Rectangle containing the sequence and (optionally) colorspace bands
            Rectangle untranslatedSequenceRect = new Rectangle(trackRectangle.x, trackRectangle.y,
                    (int) trackRectangle.getWidth(), untranslatedSequenceHeight);


            byte[] seqCS = null;
            if (showColorSpace) {
                seqCS = SOLIDUtils.convertToColorSpace(seq);
            }

            if (seq != null && seq.length > 0) {
                int hCS = (showColorSpace ? untranslatedSequenceRect.height / 2 : 0);
                int yBase = hCS + untranslatedSequenceRect.y + 2;
                int yCS = untranslatedSequenceRect.y + 2;
                int dY = (showColorSpace ? hCS : untranslatedSequenceRect.height) - 4;
                int dX = (int) (1.0 / locScale);
                // Create a graphics to use
                Graphics2D g = (Graphics2D) context.getGraphics().create();
                if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
                    g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
                }

                //dhmay adding check for adequate track height
                int fontSize = Math.min(untranslatedSequenceRect.height, Math.min(dX, 12));
                if (fontSize >= 8) {
                    Font f = FontManager.getFont(Font.BOLD, fontSize);
                    g.setFont(f);
                }

                // Loop through base pair coordinates
                int firstVisibleNucleotideStart = start;
                int lastVisibleNucleotideEnd = Math.min(end, seq.length + start);
                int lastPx0 = -1;
                int scale = Math.max(1, (int) context.getScale());
                for (int loc = firstVisibleNucleotideStart; loc < lastVisibleNucleotideEnd; loc += scale) {
                    for (; loc < lastVisibleNucleotideEnd; loc++) {
                        int idx = loc - start;
                        int pX0 = (int) ((loc - origin) / locScale);
                        if (pX0 > lastPx0) {
                            lastPx0 = pX0;
                            char c = (char) seq[idx];
                            if (Strand.NEGATIVE.equals(strand))
                                c = complementChar(c);
                            Color color = nucleotideColors.get(c);
                            if (fontSize >= 8) {
                                if (color == null) {
                                    color = Color.black;
                                }
                                g.setColor(color);
                                drawCenteredText(g, new char[]{c}, pX0, yBase + 2, dX, dY - 2);
                                if (showColorSpace) {
                                    // draw color space #.  Color space is shifted to be between bases as it represents
                                    // two bases.
                                    g.setColor(Color.black);
                                    String cCS = String.valueOf(seqCS[idx]);
                                    drawCenteredText(g, cCS.toCharArray(), pX0 - dX / 2, yCS + 2, dX, dY - 2);
                                }
                            } else {
                                int bw = Math.max(1, dX - 1);
                                if (color != null) {
                                    g.setColor(color);
                                    g.fillRect(pX0, yBase, bw, dY);
                                }
                            }
                        }
                        break;
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

    public boolean hasSequence() {
        return this.hasSequence;
    }


    /**
     * @author Damon May
     *         This class draws three amino acid bands representing the 3-frame translation of one strand
     *         of the associated SequenceTrack
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

        protected static final Color NUCLEOTIDE_SEPARATOR_COLOR = new Color(150, 150, 150, 120);

        /**
         * @param context
         * @param start          Must be the first base involved in any codon that's even partially visible
         * @param trackRectangle
         * @param seq
         */
        public void draw(RenderContext context, int start, Rectangle trackRectangle, byte[] seq, Strand strand) {
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

            //only draw nucleotide lines the last time this is called
            drawOneTranslation(context, start, bandRectangle, 0, shouldDrawLetters, fontSize,
                    nucleotideLineXPositions, seq, strand);

            //rf 1
            bandRectangle.y = trackRectangle.y + heightAlreadyUsed;
            bandRectangle.height = idealHeightPerBand;
            heightAlreadyUsed += bandRectangle.height;
            drawOneTranslation(context, start, bandRectangle, 1, shouldDrawLetters, fontSize,
                    nucleotideLineXPositions, seq, strand);

            //rf 2
            bandRectangle.y = trackRectangle.y + heightAlreadyUsed;
            bandRectangle.height = trackRectangle.height - heightAlreadyUsed;
            drawOneTranslation(context, start, bandRectangle, 2, shouldDrawLetters, fontSize,
                    nucleotideLineXPositions, seq, strand);

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
         * @param start                    the index of the first base in seq.  Should be the first nucleotide that's in a codon
         *                                 that's even partially visible, in any frame
         * @param bandRectangle
         * @param readingFrame
         * @param shouldDrawLetters
         * @param fontSize
         * @param nucleotideLineXPositions a Set that will accrue all of the x positions that we define here
         * @param seq                      nucleotide sequence starting at start
         *                                 for the beginning and end of aminoacid boxes
         */
        protected void drawOneTranslation(RenderContext context, int start,
                                          Rectangle bandRectangle, int readingFrame,
                                          boolean shouldDrawLetters, int fontSize,
                                          Set<Integer> nucleotideLineXPositions, byte[] seq,
                                          Strand strand) {

            double locScale = context.getScale();
            double origin = context.getOrigin();

            Graphics2D fontGraphics = (Graphics2D) context.getGraphic2DForColor(AA_FONT_COLOR).create();
            if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
                fontGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
            }

            //The start location of the first codon that overlaps this region
            int readingFrameOfFullSeq = start % 3;
            int indexOfFirstCodonStart = readingFrame - readingFrameOfFullSeq;
            if (indexOfFirstCodonStart < 0)
                indexOfFirstCodonStart += 3;

            if (seq != null && seq.length > 0) {
                Graphics2D g = (Graphics2D) context.getGraphics().create();
                if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.ENABLE_ANTIALISING)) {
                    g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
                }

                String nucSequence = new String(seq, indexOfFirstCodonStart, seq.length - indexOfFirstCodonStart);
                AminoAcidSequence aaSequence = AminoAcidManager.getInstance().
                        getAminoAcidSequence(strand, start + indexOfFirstCodonStart, nucSequence);

                if ((aaSequence != null) && aaSequence.hasNonNullSequence()) {
                    //This rectangle holds a single AA glyph. x and width will be updated in the for loop
                    Rectangle aaRect = new Rectangle(0, bandRectangle.y, 1, bandRectangle.height);

                    //start position for this amino acid. Will increment in for loop below
                    int aaSeqStartPosition = aaSequence.getStartPosition();

                    //calculated oddness or evenness of first amino acid
                    int firstFullAcidIndex = (int) Math.floor((aaSeqStartPosition - readingFrame) / 3);
                    boolean odd = (firstFullAcidIndex % 2) == 1;

                    if (shouldDrawLetters) {
                        Font f = FontManager.getFont(Font.BOLD, fontSize);
                        g.setFont(f);
                    }

                    for (AminoAcid acid : aaSequence.getSequence()) {
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

                                Graphics2D bgGraphics =
                                        context.getGraphic2DForColor(getColorForAminoAcid(acid.getSymbol(), odd));

                                bgGraphics.fill(aaRect);

                                if (shouldDrawLetters) {
                                    String acidString = new String(new char[]{acid.getSymbol()});
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
        }

        protected Color getColorForAminoAcid(char acidSymbol, boolean odd) {
            switch (acidSymbol) {
                case 'M':
                    return METHIONINE_COLOR;
                case '*':
                    return STOP_CODON_COLOR;
                default:
                    return odd ? AA_COLOR_1 : AA_COLOR_2;
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
