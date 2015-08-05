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
package org.broad.igv.sam;

import org.broad.igv.feature.LocusScore;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.renderer.SequenceRenderer;
import org.broad.igv.tdf.TDFDataSource;
import org.broad.igv.tdf.TDFReader;
import org.broad.igv.track.*;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Developmental track to explore display of conservation scores.
 *
 * @author jrobinso
 */
public class EWigTrack extends AbstractTrack {

    //    double dataMax;
    char[] nucleotides = {'A', 'C', 'G', 'T'};
    public static Color grey1 = new Color(230, 230, 230);
    Map<Character, TDFDataSource> baseSources;
    TDFDataSource scoreSource;

    public EWigTrack(ResourceLocator locator, Genome genome) {
        super(locator);

        TDFReader reader = TDFReader.getReader(locator.getPath());
        scoreSource = new TDFDataSource(reader, 4, "Pi", genome);

        setDataRange(new DataRange(0, 0, 10));
        baseSources = new HashMap();
        for (int i = 0; i < 4; i++) {
            TDFDataSource src = new TDFDataSource(reader, i, Character.toString(nucleotides[i]), genome);
            baseSources.put(nucleotides[i], src);
        }
    }

    public void render(RenderContext context, Rectangle rect) {
        paint(context, rect);
    }

    public boolean isLogNormalized() {
        return false;
    }

    public float getRegionScore(String chr, int start, int end, int zoom, RegionScoreType type, String frameName) {
        return 0;
    }

    private void paint(RenderContext context, Rectangle rect) {

        // The total score
        List<LocusScore> scores = scoreSource.getSummaryScoresForRange(context.getChr(),
                (int) context.getOrigin(),
                (int) context.getEndLocation(),
                context.getZoom());
        Map<Character, List<LocusScore>> nScores = new HashMap();

        for (Character c : nucleotides) {
            nScores.put(c, baseSources.get(c).getSummaryScoresForRange(context.getChr(),
                    (int) context.getOrigin(),
                    (int) context.getEndLocation(),
                    context.getZoom()));
        }

        Map<Character, Color> nucleotideColors = SequenceRenderer.getNucleotideColors();

        for (int idx = 0; idx < scores.size(); idx++) {

            LocusScore score = scores.get(idx);
            int startPosition = score.getStart();
            int endPosition = score.getEnd();

            int pX = (int) (rect.getX() + (startPosition - context.getOrigin()) / context.getScale());
            int dX = Math.max(1,
                    (int) (rect.getX() + (endPosition - context.getOrigin()) / context.getScale()) - pX);
            if (dX > 4) {
                dX -= 2;
                pX++;
            }

            if (pX + dX < 0) {
                continue;
            } else if (pX > context.getVisibleRect().getMaxX()) {
                break;
            }


            float totalCount = score.getScore();
            int pY = (int) rect.getMaxY() - 1;

            float dataMax = getDataRange().getMaximum();

            double height = (totalCount * rect.getHeight() / dataMax);
            height = Math.min(height, rect.height - 1);

            for (char c : nucleotides) {
                try {
                    LocusScore ns = nScores.get(c).get(idx);
                    float count = ns.getScore() * totalCount;

                    //pY = drawBar(context, idx, rect, dataMax, pY, pX, dX, c, count);
                    Graphics2D tGraphics = context.getGraphic2DForColor(nucleotideColors.get(c));

                    int h = (int) Math.round(count * height / totalCount);
                    h = Math.min(pY - rect.y, h);
                    int baseY = (int) (pY - h);

                    if (h > 0) {
                        tGraphics.fillRect(pX, baseY, dX, h);
                    }
                    pY = baseY;
                } catch (Exception e) {
                    e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }

        }

        // Draw border
        context.getGraphic2DForColor(Color.gray).drawLine(
                rect.x, rect.y + rect.height,
                rect.x + rect.width, rect.y + rect.height);

        // Draw scale
        /*
        DataRange range = getDataRange();
        if (range != null) {
            Graphics2D g = context.getGraphic2DForColor(Color.black);
            Font font = g.getFont();
            Font smallFont = FontManager.getScalableFont(8);
            try {
                g.setFont(smallFont);
                String scale = "Scale: " + (int) range.getMinimum() + " - " +
                        (int) range.getMaximum();
                g.drawString(scale, rect.x + 10, rect.y + 10);

            } finally {
                g.setFont(font);
            }
        }*/

    }

    public JPopupMenu getPopupMenu(
            final MouseEvent evt) {

        JPopupMenu popupMenu = new IGVPopupMenu();

        JLabel popupTitle = new JLabel("  " + getName(), JLabel.CENTER);

        Font newFont = popupMenu.getFont().deriveFont(Font.BOLD, 12);
        popupTitle.setFont(newFont);
        if (popupTitle != null) {
            popupMenu.add(popupTitle);
        }

        popupMenu.addSeparator();

        // addSortMenuItem(popupMenu);
        // addPackMenuItem(popupMenu);
        // addShadeBaseMenuItem(popupMenu);
        // addCopyToClipboardItem(popupMenu, evt);
        // addGoToMate(popupMenu, evt);
        // popupMenu.addSeparator();


        //JLabel trackSettingsHeading = new JLabel("  Track Settings",
        //        JLabel.LEFT);
        //trackSettingsHeading.setFont(newFont);

        //popupMenu.add(trackSettingsHeading);

        ArrayList<Track> tmp = new ArrayList();
        tmp.add(this);
        popupMenu.add(TrackMenuUtils.getTrackRenameItem(tmp));
        popupMenu.add(TrackMenuUtils.getChangeTrackHeightItem(tmp));
        popupMenu.add(TrackMenuUtils.getDataRangeItem(tmp));


        return popupMenu;
    }

//    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
//        return null;
//    }
}
