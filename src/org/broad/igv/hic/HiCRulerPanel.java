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

package org.broad.igv.hic;

import org.apache.log4j.Logger;

import org.broad.igv.feature.Chromosome;
import org.broad.igv.hic.track.HiCGridAxis;
import org.broad.igv.ui.FontManager;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.io.Serializable;
import java.text.DecimalFormat;

/**
 * @author jrobinso
 */
public class HiCRulerPanel extends JPanel implements Serializable {

    private static Logger log = Logger.getLogger(HiCRulerPanel.class);

    enum Orientation {HORIZONTAL, VERTICAL}


    private HiC hic;
    private Orientation orientation;

    private Font tickFont = FontManager.getFont(Font.BOLD, 9);
    private Font spanFont = FontManager.getFont(Font.BOLD, 12);


    Context context;

    /**
     * Empty constructor for form builder
     */
    public HiCRulerPanel() {
    }

    public HiCRulerPanel(HiC hic) {
        this.hic = hic;
    }

    public void setFrame(Context frame, Orientation orientation) {
        this.context = frame;
        this.orientation = orientation;
    }


    @Override
    protected void paintComponent(Graphics g) {

        super.paintComponent(g);

        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        Graphics2D g2D = (Graphics2D) g;


        if (context == null || hic.zd == null) return;

        g.setColor(Color.black);

        AffineTransform t = g2D.getTransform();

        if (orientation == Orientation.VERTICAL) {
            AffineTransform rotateTransform = new AffineTransform();
            rotateTransform.quadrantRotate(-1);
            g2D.transform(rotateTransform);
        }

        // Clear panel
        drawTicks(g2D);
        drawChr(g2D);

        g2D.setTransform(t);


    }

    private void drawChr(Graphics g) {
        int w = isHorizontal() ? getWidth() : getHeight();
        int h = isHorizontal() ? getHeight() : getWidth();

        g.setFont(spanFont);

        Chromosome chromosome = context.getChromosome();

        if (chromosome != null) {
            if (chromosome.getName().equals("All")) {

            } else {
                String rangeString = chromosome.getName();
                int strWidth = g.getFontMetrics().stringWidth(rangeString);
                int strPosition = (w - strWidth) / 2;

                if (!isHorizontal()) strPosition = -strPosition;

                int vPos = h - 35;
                g.drawString(rangeString, strPosition, vPos);
            }
        }

    }

    private boolean isHorizontal() {
        return orientation == Orientation.HORIZONTAL;
    }


    private void drawTicks(Graphics g) {

        int w = isHorizontal() ? getWidth() : getHeight();
        int h = isHorizontal() ? getHeight() : getWidth();

        if (w < 50 || context.getScaleFactor() == 0) {
            return;
        }

        g.setFont(tickFont);

        Chromosome chromosome = context.getChromosome();

        if (chromosome == null) return;

        if (chromosome.getName().equals("All")) {
            int x1 = 0;
            Chromosome[] chromosomes = hic.getChromosomes();
            // Index 0 is whole genome
            int genomeCoord = 0;
            for (int i = 1; i < chromosomes.length; i++) {
                Chromosome c = chromosomes[i];
                genomeCoord += (c.getLength() / 1000);

                int xBin = hic.zd.getxGridAxis().getBinNumberForGenomicPosition(genomeCoord);
                int x2 = (int) (xBin * hic.xContext.getScaleFactor());

                int x = (x1 + x2) / 2;
                int strWidth = g.getFontMetrics().stringWidth(c.getName());
                int strPosition = isHorizontal() ? x - strWidth / 2 : -x - strWidth / 2;
                g.drawString(c.getName(), strPosition, h - 15);

                int xpos = (orientation == Orientation.HORIZONTAL ? x2 : -x2);
                g.drawLine(xpos, h - 10, xpos, h - 2);

                x1 = x2;
            }
        } else {

            HiCGridAxis axis = isHorizontal() ? hic.zd.getxGridAxis() : hic.zd.getyGridAxis();

            int binRange = (int) (w / context.getScaleFactor());
            int binOrigin = context.getBinOrigin();     // <= by definition at left/top of panel

            int genomeOrigin = axis.getGenomicStart(binOrigin);

            int genomeEnd = axis.getGenomicEnd(binOrigin + binRange);

            int range = genomeEnd - genomeOrigin;


            TickSpacing ts = findSpacing(range, false);
            double spacing = ts.getMajorTick();

            // Find starting point closest to the current origin
            int maxX = context.getChromosome().getLength();
            int nTick = (int) (genomeOrigin / spacing) - 1;
            int genomePosition = (int) (nTick * spacing);

            //int x = frame.getScreenPosition(genomeTickNumber);
            int binNUmber = axis.getBinNumberForGenomicPosition(genomePosition);

            int x = (int) ((binNUmber - binOrigin) * context.getScaleFactor());

            while (genomePosition < maxX && x < w) {
                genomePosition = (int) (nTick * spacing);

                // x = frame.getScreenPosition(genomeTickNumber);
                binNUmber = axis.getBinNumberForGenomicPosition(genomePosition);

                x = (int) ((binNUmber - binOrigin) * context.getScaleFactor());

                String chrPosition = formatNumber((double) genomePosition / ts.getUnitMultiplier()) + " " + ts.getMajorUnit();
                int strWidth = g.getFontMetrics().stringWidth(chrPosition);
                int strPosition = isHorizontal() ? x - strWidth / 2 : -x - strWidth / 2;
                //if (strPosition > strEnd) {

                if (nTick % 2 == 0) {
                    g.drawString(chrPosition, strPosition, h - 15);
                }
                //strEnd = strPosition + strWidth;
                //}

                int xpos = (orientation == Orientation.HORIZONTAL ? x : -x);
                g.drawLine(xpos, h - 10, xpos, h - 2);
                nTick++;
            }
        }
    }

    public static String formatNumber(double position) {

        //NumberFormatter f = new NumberFormatter();
        DecimalFormat formatter = new DecimalFormat();
        return formatter.format((int) position);
        //return f.valueToString(position);

    }

    public static TickSpacing findSpacing(long maxValue, boolean scaleInKB) {

        if (maxValue < 10) {
            System.out.println("max value = " + maxValue);
            return new TickSpacing(1, "bp", 1);
        }


        // Now man zeroes?
        int nZeroes = (int) Math.log10(maxValue);
        String majorUnit = scaleInKB ? "kb" : "bp";
        int unitMultiplier = 1;
        if (nZeroes > 9) {
            majorUnit = scaleInKB ? "tb" : "gb";
            unitMultiplier = 1000000000;
        }
        if (nZeroes > 6) {
            majorUnit = scaleInKB ? "gb" : "mb";
            unitMultiplier = 1000000;
        } else if (nZeroes > 3) {
            majorUnit = scaleInKB ? "mb" : "kb";
            unitMultiplier = 1000;
        }

        double nMajorTicks = maxValue / Math.pow(10, nZeroes - 1);
        if (nMajorTicks < 25) {
            return new TickSpacing(Math.pow(10, nZeroes - 1), majorUnit, unitMultiplier);
        } else {
            return new TickSpacing(Math.pow(10, nZeroes) / 2, majorUnit, unitMultiplier);
        }
    }


    public static class TickSpacing {

        private double majorTick;
        private double minorTick;
        private String majorUnit = "";
        private int unitMultiplier = 1;

        TickSpacing(double majorTick, String majorUnit, int unitMultiplier) {
            this.majorTick = majorTick;
            this.minorTick = majorTick / 10.0;
            this.majorUnit = majorUnit;
            this.unitMultiplier = unitMultiplier;
        }

        public double getMajorTick() {
            return majorTick;
        }

        public double getMinorTick() {
            return minorTick;
        }

        public String getMajorUnit() {
            return majorUnit;
        }

        public void setMajorUnit(String majorUnit) {
            this.majorUnit = majorUnit;
        }

        public int getUnitMultiplier() {
            return unitMultiplier;
        }

        public void setUnitMultiplier(int unitMultiplier) {
            this.unitMultiplier = unitMultiplier;
        }
    }

// TODO -- possibly generalize?

    class ClickLink {

        Rectangle region;
        String value;
        String tooltipText;

        ClickLink(Rectangle region, String value, String tooltipText) {
            this.region = region;
            this.value = value;
            this.tooltipText = tooltipText;
        }
    }
}
