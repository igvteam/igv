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

package org.broad.igv.peaks;

import org.broad.igv.data.DataSource;
import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.BarChartRenderer;
import org.broad.igv.renderer.Renderer;
import org.broad.igv.track.RenderContext;
import org.broad.igv.track.Track;
import org.broad.igv.ui.color.ColorUtilities;

import java.awt.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Apr 23, 2011
 */
public class PeakRenderer implements Renderer<LocusScore> {

    BarChartRenderer chartRenderer = new BarChartRenderer();
    public static final int SIGNAL_CHART_HEIGHT = 15;

    public void render(List<LocusScore> peakList, RenderContext context, Rectangle rect, Track t) {


        PeakTrack track = (PeakTrack) t;

        final double locScale = context.getScale();
        final double origin = context.getOrigin();

        int lastPX = -1;
        double pXMin = rect.getMinX();
        double pXMax = rect.getMaxX();

        Color fgColor = track.getColor();
        Graphics2D borderGraphics = context.getGraphic2DForColor(Color.black);

        String chr = context.getChr();
        int contextStart = (int) context.getOrigin();
        int contextEnd = (int) context.getEndLocation();
        int zoom = context.getZoom();

        if (PeakTrack.isShowPeaks()) {
            int h = track.bandHeight;
            int peakHeight = PeakTrack.isShowSignals() ? track.peakHeight : h;

            for (LocusScore ls : peakList) {
                Peak peak = (Peak) ls;

                int start = peak.getStart();
                int end = peak.getEnd();
                int pX = (int) ((start - origin) / locScale);
                int dX = (int) Math.max(2, (end - start) / locScale);

                if (pX + dX < pXMin) continue;
                if (pX > pXMax) break;

                float score = peak.getCombinedScore();
                if (PeakTrack.getColorOption() == PeakTrack.ColorOption.FOLD_CHANGE) {
                    score = peak.getDynamicScore();
                }

                int top = rect.y + 1;
                if (PeakTrack.isShowSignals()) {
                    top += track.signalHeight;
                }

                if (PeakTrack.animate) {
                    int step = track.getTimeStep();
                    float[] timeScores = peak.getTimeScores();
                    int idx = Math.min(step, timeScores.length - 1);
                    score = timeScores[idx];
                    drawPeak(context, fgColor, pX, dX, top, peakHeight, score, PeakTrack.getColorOption());
                } else {
                    if (track.getDisplayMode() == Track.DisplayMode.EXPANDED) {
                        float[] timeScores = peak.getTimeScores();
                        for (int i = 0; i < timeScores.length; i++) {
                            score = timeScores[i];
                            drawPeak(context, fgColor, pX, dX, top, peakHeight, score, PeakTrack.getColorOption());
                            top += h;
                        }
                    } else {
                        drawPeak(context, fgColor, pX, dX, top, peakHeight, score, PeakTrack.getColorOption());

//                    float[] timeScores = peak.getTimeScores();
//                    if (dX > 10 && timeScores.length > 1) {
//                        Graphics2D gLine = context.getGraphic2DForColor(Color.black);
//                        gLine.setRenderingHint(RenderingHints.KEY_ANTIALIASING,PreferenceManager.getInstance().getAntiAliasingHint());
//                        double deltaX = ((double) dX / (timeScores.length - 1));
//                        double scaleY = peakHeight / 200.0;
//                        int bottomY = top + peakHeight;
//                        int lastX = pX;
//                        int lastY = (bottomY - (int) (scaleY * timeScores[0]));
//                        for (int i = 1; i < timeScores.length; i++) {
//                            score = timeScores[i];
//                            int x = pX + (int) (i * deltaX);
//                            int y = Math.max(top + 2, (int) (bottomY - (int) (scaleY * timeScores[i])));
//                            gLine.drawLine(lastX, lastY, x, y);
//                            lastX = x;
//                            lastY = y;
//
//                        }
//                    }
                    }
                }
            }
        }

        if (PeakTrack.isShowSignals() && track.getSignalPath() != null) {
            int h = track.bandHeight;
            int signalHeight = PeakTrack.isShowPeaks() ? track.signalHeight : h;

            if (PeakTrack.animate) {
                int step = track.getTimeStep();
                DataSource[] timeSignalSources = track.getTimeSignalSources();
                int idx = Math.min(step, timeSignalSources.length - 1);
                DataSource src = timeSignalSources[idx];
                if (src != null) {
                    int top = rect.y + 2;List<LocusScore> timeSignals = src.getSummaryScoresForRange(chr, contextStart, contextEnd, zoom);
                    Rectangle timeSignalRect = new Rectangle(rect.x, top, rect.width, signalHeight - 1);
                    chartRenderer.render(timeSignals, context, timeSignalRect, track);
                }
            } else {
                if (track.getDisplayMode() == Track.DisplayMode.EXPANDED) {

                    DataSource[] timeSignalSources = track.getTimeSignalSources();
                    if (timeSignalSources != null) {
                        int top = rect.y + 2;
                        for (int i = 0; i < timeSignalSources.length; i++) {
                            DataSource src = timeSignalSources[i];
                            if (src != null) {
                                List<LocusScore> timeSignals = src.getSummaryScoresForRange(chr, contextStart, contextEnd, zoom);
                                Rectangle timeSignalRect = new Rectangle(rect.x, top, rect.width, signalHeight - 1);
                                chartRenderer.render(timeSignals, context, timeSignalRect, track);
                            }
                            top += h;

                        }
                    }
                } else {
                    final PeakTrack.WrappedDataSource signalSource = track.getSignalSource(chr, contextStart, contextEnd, zoom);
                    if (signalSource != null) {
                        List<LocusScore> signals = signalSource.getSummaryScoresForRange(chr, contextStart, contextEnd, zoom);
                        Rectangle signalRect = new Rectangle(rect.x, rect.y + 1, rect.width, signalHeight - 1);
                        chartRenderer.render(signals, context, signalRect, track);
                    }
                }
            }
        }


        if (track.getDisplayMode() == Track.DisplayMode.EXPANDED) {
            borderGraphics.drawLine(rect.x, rect.y, rect.x + rect.width, rect.y);
            borderGraphics.drawLine(rect.x, rect.y + rect.height, rect.x + rect.width, rect.y + rect.height);
        }
    }


    private void drawPeak(RenderContext context, Color fgColor,
                          int pX, int dX, int top, int h, float score, PeakTrack.ColorOption option) {
        Color c = getColor(fgColor, score, option);
        Graphics2D g = context.getGraphic2DForColor(c);
        g.fillRect(pX, top + 1, dX, h - 2);

    }

    private Color getColor(Color fgColor, float score, PeakTrack.ColorOption option) {
        Color c = null;
        float alpha = 1.0f;
        if (option == PeakTrack.ColorOption.SCORE) {
            // scale is 1 -> 100
            int shadeStep = (int) (score / 10);
            alpha = Math.max(0.2f, (Math.min(1.0f, shadeStep * 0.1f)));

        } else if (option == PeakTrack.ColorOption.FOLD_CHANGE) {

            // Scale is 1.58 -> 3.3 log2
            if (Math.abs(score) < 1.58) {
                alpha = 1f;
                fgColor = Color.gray;
            } else {
                // vary alpha from .3 -> 1 over range 3.3 -> 1.58 in steps of 10
                int shadeStep = (int) ((Math.abs(score) - 1.58) / (3.3 - 1.58) * 10);
                alpha = Math.max(0.3f, (Math.min(1f, 0.3f + shadeStep * 0.1f)));
                fgColor = (score < 0 ? Color.blue : Color.red);
            }
        }

        c = ColorUtilities.getCompositeColor(fgColor, alpha);
        //}
        return c;
    }
}
