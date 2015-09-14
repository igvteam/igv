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

package org.broad.igv.plugin.mongovariant;

import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.FontManager;
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.VariantRenderer;
import org.broad.igv.variant.VariantTrack;
import org.broadinstitute.gatk.tools.walkers.na12878kb.core.TruthStatus;

import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;

/**
 * User: jacob
 * Date: 2013-Jan-28
 */
public class VariantReviewRenderer extends VariantRenderer {

    private static HashMap<TruthStatus, String> symbolMap;
    private static final String REVIEWED_STRING = "*";

    private static Graphics2D stringGraphics;
    private static final int FONT_SIZE = 10;


    static{
        symbolMap = new HashMap<TruthStatus, String>(4);
        symbolMap.put(TruthStatus.TRUE_POSITIVE, "TP");
        symbolMap.put(TruthStatus.FALSE_POSITIVE, "FP");
        symbolMap.put(TruthStatus.UNKNOWN, "?");
        symbolMap.put(TruthStatus.SUSPECT, "S");
    }

    public VariantReviewRenderer(VariantTrack track) {
        super(track);
    }

    @Override
    protected boolean defaultUseAlpha() {
        return false;
    }

    @Override
    public void renderSiteBand(Variant variant, Rectangle bandRectangle, int pixelX, int xWidth, RenderContext context) {
        super.renderSiteBand(variant, bandRectangle, pixelX, xWidth, context);

        if(variant instanceof MongoVCFVariant){
            MongoVCFVariant mvc = (MongoVCFVariant) variant;
            TruthStatus truthStatus = mvc.getTruthStatus();

            int bandY = calculateBottomYSiteBand(bandRectangle);
            int bandHeight = calculateBarHeightSiteBand(bandRectangle);

            Graphics2D g = context.getGraphic2DForColor(Color.black);
            g.setFont(FontManager.getFont(FONT_SIZE).deriveFont(Font.BOLD));

            String symbol = symbolMap.get(truthStatus);
            symbol = symbol != null ? symbol : "";
            GraphicUtils.drawCenteredText(symbol, pixelX, bandY - bandHeight, xWidth, bandHeight, g, Color.white);

            if(mvc.isReviewed()){
                g = context.getGraphic2DForColor(Color.red);
                g.setFont(FontManager.getFont(FONT_SIZE + 4).deriveFont(Font.BOLD));
                FontMetrics fontMetrics = g.getFontMetrics();
                Rectangle2D textBounds = fontMetrics.getStringBounds(REVIEWED_STRING, g);
                g.drawString(REVIEWED_STRING, pixelX - (int) (textBounds.getWidth() / 2), bandY - bandHeight + (int) (textBounds.getHeight() / 2));
            }
        }
    }

}
