/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.plugin.mongovariant;

import org.broad.igv.renderer.GraphicUtils;
import org.broad.igv.track.RenderContext;
import org.broad.igv.ui.FontManager;
import org.broad.igv.variant.Variant;
import org.broad.igv.variant.VariantRenderer;
import org.broad.igv.variant.VariantTrack;
import org.broadinstitute.sting.gatk.walkers.na12878kb.core.TruthStatus;

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
        return true;
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
            GraphicUtils.drawCenteredText(symbol, pixelX, bandY - bandHeight, xWidth, bandHeight, g);//Color.white);

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
