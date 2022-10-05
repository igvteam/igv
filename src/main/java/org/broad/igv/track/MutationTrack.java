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

package org.broad.igv.track;

import htsjdk.tribble.Feature;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.Mutation;
import org.broad.igv.prefs.Constants;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.TooltipTextFrame;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.util.UIUtilities;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.net.URL;

/**
 * @author Jim Robinson
 * @date 11/8/11
 */
public class MutationTrack extends FeatureTrack {

    private static Logger log = LogManager.getLogger(MutationTrack.class);

    public MutationTrack(ResourceLocator locator, String id, FeatureSource source) {
        super(locator, id, source);
        setSortable(true);
    }

    public MutationTrack() {
    }

    @Override
    public void overlay(RenderContext context, Rectangle rect) {
        if (!context.getChr().equals(Globals.CHR_ALL) ||
                IGV.getInstance().getSession().getPreferenceAsBoolean(Constants.OVERLAY_MUTATIONS_WHOLE_GENOME)) {
            renderFeatures(context, rect);
        }
    }


    @Override
    public boolean handleDataClick(TrackClickEvent te) {
        super.openTooltipWindow(te);
        return true;
    }


}
