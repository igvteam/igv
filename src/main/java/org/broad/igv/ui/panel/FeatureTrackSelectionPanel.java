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

package org.broad.igv.ui.panel;

import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import java.awt.*;
import java.util.Collection;

/**
 * Created with IntelliJ IDEA.
 * User: jrobinso
 * Date: 6/14/12
 * Time: 8:34 AM
 * To change this template use File | Settings | File Templates.
 */
public class FeatureTrackSelectionPanel extends JPanel {


    public FeatureTrackSelectionPanel() {
        init();
    }

    void init() {

        if(IGV.hasInstance()) {

            LayoutManager lm = new BoxLayout(this, BoxLayout.Y_AXIS);
            this.setLayout(lm);

            Collection<Track> tracks = IGV.getInstance().getAllTracks();
            for(Track t : tracks) {
                if(t instanceof FeatureTrack) {

                    JCheckBox cb = new JCheckBox(t.getName());
                    this.add(cb);


                }
            }
        }

    }

}
