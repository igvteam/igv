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

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.broad.igv.sam;

import org.broad.igv.track.RenderContext;

import java.awt.*;
import java.util.List;
import java.util.Map;

/**
 * @author jrobinso
 */
public interface FeatureRenderer {

    /**
     * @param alignments
     * @param context
     * @param rowRect           Rectangle in which alignments should be drawn
     * @param trackRect         Rectangle of entire alignment track. Useful if drawing vertical features.
     * @param renderOptions
     * @param leaveMargin
     * @param selectedReadNames
     */
    public void renderAlignments(List<Alignment> alignments, RenderContext context,
                                 Rectangle rowRect, Rectangle trackRect, AlignmentTrack.RenderOptions renderOptions,
                                 boolean leaveMargin, Map<String, Color> selectedReadNames);
}
