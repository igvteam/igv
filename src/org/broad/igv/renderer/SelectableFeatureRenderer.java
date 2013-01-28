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

package org.broad.igv.renderer;

import org.broad.igv.feature.Exon;
import org.broad.igv.feature.IExon;

import java.awt.*;
import java.util.Collections;
import java.util.Set;

/**
 * User: jacob
 * Date: 2013-Jan-28
 */
public class SelectableFeatureRenderer extends IGVFeatureRenderer {

    private static final Stroke lineStroke;
    private static final Color selectedBorderColor;
    private Set<IExon> selectedExons = Collections.emptySet();



    static{
        lineStroke = new BasicStroke(2.0f);
        selectedBorderColor = Color.blue;
    }

    public SelectableFeatureRenderer(){
        AA_COLOR_1 = new Color(AA_COLOR_1.getRed(), AA_COLOR_1.getGreen(), AA_COLOR_1.getBlue(), 120);
        AA_COLOR_2 = new Color(AA_COLOR_2.getRed(), AA_COLOR_2.getGreen(), AA_COLOR_2.getBlue(), 120);
    }

    @Override
    protected void drawExonRect(Graphics blockGraphics, Exon exon, int x, int y, int width, int height) {
        boolean isSelected = selectedExons.contains(exon);
//        for(IExon selEx: selectedExons){
//            if(selEx.contains(exon.getStart()) || selEx.contains(exon.getEnd())){
//                isSelected = true;
//                break;
//            }
//        }

        if(isSelected){
            blockGraphics.clearRect(x, y, width, height);

            Graphics2D borderGraphics = (Graphics2D) blockGraphics.create();
            borderGraphics.setColor(selectedBorderColor);
            borderGraphics.setStroke(lineStroke);
            borderGraphics.drawRect(x, y, width, height);
        }else{
            super.drawExonRect(blockGraphics, exon, x, y, width, height);
        }
    }

    public void setSelectedExons(Set<IExon> selectedExons) {
        this.selectedExons = selectedExons;
    }
}
