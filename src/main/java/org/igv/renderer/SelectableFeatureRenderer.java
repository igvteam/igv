package org.igv.renderer;

import com.google.common.base.Predicate;
import org.igv.feature.Exon;
import org.igv.feature.IExon;
import org.igv.feature.IGVFeature;
import org.igv.track.RenderContext;
import org.igv.track.Track;
import org.igv.util.collections.CollUtils;

import java.awt.*;
import java.util.Collections;
import java.util.List;
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
        for(IExon selEx: selectedExons){
            if(selEx.contains(exon.getStart()) || selEx.contains(exon.getEnd())){
                isSelected = true;
                break;
            }
        }

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
