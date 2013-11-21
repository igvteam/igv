/*
 * Copyright (c) 2007-2013 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.track;


import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.TrackPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Track to serve as a container for several tracks rendered on top of each other
 *
 * @author jacob
 * @date 2013-Nov-05
 */
public class MergedTracks extends DataTrack{

    private Collection<DataTrack> trackList;

    public MergedTracks(String id, String name, Collection<DataTrack> trackList){
        super(null, id, name);
        this.trackList = trackList;
        this.autoScale = this.getAutoScale();
        setTrackAlphas(150);
    }

    private void setTrackAlphas(int alpha) {
        for(Track track: trackList){
            track.setColor(ColorUtilities.modifyAlpha(track.getColor(), alpha));
            track.setAltColor(ColorUtilities.modifyAlpha(track.getAltColor(), alpha));
        }
    }


    @Override
    public void render(RenderContext context, Rectangle rect) {
        resetLastY();
        for(Track track: trackList){
            if(isRepeatY(rect)){
                track.overlay(context, rect);
            }else{
                track.render(context, rect);
                lastRenderY = rect.y;
            }

        }
    }

    @Override
    public int getHeight() {
        int height = super.getHeight();
        for(Track track: trackList){
            height = Math.max(height, track.getHeight());
        }
        return height;
    }

    @Override
    public void setHeight(int height) {
        super.setHeight(height);
        for(Track track: trackList){
            track.setHeight(height);
        }
    }

    @Override
    public void setDataRange(DataRange axisDefinition) {
        super.setDataRange(axisDefinition);
        for(Track track: trackList){
            track.setDataRange(axisDefinition);
        }
    }

    @Override
    public DataRange getDataRange() {
        if(this.dataRange == null){
            this.dataRange = DataRange.getFromTracks(trackList);
        }
        return this.dataRange;
    }

    @Override
    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        StringBuilder builder = new StringBuilder(trackList.size() + 2);
        builder.append(getName());
        builder.append("<br/>--------------<br/>");
        for(Track track: trackList){
            String curS = track.getValueStringAt(chr, position, y, frame);
            if (curS != null) {
                builder.append(curS);
                builder.append("<br/>--------------<br/>");
            }
        }
        return builder.toString();
    }

    @Override
    public List<LocusScore> getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        return null;
    }

    @Override
    public void setRendererClass(Class rc) {
        super.setRendererClass(rc);
        for(Track track: trackList){
            track.setRendererClass(rc);
        }
    }

    @Override
    public boolean getAutoScale() {
        boolean autoScale = super.getAutoScale();
        for(Track track: trackList){
            autoScale &= track.getAutoScale();
        }
        return autoScale;
    }

    @Override
    public void setAutoScale(boolean autoScale) {
        super.setAutoScale(autoScale);
        for(Track track: trackList){
            track.setAutoScale(autoScale);
        }
    }

    @Override
    public ContinuousColorScale getColorScale() {
        return super.getColorScale();
    }

    @Override
    public void setColorScale(ContinuousColorScale colorScale) {
        super.setColorScale(colorScale);
        for(Track track: trackList){
            track.setColorScale(colorScale);
        }
    }

    @Override
    public void setWindowFunction(WindowFunction type) {
        super.setWindowFunction(type);
        for(Track track: trackList){
            track.setWindowFunction(type);
        }
    }

    @Override
    public void setShowDataRange(boolean showDataRange) {
        super.setShowDataRange(showDataRange);
        for(DataTrack track: trackList){
            track.setShowDataRange(showDataRange);
        }
    }

    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        List<Track> selfAsList = Arrays.asList((Track) this);
        menu.add(TrackMenuUtils.getTrackRenameItem(selfAsList));

        //Give users the ability to set the color of each track individually
        JMenu setPosColorMenu = new JMenu("Change Track Color (Positive Values)");
        JMenu setNegColorMenu = new JMenu("Change Track Color (Negative Values)");
        for(DataTrack track: trackList){

            Icon posColorIcon = new ColorIcon(track.getColor());
            JMenuItem posItem = new JMenuItem(track.getName(), posColorIcon);
            posItem.addActionListener(new ChangeTrackColorActionListener(track, ChangeTrackMethod.POSITIVE));
            setPosColorMenu.add(posItem);

            Icon negColorIcon = new ColorIcon(track.getAltColor());
            JMenuItem negItem = new JMenuItem(track.getName(), negColorIcon);
            negItem.addActionListener(new ChangeTrackColorActionListener(track, ChangeTrackMethod.NEGATIVE));
            setNegColorMenu.add(negItem);
        }
        menu.add(setPosColorMenu);
        menu.add(setNegColorMenu);

        menu.add(TrackMenuUtils.getChangeTrackHeightItem(selfAsList));
        menu.add(TrackMenuUtils.getChangeFontSizeItem(selfAsList));

        menu.addSeparator();
        TrackMenuUtils.addDataItems(menu, selfAsList);
        for(Component c: menu.getComponents()){
            if(c instanceof JMenuItem){
                String text = ((JMenuItem) c).getText();
                text = text != null ? text.toLowerCase() : "null";
                if(text.contains("heatmap")){
                    c.setEnabled(false);
                }    
            }

        }

        menu.addSeparator();

        JMenuItem unmergeItem = new JMenuItem("Separate Tracks");
        unmergeItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                TrackPanel panel = TrackPanel.getParentPanel(MergedTracks.this);
                setTrackAlphas(255);
                panel.addTracks(trackList);
                IGV.getInstance().removeTracks(Arrays.<Track>asList(MergedTracks.this));
            }
        });
        menu.add(unmergeItem);

        menu.add(TrackMenuUtils.getRemoveMenuItem(selfAsList));


        return menu;
    }

    private enum ChangeTrackMethod{
        POSITIVE, NEGATIVE
    }

    private static class ChangeTrackColorActionListener implements ActionListener{

        private Track mTrack;
        private ChangeTrackMethod method;

        private ChangeTrackColorActionListener(Track track, ChangeTrackMethod method){
            this.mTrack=track;
            this.method = method;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            switch(this.method){
                case POSITIVE:
                    TrackMenuUtils.changeTrackColor(Arrays.asList(this.mTrack));
                    break;
                case NEGATIVE:
                    TrackMenuUtils.changeAltTrackColor(Arrays.asList(this.mTrack));
                    break;
                default:
                    throw new IllegalStateException("Method not understood: " + this.method);
            }
        }
    }

    /**
     * A square, solid color Icon
     */
    private static class ColorIcon implements Icon{

        private Color color;
        private int iconSize;

        ColorIcon(Color color){
            this(color, 16);
        }

        ColorIcon(Color color, int iconSize){
            this.color = color;
            this.iconSize = iconSize;
        }

        @Override
        public void paintIcon(Component c, Graphics g, int x, int y) {
            Graphics cg = g.create();
            cg.setColor(this.color);
            if(this.iconSize > c.getHeight()){
                this.iconSize = c.getHeight();
            }
            cg.fillRect(x, y, this.iconSize, this.iconSize);
        }

        @Override
        public int getIconWidth() {
            return this.iconSize;
        }

        @Override
        public int getIconHeight() {
            return this.iconSize;
        }
    }
}
