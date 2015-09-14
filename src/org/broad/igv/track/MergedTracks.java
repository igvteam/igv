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

import org.broad.igv.feature.LocusScore;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.renderer.DataRange;
import org.broad.igv.session.IGVSessionReader;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.TrackPanel;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Track to serve as a container for several tracks rendered on top of each other
 *
 * @author jacob
 * @date 2013-Nov-05
 */
@XmlType(factoryMethod = "getNextTrack")
public class MergedTracks extends DataTrack{

    @XmlAttribute
    protected Class clazz = MergedTracks.class;

    private Collection<DataTrack> memberTracks;

    /**
     * This is a session tag, changing it will break backwards compatibility
     * with sessions!
     */
    public static final String MEMBER_TRACK_TAG_NAME = "Track";

    public MergedTracks(String id, String name, Collection<DataTrack> inputTracks){
        super(null, id, name);
        initTrackList(inputTracks);
        this.autoScale = this.getAutoScale();
        setTrackAlphas(120);
    }

    private void initTrackList(Collection<DataTrack> inputTracks) {
        this.memberTracks = new ArrayList<DataTrack>(inputTracks.size());
        for(DataTrack inputTrack: inputTracks){
            if(inputTrack instanceof MergedTracks){
                this.memberTracks.addAll(((MergedTracks) inputTrack).getMemberTracks());
            }else{
                this.memberTracks.add(inputTrack);
            }
        }
    }

    @Override
    public Collection<ResourceLocator> getResourceLocators(){
        Collection<ResourceLocator> locators = new ArrayList<ResourceLocator>(memberTracks.size());
        for(DataTrack memTrack: memberTracks){
            locators.addAll(memTrack.getResourceLocators());
        }
        return locators;
    }

    @XmlElement(name = MEMBER_TRACK_TAG_NAME)
    Collection<DataTrack> getMemberTracks(){
        return this.memberTracks;
    }

    public void setTrackAlphas(int alpha) {
        for(Track track: memberTracks){
            track.setColor(ColorUtilities.modifyAlpha(track.getColor(), alpha));
            track.setAltColor(ColorUtilities.modifyAlpha(track.getAltColor(), alpha));
        }
    }


    @Override
    public void render(RenderContext context, Rectangle rect) {
        resetLastY();
        for(Track track: memberTracks){
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
        for(Track track: memberTracks){
            height = Math.max(height, track.getHeight());
        }
        return height;
    }

    @Override
    public void setHeight(int height) {
        super.setHeight(height);
        for(Track track: memberTracks){
            track.setHeight(height);
        }
    }

    @Override
    public void setDataRange(DataRange axisDefinition) {
        super.setDataRange(axisDefinition);
        for(Track track: memberTracks){
            track.setDataRange(axisDefinition);
        }
    }

    @Override
    public DataRange getDataRange() {
        if(this.dataRange == null){
            this.dataRange = DataRange.getFromTracks(memberTracks);
        }
        return this.dataRange;
    }

    @Override
    public String getValueStringAt(String chr, double position, int y, ReferenceFrame frame) {
        StringBuilder builder = new StringBuilder(memberTracks.size() + 2);
        builder.append(getName());
        builder.append("<br/>--------------<br/>");
        for(Track track: memberTracks){
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
        for(Track track: memberTracks){
            track.setRendererClass(rc);
        }
    }

    @Override
    public boolean getAutoScale() {
        boolean autoScale = super.getAutoScale();
        for(Track track: memberTracks){
            autoScale &= track.getAutoScale();
        }
        return autoScale;
    }

    @Override
    public void setAutoScale(boolean autoScale) {
        super.setAutoScale(autoScale);
        for(Track track: memberTracks){
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
        for(Track track: memberTracks){
            track.setColorScale(colorScale);
        }
    }

    @Override
    public void setWindowFunction(WindowFunction type) {
        super.setWindowFunction(type);
        for(Track track: memberTracks){
            track.setWindowFunction(type);
        }
    }

    @Override
    public void setShowDataRange(boolean showDataRange) {
        super.setShowDataRange(showDataRange);
        for(DataTrack track: memberTracks){
            track.setShowDataRange(showDataRange);
        }
    }

    @Override
    public IGVPopupMenu getPopupMenu(TrackClickEvent te) {

        IGVPopupMenu menu = new IGVPopupMenu();

        final List<Track> selfAsList = Arrays.asList((Track) this);
        menu.add(TrackMenuUtils.getTrackRenameItem(selfAsList));

        //Give users the ability to set the color of each track individually
        JMenu setPosColorMenu = new JMenu("Change Track Color (Positive Values)");
        JMenu setNegColorMenu = new JMenu("Change Track Color (Negative Values)");
        for(DataTrack track: memberTracks){

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

        menu.add(TrackMenuUtils.getRemoveMenuItem(selfAsList));


        return menu;
    }

    //Not clear whether this is necessary, don't think it is, the real problem is loading not saving
//    public void marshalMemberTracks(Marshaller m, Element trackElement) throws JAXBException{
//        for(DataTrack memberTrack: memberTracks){
//            JAXBElement element = new JAXBElement<DataTrack>(new QName("", DataTrack.DATA_TRACK), DataTrack.class, memberTrack);
//            m.marshal(element, trackElement);
//        }
//    }

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

    @SubtlyImportant
    private MergedTracks(){
        super(null, null, null);
    }

    @SubtlyImportant
    private static MergedTracks getNextTrack(){
        MergedTracks out = (MergedTracks) IGVSessionReader.getNextTrack();
        if (out == null){
            out = new MergedTracks();
        }
        return out;
    }
}
