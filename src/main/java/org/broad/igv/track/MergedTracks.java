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
import org.broad.igv.session.SessionElement;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
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

public class MergedTracks extends DataTrack implements ScalableTrack {

    private static double DEFAULT_ALPHA = 0.5;

    private Collection<DataTrack> memberTracks;
    private double alpha;


    public MergedTracks(String id, String name, Collection<DataTrack> inputTracks) {
        super(null, id, name);
        initTrackList(inputTracks);
        this.autoScale = this.getAutoScale();
        setTrackAlphas(DEFAULT_ALPHA);
    }

    public MergedTracks() {
        this.alpha = DEFAULT_ALPHA;
    }

    public void setMemberTracks(Collection<DataTrack> inputTracks) {
        initTrackList(inputTracks);
    }

    private void initTrackList(Collection<DataTrack> inputTracks) {
        this.memberTracks = new ArrayList<DataTrack>(inputTracks.size());
        for (DataTrack inputTrack : inputTracks) {
            if (inputTrack instanceof MergedTracks) {
                this.memberTracks.addAll(((MergedTracks) inputTrack).getMemberTracks());
            } else {
                this.memberTracks.add(inputTrack);
            }
        }

        // Set the group autoscale attribute only IF all tracks are in the same group
        if (memberTracks != null && memberTracks.size() > 1) {
            String group = memberTracks.iterator().next().getAttributeValue(AttributeManager.GROUP_AUTOSCALE);
            if (group != null) {
                for (Track t : memberTracks) {
                    if (!group.equals(t.getAttributeValue(AttributeManager.GROUP_AUTOSCALE))) {
                        this.removeAttribute(AttributeManager.GROUP_AUTOSCALE);
                        group = null;
                        break;
                    }
                }
                if (group != null) {
                    this.setAttributeValue(AttributeManager.GROUP_AUTOSCALE, group);
                }
            }
        }


        setTrackAlphas(alpha);
    }

    @Override
    public boolean isReadyToPaint(ReferenceFrame frame) {
        return this.memberTracks.stream().allMatch((t) -> t.isReadyToPaint(frame));
    }


    @Override
    public synchronized void load(ReferenceFrame referenceFrame) {

        for (DataTrack t : memberTracks) {
            if (!t.isReadyToPaint(referenceFrame)) {
                t.load(referenceFrame);
            }
        }
    }

    @Override
    public Collection<ResourceLocator> getResourceLocators() {
        Collection<ResourceLocator> locators = new ArrayList<ResourceLocator>(memberTracks.size());
        for (DataTrack memTrack : memberTracks) {
            locators.addAll(memTrack.getResourceLocators());
        }
        return locators;
    }

    public Collection<DataTrack> getMemberTracks() {
        return this.memberTracks;
    }

    public void setTrackAlphas(double alpha) {
        this.alpha = alpha;
        if (memberTracks != null) {
            int iAlpha = (int) Math.floor(alpha * 255);
            for (Track track : memberTracks) {
                track.setColor(ColorUtilities.modifyAlpha(track.getColor(), iAlpha));
                track.setAltColor(ColorUtilities.modifyAlpha(track.getAltColor(), iAlpha));
            }
        }
    }

    public double getTrackAlpha() {
        return alpha;
    }


    @Override
    public void render(RenderContext context, Rectangle rect) {
        for (Track track : memberTracks) {
            track.render(context, rect);
        }
    }

    @Override
    public int getHeight() {
        int height = super.getHeight();
        for (Track track : memberTracks) {
            height = Math.max(height, track.getHeight());
        }
        return height;
    }

    @Override
    public void setHeight(int height) {
        super.setHeight(height);
        for (Track track : memberTracks) {
            track.setHeight(height);
        }
    }

    @Override
    public void setDataRange(DataRange axisDefinition) {
        super.setDataRange(axisDefinition);
        for (Track track : memberTracks) {
            track.setDataRange(axisDefinition);
        }
    }

    @Override
    public DataRange getDataRange() {
        if (this.dataRange == null) {
            this.dataRange = DataRange.getFromTracks(memberTracks);
        }
        return this.dataRange;
    }

    @Override
    public String getValueStringAt(String chr, double position, int mouseX, int mouseY, ReferenceFrame frame) {
        StringBuilder builder = new StringBuilder(memberTracks.size() + 2);
        builder.append(getName());
        builder.append("<br/>--------------<br/>");
        for (Track track : memberTracks) {
            String curS = track.getValueStringAt(chr, position, mouseX, mouseY, frame);
            if (curS != null) {
                builder.append(curS);
                builder.append("<br/>--------------<br/>");
            }
        }
        return builder.toString();
    }

    @Override
    public LoadedDataInterval<List<LocusScore>> getSummaryScores(String chr, int startLocation, int endLocation, int zoom) {
        return null;
    }

    @Override
    public void setRendererClass(Class rc) {
        super.setRendererClass(rc);
        for (Track track : memberTracks) {
            track.setRendererClass(rc);
        }
    }

    @Override
    public boolean getAutoScale() {
        boolean autoScale = true;
        for (Track track : memberTracks) {
            autoScale &= track.getAutoScale();
        }
        return autoScale;
    }

    @Override
    public void setAutoScale(boolean autoScale) {
        for (Track track : memberTracks) {
            track.setAutoScale(autoScale);
        }
    }

    @Override
    public void setAttributeValue(String name, String value) {
        super.setAttributeValue(name, value);
        if (name.equals(AttributeManager.GROUP_AUTOSCALE) && memberTracks != null) {
            for (Track track : memberTracks) {
                track.setAttributeValue(name, value);
            }
        }
    }

    @Override
    public void removeAttribute(String name) {
        super.removeAttribute(name);
        if (name.equals(AttributeManager.GROUP_AUTOSCALE)) {
            for (Track track : memberTracks) {
                track.removeAttribute(name);
            }
        }
    }

    @Override
    public String getAttributeValue(String attributeName) {
        if (super.getAttributeValue(attributeName) != null) {
            return super.getAttributeValue(attributeName);
        } else if (memberTracks.size() > 0) {
            String value = null;
            for (Track track : memberTracks) {
                String nm = track.getAttributeValue(attributeName);
                if (nm != null) {
                    if (value == null) {
                        value = nm;
                    } else if (!value.equals(nm)) {
                        value += " / " + nm;
                    }
                }
            }
            return value;
        } else {
            // Should be impossible
            return null;
        }
    }

    @Override
    public ContinuousColorScale getColorScale() {
        return super.getColorScale();
    }

    @Override
    public void setColorScale(ContinuousColorScale colorScale) {
        super.setColorScale(colorScale);
        for (Track track : memberTracks) {
            track.setColorScale(colorScale);
        }
    }

    @Override
    public void setWindowFunction(WindowFunction type) {
        super.setWindowFunction(type);
        for (Track track : memberTracks) {
            track.setWindowFunction(type);
        }
    }

    @Override
    public void setShowDataRange(boolean showDataRange) {
        super.setShowDataRange(showDataRange);
        for (DataTrack track : memberTracks) {
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
        for (DataTrack track : memberTracks) {

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
        TrackMenuUtils.addDataItems(menu, selfAsList, true);
        for (Component c : menu.getComponents()) {
            if (c instanceof JMenuItem) {
                String text = ((JMenuItem) c).getText();
                text = text != null ? text.toLowerCase() : "null";
                if (text.contains("heatmap")) {
                    c.setEnabled(false);
                }
            }

        }

        return menu;
    }

    @Override
    public Range getInViewRange(ReferenceFrame referenceFrame) {

        List<LocusScore> scores = new ArrayList<LocusScore>();
        for (DataTrack track : memberTracks) {
            scores.addAll(track.getInViewScores(referenceFrame));
        }

        if (scores.size() > 0) {
            float min = Float.MAX_VALUE;
            float max = -Float.MAX_VALUE;
            for (LocusScore score : scores) {
                float value = score.getScore();
                if (!Float.isNaN(value)) {
                    min = Math.min(value, min);
                    max = Math.max(value, max);
                }
            }
            return new Range(min, max);
        } else {
            return null;
        }
    }

    public void marshalXML(Document document, Element element) {

        super.marshalXML(document, element);

        if (alpha != DEFAULT_ALPHA) {
            element.setAttribute("alpha", String.valueOf(alpha));
        }

        for (DataTrack track : memberTracks) {
            Element trackElement = document.createElement(SessionElement.TRACK);
            track.marshalXML(document, trackElement);
            element.appendChild(trackElement);
        }

    }

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);

        if (element.hasAttribute("alpha")) {
            this.alpha = Double.valueOf(element.getAttribute("alpha"));
        }
        // Un-marshalling handled in IGVSessionReader

    }


    private enum ChangeTrackMethod {
        POSITIVE, NEGATIVE
    }

    private static class ChangeTrackColorActionListener implements ActionListener {

        private Track mTrack;
        private ChangeTrackMethod method;

        private ChangeTrackColorActionListener(Track track, ChangeTrackMethod method) {
            this.mTrack = track;
            this.method = method;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            switch (this.method) {
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
    private static class ColorIcon implements Icon {

        private Color color;
        private int iconSize;

        ColorIcon(Color color) {
            this(color, 16);
        }

        ColorIcon(Color color, int iconSize) {
            this.color = color;
            this.iconSize = iconSize;
        }

        @Override
        public void paintIcon(Component c, Graphics g, int x, int y) {
            Graphics cg = g.create();
            cg.setColor(this.color);
            if (this.iconSize > c.getHeight()) {
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


    public JDialog getAlphaDialog() {

        ChangeListener changeListener = new ChangeListener() {
            public void stateChanged(ChangeEvent changeEvent) {
                JSlider theSlider = (JSlider) changeEvent.getSource();
                double alpha = theSlider.getValue() / 100.0;
                MergedTracks.this.setTrackAlphas(alpha);
                IGV.getInstance().getContentPane().repaint();
            }
        };

        JOptionPane optionPane = new JOptionPane();

        int initialValue = (int) Math.floor(MergedTracks.this.getTrackAlpha() * 100);
        JSlider slider = new JSlider(0, 100, initialValue);
        slider.setMajorTickSpacing(10);
        slider.setPaintTicks(false);
        slider.setPaintLabels(false);
        slider.addChangeListener(changeListener);

        optionPane.setMessage(new Object[]{"Adjust transparency: ", slider});
        optionPane.setMessageType(JOptionPane.QUESTION_MESSAGE);
        optionPane.setOptionType(JOptionPane.OK_CANCEL_OPTION);
        JDialog dialog = optionPane.createDialog(IGV.getMainFrame(), "Transparency");
        return dialog;
    }


}
