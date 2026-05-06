package org.igv.track;

import org.igv.feature.LocusScore;
import org.igv.renderer.ContinuousColorScale;
import org.igv.renderer.DataRange;
import org.igv.session.SessionElement;
import org.igv.ui.IGV;
import org.igv.ui.action.OverlayTracksMenuAction;
import org.igv.ui.color.ColorUtilities;
import org.igv.ui.panel.IGVPopupMenu;
import org.igv.ui.panel.ReferenceFrame;
import org.igv.util.ResourceLocator;
import org.json.JSONArray;
import org.json.JSONObject;
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

    private List<DataTrack> memberTracks;
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

    @Override
    public TrackType getType() {
        return TrackType.merged;
    }

    public void setMemberTracks(Collection<DataTrack> inputTracks) {
        initTrackList(inputTracks);
    }

    private void initTrackList(Collection<DataTrack> inputTracks) {
        this.memberTracks = new ArrayList<>(inputTracks.size());
        for (DataTrack inputTrack : inputTracks) {
            if (inputTrack instanceof MergedTracks) {
                this.memberTracks.addAll(((MergedTracks) inputTrack).getMemberTracks());
            } else {
                this.memberTracks.add(inputTrack);
            }
        }

        // Set the height all tracks to the same value as the merged track.
        for (Track track : memberTracks) {
            track.setHeight(this.getHeight());
        }

        // Set the group autoscale attribute only IF all tracks are in the same group
        if (memberTracks.size() > 1) {
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

    public List<DataTrack> getMemberTracks() {
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
    public void render(RenderContext context) {
        for (Track track : memberTracks) {
            track.render(context);
        }
    }

    @Override
    public int getContentHeight() {
        int height = super.getContentHeight();
        for (Track track : memberTracks) {
            height = Math.max(height, track.getContentHeight());
        }
        return height;
    }

    @Override
    public void setHeight(int height) {
        super.setHeight(height);
        if(memberTracks != null) {
            for (Track track : memberTracks) {
                track.setHeight(height);
            }
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
    public List<Component> getPopupMenuItems(TrackClickEvent te) {

        List<Component> items = new ArrayList<>();

        final List<Track> selfAsList = Arrays.asList((Track) this);
        items.add(TrackMenuUtils.getTrackRenameItem(selfAsList));

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
        items.add(setPosColorMenu);
        items.add(setNegColorMenu);

        items.add(TrackMenuUtils.getChangeTrackHeightItem(selfAsList));

        items.add(new JPopupMenu.Separator());
        JMenuItem alphaItem = new JMenuItem("Adjust Transparency");
        alphaItem.addActionListener(e -> {
            JDialog alphaDialog = getAlphaDialog();
            alphaDialog.setVisible(true);
        });
        items.add(alphaItem);

        items.add(new JPopupMenu.Separator());
        JMenuItem unmergeItem = new JMenuItem("Separate Tracks");

        unmergeItem.addActionListener(e -> {
            long order = getOrder();
            setTrackAlphas(1.0);
            // Set the order of member tracks to match the merged track's order
            for (Track memberTrack : getMemberTracks()) {
                memberTrack.setOrder(order);
            }
            IGV.getInstance().deleteTracks(List.of(this));
            IGV.getInstance().addTracks(new ArrayList<>(getMemberTracks()));
            IGV.getInstance().repaint();
        });
        items.add(unmergeItem);

        return items;
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

    @Override
    public void unmarshalXML(Element element, Integer version) {

        super.unmarshalXML(element, version);
        if (element.hasAttribute("alpha")) {
            this.alpha = Double.valueOf(element.getAttribute("alpha"));
        }
        // Un-marshalling handled in IGVSessionReader
    }

    @Override
    public void marshalJSON(JSONObject json) {
        super.marshalJSON(json);
        json.put("alpha", this.alpha);

        JSONArray tracksArray = new JSONArray();
        for (DataTrack track : memberTracks) {
            JSONObject trackJson = new JSONObject();
            track.marshalJSON(trackJson);
            tracksArray.put(trackJson);
        }
        json.put("tracks", tracksArray);
    }

    @Override
    public void unmarshalJSON(JSONObject jsonObject) {
        super.unmarshalJSON(jsonObject);
        if (jsonObject.has("alpha")) {
            this.alpha = jsonObject.getDouble("alpha");
        }
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
        JDialog dialog = optionPane.createDialog(IGV.getInstance().getMainFrame(), "Transparency");
        return dialog;
    }


}
