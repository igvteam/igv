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

import com.jidesoft.swing.JideSplitPane;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.lang.reflect.InvocationTargetException;
import java.util.List;
import java.util.*;

import static org.broad.igv.prefs.Constants.*;

/**
 * @author jrobinso
 * @date Sep 10, 2010
 */
public class MainPanel extends JPanel implements Paintable {

    private static Logger log = LogManager.getLogger(MainPanel.class);

    IGV igv;

    // private static final int DEFAULT_NAME_PANEL_WIDTH = 160;

    private int namePanelX;
    private int namePanelWidth = PreferencesManager.getPreferences().getAsInt(NAME_PANEL_WIDTH);
    private int attributePanelX;
    private int attributePanelWidth;
    private int dataPanelX;
    private int dataPanelWidth;


    public IGVPanel applicationHeaderPanel;
    public HeaderPanelContainer headerPanelContainer;
    private TrackPanelScrollPane dataTrackScrollPane;
    private TrackPanelScrollPane featureTrackScrollPane;
    private JideSplitPane centerSplitPane;
    private NameHeaderPanel nameHeaderPanel;
    private AttributeHeaderPanel attributeHeaderPanel;

    private int hgap = 5;
    private JScrollPane headerScrollPane;


    public MainPanel(IGV igv) {
        this.igv = igv;
        initComponents();

        //Load IGV logo
//        try {
//            BufferedImage logo = ImageIO.read(getClass().getResource("resources/IGV_64.png"));
//            JLabel picLabel = new JLabel(new ImageIcon(logo));
//            picLabel.setVerticalAlignment(SwingConstants.CENTER);
//            nameHeaderPanel.add(picLabel);
//        } catch (IOException e) {
//            //pass
//        }

        addComponentListener(new ComponentListener() {

            public void componentResized(ComponentEvent componentEvent) {
                revalidateTrackPanels();
                igv.repaint();
            }

            public void componentMoved(ComponentEvent componentEvent) {
                //To change body of implemented methods use File | Settings | File Templates.
            }

            public void componentShown(ComponentEvent componentEvent) {
                //To change body of implemented methods use File | Settings | File Templates.
            }

            public void componentHidden(ComponentEvent componentEvent) {
                //To change body of implemented methods use File | Settings | File Templates.
            }
        });
    }

    public void setDividerFractions(double[] fractions) {
        int[] dividerLocations = new int[fractions.length];
        double h = centerSplitPane.getHeight();
        for (int i = 0; i < fractions.length; i++) {
            dividerLocations[i] = (int) Math.round(h * fractions[i]);
        }
        centerSplitPane.setDividerLocations(dividerLocations);
    }

    public double[] getDividerFractions() {
        int[] dividerLocations = centerSplitPane.getDividerLocations();
        double h = centerSplitPane.getHeight();
        double[] dividerFractions = new double[dividerLocations.length];
        for (int i = 0; i < dividerLocations.length; i++) {
            dividerFractions[i] = dividerLocations[i] / h;
        }
        return dividerFractions;
    }


    public void collapseNamePanel() {
        namePanelWidth = 0;
        revalidateTrackPanels();
    }

    public void expandNamePanel() {
        namePanelWidth = PreferencesManager.getPreferences().getAsInt(NAME_PANEL_WIDTH);
        revalidateTrackPanels();
    }

    public void setNamePanelWidth(int width) {
        this.namePanelWidth = width;
        revalidateTrackPanels();
    }

    public void revalidateTrackPanels() {
        updatePanelDimensions();
        UIUtilities.invokeOnEventThread(() -> {
            this.applicationHeaderPanel.invalidate();
            for (TrackPanel tp : this.getTrackPanels()) {
                tp.invalidate();
            }
            this.invalidate(); // this should not be neccessary, but is harmless
            this.validate();
        });
    }

    public void removeHeader() {
        remove(headerScrollPane);
        revalidate();
    }

    public void restoreHeader() {
        add(headerScrollPane, BorderLayout.NORTH);
        revalidate();
    }


    @Override
    public void doLayout() {
        super.doLayout();
        applicationHeaderPanel.doLayout();
        for (TrackPanel tp : getTrackPanels()) {
            tp.getScrollPane().doLayout();
        }
    }

    @Override
    public void setBackground(Color color) {
        super.setBackground(color);
        if (headerPanelContainer != null) {
            applicationHeaderPanel.setBackground(color);
            nameHeaderPanel.setBackground(color);
            attributeHeaderPanel.setBackground(color);
            headerPanelContainer.setBackground(color);
            nameHeaderPanel.setBackground(color);
            attributeHeaderPanel.setBackground(color);
            for (TrackPanel tsp : getTrackPanels()) {
                tsp.setBackground(color);
            }
        }

    }

    private void initComponents() {

        setPreferredSize(new java.awt.Dimension(1021, 510));
        setLayout(new java.awt.BorderLayout());

        nameHeaderPanel = new NameHeaderPanel();
        nameHeaderPanel.setBackground(new java.awt.Color(255, 255, 255));
        nameHeaderPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        nameHeaderPanel.setPreferredSize(new java.awt.Dimension(0, 0));

        attributeHeaderPanel = new AttributeHeaderPanel();
        attributeHeaderPanel.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        attributeHeaderPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        attributeHeaderPanel.setPreferredSize(new java.awt.Dimension(0, 0));


        headerPanelContainer = new HeaderPanelContainer();
        headerScrollPane = new JScrollPane();
        headerScrollPane.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(102, 102, 102)));
        headerScrollPane.setForeground(new java.awt.Color(153, 153, 153));
        headerScrollPane.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        headerScrollPane.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        headerScrollPane.setPreferredSize(new java.awt.Dimension(1021, 130));
        add(headerScrollPane, java.awt.BorderLayout.NORTH);

        applicationHeaderPanel = new IGVPanel(this);
        applicationHeaderPanel.add(nameHeaderPanel);
        applicationHeaderPanel.add(attributeHeaderPanel);
        applicationHeaderPanel.add(headerPanelContainer);
        headerScrollPane.setViewportView(applicationHeaderPanel);

        dataTrackScrollPane = new TrackPanelScrollPane();
        dataTrackScrollPane.setPreferredSize(new java.awt.Dimension(1021, 349));

        final TrackPanel dataTrackPanel = new TrackPanel(IGV.DATA_PANEL_NAME, this);
        dataTrackScrollPane.setViewportView(dataTrackPanel);

        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {
            featureTrackScrollPane = new TrackPanelScrollPane();
            featureTrackScrollPane.setPreferredSize(new java.awt.Dimension(1021, 50));
            featureTrackScrollPane.setViewportView(new TrackPanel(IGV.FEATURE_PANEL_NAME, this));
            // add(featureTrackScrollPane, java.awt.BorderLayout.SOUTH);
        }


        centerSplitPane = new SplitPane() {

            @Override
            public Insets getInsets(Insets insets) {
                return new Insets(0, 0, 0, 0);
            }
        };
        centerSplitPane.setDividerSize(3);
        //centerSplitPane.setResizeWeight(0.5d);
        centerSplitPane.setOrientation(JSplitPane.VERTICAL_SPLIT);

        centerSplitPane.add(dataTrackScrollPane, JSplitPane.TOP);
        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {
            centerSplitPane.add(featureTrackScrollPane, JSplitPane.BOTTOM);
        }

        add(centerSplitPane, BorderLayout.CENTER);

        setBackground(PreferencesManager.getPreferences().getAsColor(BACKGROUND_COLOR));


    }


    /**
     * Reset panels to initial state
     */
    public void resetPanels() {
        for (TrackPanel tp : getTrackPanels()) {
            tp.reset();
        }
        for (TrackPanel tp : getTrackPanels()) {
            final TrackPanelScrollPane tsp = tp.getScrollPane();
            if (tsp == dataTrackScrollPane || tsp == featureTrackScrollPane) {
                continue;
            }
            centerSplitPane.remove(tsp);
            TrackNamePanel.removeDropListenerFor(tsp.getNamePanel());
        }
    }

    /**
     * Add a new data panel set
     */
    public synchronized TrackPanelScrollPane addDataPanel(String name) {

        final TrackPanel trackPanel = new TrackPanel(name, this);
        final TrackPanelScrollPane sp = new TrackPanelScrollPane();
        Runnable runnable = () -> {

            sp.setViewportView(trackPanel);

            for (TrackPanel tp : getTrackPanels()) {
                tp.getScrollPane().minimizeHeight();
            }

            // Insert the new panel just before the feature panel, or at the end if there is no feature panel.
            int featurePaneIdx = centerSplitPane.indexOfPane(featureTrackScrollPane);
            if (featurePaneIdx > 0) {
                centerSplitPane.insertPane(sp, featurePaneIdx);
            } else {
                centerSplitPane.add(sp);
            }

            if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_SINGLE_TRACK_PANE_KEY)) {
                if (sp.getTrackPanel().getTracks().size() == 0) {
                    //If the igv window is too small the divider won't exist and this causes an exception
                    //We solved by setting a minimum size
                    centerSplitPane.setDividerLocation(0, 3);
                }
            }
        };

        UIUtilities.invokeAndWaitOnEventThread(runnable);

        return sp;
    }

    /**
     * Return an ordered list of TrackPanels.  This method is provided primarily for storing sessions, where
     * TrackPanels need to be stored in proper order
     *
     * @return
     */
    public java.util.List<TrackPanel> getTrackPanels() {
        ArrayList<TrackPanel> panels = new ArrayList<TrackPanel>();
        for (Component c : centerSplitPane.getComponents()) {
            if (c instanceof TrackPanelScrollPane) {
                panels.add(((TrackPanelScrollPane) c).getTrackPanel());
            }
        }

        return panels;
    }

    public void reorderPanels(java.util.List<String> names) {

        // First get visibile "heights" (distance between split pane dividers)
        int h = centerSplitPane.getHeight();
        int[] dividerLocations = centerSplitPane.getDividerLocations();
        Map<String, Integer> panelHeights = new HashMap();
        int idx = 0;

        Map<String, TrackPanelScrollPane> panes = new HashMap();
        for (Component c : centerSplitPane.getComponents()) {
            if (c instanceof TrackPanelScrollPane) {
                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;
                panes.put(tsp.getTrackPanelName(), tsp);
                int top = idx == 0 ? 0 : dividerLocations[idx - 1];
                int bottom = idx < dividerLocations.length ? dividerLocations[idx] : h;
                panelHeights.put(tsp.getTrackPanelName(), (bottom - top));
                idx++;
            }
        }

        //
        centerSplitPane.removeAll();
        idx = 0;
        int divLoc = 0;
        for (String name : names) {
            centerSplitPane.add(panes.get(name));
            if (idx < dividerLocations.length) {
                divLoc += panelHeights.get(name);
                dividerLocations[idx] = divLoc;
                idx++;
            }
        }
        centerSplitPane.setDividerLocations(dividerLocations);
        centerSplitPane.invalidate();
    }


    public void tweakPanelDivider() {
        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                // TODO Resize the data panel to make as much space as possible
                int h = centerSplitPane.getHeight();
                int nPanes = centerSplitPane.getPaneCount();

                double prefHeight = 0;
                for (int i = 0; i < nPanes; i++) {
                    prefHeight += centerSplitPane.getPaneAt(i).getPreferredSize().getHeight();
                }
                double ratio = h / prefHeight;
                int pos = 0;
                for (int i = 0; i < nPanes - 1; i++) {

                    pos += (int) (ratio * centerSplitPane.getPaneAt(i).getPreferredSize().getHeight());
                    centerSplitPane.setDividerLocation(i, pos);
                }
            }
        });

    }

    public void removeEmptyDataPanels() {
        List<TrackPanelScrollPane> emptyPanels = new ArrayList();
        for (TrackPanel tp : getTrackPanels()) {
            if (tp.getTracks().isEmpty()) {
                emptyPanels.add(tp.getScrollPane());
            }
        }
        for (TrackPanelScrollPane panel : emptyPanels) {
            if (panel != null) {
                centerSplitPane.remove(panel);
                TrackNamePanel.removeDropListenerFor(panel.getNamePanel());
            }

        }

    }

    public void removeDataPanel(String name) {

        for (TrackPanel tp : getTrackPanels()) {
            if (name.equals(tp.getName())) {
                removeTrackPanel(tp);
                return;
            }
        }
    }

    public void removeTrackPanel(TrackPanel trackPanel) {
        // Don't remove the "special" panes
        if (panelIsRemovable(trackPanel)) {
            TrackPanelScrollPane sp = trackPanel.getScrollPane();
            if (sp != null) {
                centerSplitPane.remove(sp);
                TrackNamePanel.removeDropListenerFor(sp.getNamePanel());
                centerSplitPane.revalidate();
            }
        }
    }

    public boolean panelIsRemovable(TrackPanel trackPanel) {
        return trackPanel.getScrollPane() != dataTrackScrollPane && trackPanel.getScrollPane() != featureTrackScrollPane;
    }

    public void updatePanelDimensions() {
        Insets insets = applicationHeaderPanel.getInsets();
        namePanelX = insets.left;
        attributePanelX = namePanelX + namePanelWidth + hgap;
        attributePanelWidth = calculateAttributeWidth();
        dataPanelX = attributePanelX + attributePanelWidth + hgap;
        dataPanelWidth = applicationHeaderPanel.getWidth() - insets.right - dataPanelX;
    }

    public int calculateAttributeWidth() {

        if (!PreferencesManager.getPreferences().getAsBoolean(SHOW_ATTRIBUTE_VIEWS_KEY)) {
            return 0;
        }
        Collection<String> attributeKeys = AttributeManager.getInstance().getVisibleAttributes();
        int attributeCount = attributeKeys.size();
        int packWidth = (attributeCount) * (AttributeHeaderPanel.ATTRIBUTE_COLUMN_WIDTH +
                AttributeHeaderPanel.COLUMN_BORDER_WIDTH) + AttributeHeaderPanel.COLUMN_BORDER_WIDTH;
        return packWidth;
    }

    public boolean isExpanded() {
        return namePanelWidth > 0;
    }

    public int getAttributePanelWidth() {
        return attributePanelWidth;
    }

    public int getNamePanelX() {
        return namePanelX;
    }

    public int getNamePanelWidth() {
        return namePanelWidth;
    }

    public int getAttributePanelX() {
        return attributePanelX;
    }


    public int getDataPanelX() {
        return dataPanelX;
    }


    public int getDataPanelWidth() {
        return dataPanelWidth;
    }

    public JideSplitPane getCenterSplitPane() {
        return centerSplitPane;
    }

    static class SplitPane extends JideSplitPane {
//        @Override
//        public void doLayout() {
//            if (log.isTraceEnabled()) {
//                log.trace("Layout");
//            }
//            super.doLayout();    //To change body of overridden methods use File | Settings | File Templates.
//        }
    }


    public void paintOffscreen(Graphics2D g, Rectangle rect, boolean batch) {

        // Header
        int width = applicationHeaderPanel.getWidth();
        int height = applicationHeaderPanel.getHeight();

        Graphics2D headerGraphics = (Graphics2D) g.create();
        Rectangle headerRect = new Rectangle(0, 0, width, height);
        applicationHeaderPanel.paintOffscreen(headerGraphics, headerRect, batch);
        headerGraphics.dispose();

        // Now loop through track panels
        Rectangle r = centerSplitPane.getBounds();
        g.translate(0, r.y);

        // Get the components of the center pane and sort by Y position.
        Component[] components = centerSplitPane.getComponents();
        Arrays.sort(components, Comparator.comparingInt(Component::getY));

        int dy = components[0].getY();
        for (Component c : components) {

            Graphics2D g2d = (Graphics2D) g.create();
            g2d.translate(0, dy);

            if (c instanceof TrackPanelScrollPane) {
                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;

                //Skip if panel has no tracks
                if (tsp.getTrackPanel().getTracks().size() == 0) {
                    continue;
                }

                int panelHeight = tsp.getSnapshotHeight(batch);

                Rectangle tspRect = new Rectangle(0, 0, tsp.getWidth(), panelHeight);

                g2d.setClip(tspRect);
                tsp.paintOffscreen(g2d, tspRect, batch);
                dy += tspRect.height;

            } else {
                g2d.setClip(new Rectangle(0, 0, c.getWidth(), c.getHeight()));
                c.paint(g2d);
                dy += c.getHeight();
            }

            g2d.dispose();

        }

//        //super.paintBorder(g);

    }

    /**
     * Return the image height required to paint this component with current options.  This is used to size bitmap
     * images for offscreen drawing.
     *
     * @return
     */
    @Override
    public int getSnapshotHeight(boolean batch) {

        if (batch) {
            int height = applicationHeaderPanel.getHeight();

            for (Component c : centerSplitPane.getComponents()) {

                if (c instanceof TrackPanelScrollPane) {

                    TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;

                    //Skip if panel has no tracks
                    if (tsp.getTrackPanel().getTracks().size() == 0) {
                        continue;
                    }

                    height += tsp.getSnapshotHeight(batch);

                } else {
                    height += c.getHeight();
                }

            }
            return height;
        } else {
            return getHeight();
        }
    }


}
