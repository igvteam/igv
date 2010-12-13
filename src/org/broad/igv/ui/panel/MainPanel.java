/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.ui.panel;

import com.jidesoft.swing.JideSplitPane;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackManager;
import org.broad.igv.ui.util.SnapshotUtilities;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;

/**
 * @author jrobinso
 * @date Sep 10, 2010
 */
public class MainPanel extends JPanel implements Paintable {

    private static Logger log = Logger.getLogger(MainPanel.class);

    TrackManager trackManager;

    private static final int DEFAULT_NAME_PANEL_WIDTH = 160;

    private int namePanelX;
    private int namePanelWidth = DEFAULT_NAME_PANEL_WIDTH;
    private int attributePanelX;
    private int attributePanelWidth;
    private int dataPanelX;
    private int dataPanelWidth;


    public IGVPanel applicationHeaderPanel;
    public HeaderPanelContainer headerPanelContainer;
    private TrackPanelScrollPane dataTrackScrollPane;
    private TrackPanelScrollPane featureTrackScrollPane;
    private JideSplitPane centerSplitPane;

    private int hgap = 5;


    public MainPanel(TrackManager trackManager) {
        this.trackManager = trackManager;
        initComponents();
    }


    public void collapseNamePanel() {
        namePanelWidth = 0;
        revalidate();
    }

    public void expandNamePanel() {
        namePanelWidth = DEFAULT_NAME_PANEL_WIDTH;
        revalidate();
    }


    @Override
    public void doLayout() {
        layoutFrames();
        super.doLayout();    //To change body of overridden methods use File | Settings | File Templates.
        applicationHeaderPanel.doLayout();
        for (TrackPanelScrollPane tsp : trackManager.getTrackPanelScrollPanes()) {
            tsp.doLayout();
        }
    }

    private void initComponents() {

        setBackground(new java.awt.Color(204, 204, 204));
        setPreferredSize(new java.awt.Dimension(1021, 510));
        setLayout(new java.awt.BorderLayout());


        NameHeaderPanel nameHeaderPanel = new NameHeaderPanel();
        nameHeaderPanel.setBackground(new java.awt.Color(255, 255, 255));
        nameHeaderPanel.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        nameHeaderPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        nameHeaderPanel.setPreferredSize(new java.awt.Dimension(0, 0));
        nameHeaderPanel.setLayout(null);

        AttributeHeaderPanel attributeHeaderPanel = new AttributeHeaderPanel();
        attributeHeaderPanel.setBackground(new java.awt.Color(255, 255, 255));
        attributeHeaderPanel.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        attributeHeaderPanel.setDebugGraphicsOptions(javax.swing.DebugGraphics.NONE_OPTION);
        attributeHeaderPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        attributeHeaderPanel.setPreferredSize(new java.awt.Dimension(0, 0));


        headerPanelContainer = new HeaderPanelContainer();


        JScrollPane headerScrollPane = new JScrollPane();
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

        final TrackPanel dataTrackPanel = new TrackPanel(TrackManager.DATA_PANEL_NAME, this);
        dataTrackScrollPane.setViewportView(dataTrackPanel);
        trackManager.putScrollPane(TrackManager.DATA_PANEL_NAME, dataTrackScrollPane);

        if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY)) {
            featureTrackScrollPane = new TrackPanelScrollPane();
            featureTrackScrollPane.setPreferredSize(new java.awt.Dimension(1021, 50));
            featureTrackScrollPane.setViewportView(new TrackPanel(TrackManager.FEATURE_PANEL_NAME, this));
            add(featureTrackScrollPane, java.awt.BorderLayout.SOUTH);
            trackManager.putScrollPane(TrackManager.FEATURE_PANEL_NAME, featureTrackScrollPane);
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
        if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY)) {
            centerSplitPane.add(featureTrackScrollPane, JSplitPane.BOTTOM);
        }


        /*JPanel testPanel = new JPanel() {
            @Override
        public void paintComponent(Graphics g) {
               g.setColor(Color.blue);
               ((Graphics2D) g).fill(getBounds());
            }
        };*/

        add(centerSplitPane, BorderLayout.CENTER);


    }

    public void adjustPanelDivider() {
        // Adjust divider for data panel.  The data panel divider can be
        // zero if there are no data tracks loaded.
        TrackPanelScrollPane dsp = trackManager.getScrollPane(TrackManager.DATA_PANEL_NAME);
        if (dsp.getDataPanel().getAllTracks().size() > 0 &&
                centerSplitPane.getDividerLocation(0) < 10) {
            centerSplitPane.setDividerLocation(0, 40);
        }
    }

    public void resetPanels() {
        // Remove user added panels
        for (TrackPanelScrollPane tsp : trackManager.getTrackPanelScrollPanes()) {
            tsp.getTrackPanel().clearTracks();
            if (tsp == dataTrackScrollPane || tsp == featureTrackScrollPane) {
                continue;
            }
            centerSplitPane.remove(tsp);
            TrackNamePanel.removeDropListenerFor(tsp.getNamePanel());
        }

        trackManager.reset();
        trackManager.clearScrollPanes();
        trackManager.putScrollPane(TrackManager.DATA_PANEL_NAME, dataTrackScrollPane);
        Track geneTrack = trackManager.getGeneTrack();
        if (PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY)) {
            if (geneTrack != null) {
                dataTrackScrollPane.getTrackPanel().addTrack(geneTrack);
            }
        } else {
            trackManager.putScrollPane(TrackManager.FEATURE_PANEL_NAME, featureTrackScrollPane);
            if (geneTrack != null) {
                featureTrackScrollPane.getTrackPanel().addTrack(geneTrack);
            }
        }
    }

    /**
     * Add a new data panel set
     */
    public TrackPanelScrollPane addDataPanel(String name) {

        TrackPanel trackPanel = new TrackPanel(name, this);
        final TrackPanelScrollPane sp = new TrackPanelScrollPane();
        sp.setViewportView(trackPanel);
        //sp.setPreferredSize(new Dimension(700, 300));

        for (TrackPanelScrollPane tsp : trackManager.getTrackPanelScrollPanes()) {
            tsp.minimizeHeight();
        }

        trackManager.putScrollPane(name, sp);

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                // TODO Resize the data panel to make as much space as possible

                // Insert the new panel just before the feature panel, or at the end if there is no feature panel.
                int featurePaneIdx = centerSplitPane.indexOfPane(featureTrackScrollPane);
                if (featurePaneIdx > 0) {
                    centerSplitPane.insertPane(sp, featurePaneIdx);
                } else {
                    centerSplitPane.add(sp);
                }

                if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY)) {
                    if (sp.getTrackPanel().getTracks().size() == 0) {
                        centerSplitPane.setDividerLocation(0, 3);
                    }
                }

            }
        });

        return sp;
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

    public void removeDataPanel(String name) {
        TrackPanelScrollPane sp = trackManager.getScrollPane(name);
        // Don't remove the "special" panes
        if (sp == dataTrackScrollPane || sp == featureTrackScrollPane) {
            return;
        }
        if (sp != null) {
            centerSplitPane.remove(sp);
            trackManager.removeScrollPane(name);
            TrackNamePanel.removeDropListenerFor(sp.getNamePanel());
        }
    }

    public void layoutFrames() {
        synchronized (getTreeLock()) {

            Insets insets = applicationHeaderPanel.getInsets();
            namePanelX = insets.left;

            attributePanelX = namePanelX + namePanelWidth + hgap;
            attributePanelWidth = calculateAttributeWidth();

            dataPanelX = attributePanelX + attributePanelWidth + hgap;

            java.util.List<ReferenceFrame> frames = FrameManager.getFrames();
            dataPanelWidth = applicationHeaderPanel.getWidth() - insets.right - dataPanelX;

            if (frames.size() == 1) {
                frames.get(0).setBounds(0, dataPanelWidth);
            } else {

                float gap = Math.min(1, 20.0f / ((int) (1.5 * frames.size()))) * hgap;
                int x = 0;

                // Width is in floating point because we need to fill data panel,  going straight to an "int" here
                // would cause truncation
                float wc = ((float) dataPanelWidth - (frames.size() - 1) * gap) / frames.size();
                for (int i = 0; i < frames.size(); i++) {
                    ReferenceFrame frame = frames.get(i);
                    int nextX = (int) ((i + 1) * (wc + gap));
                    int w = nextX - x;
                    frame.setBounds(x, w);
                    x = nextX;
                }
            }
        }
    }


    private int calculateAttributeWidth() {

        if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_ATTRIBUTE_VIEWS_KEY)) {
            return 0;
        }

        HashSet<String> attributeKeys = new HashSet(AttributeManager.getInstance().getAttributeKeys());
        attributeKeys.removeAll(AttributeManager.getInstance().getHiddenAttributes());

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


    static class SplitPane extends JideSplitPane {
        @Override
        public void doLayout() {
            if (log.isDebugEnabled()) {
                log.debug("Layout");
            }
            super.doLayout();    //To change body of overridden methods use File | Settings | File Templates.
        }
    }


    public void paintOffscreen(Graphics2D g, Rectangle rect) {

        g.setColor(Color.lightGray);
        g.fill(rect);


        // Header
        int width = applicationHeaderPanel.getWidth();
        int height = applicationHeaderPanel.getHeight();

        Rectangle headerRect = new Rectangle(0, 0, width, height);
        applicationHeaderPanel.paintOffscreen(g, headerRect);

        // Now loop through track panel
        Rectangle r = centerSplitPane.getBounds();


        g.translate(0, r.y);

        // Get the components of the center pane and sort by Y position.
        Component[] components = centerSplitPane.getComponents();
        Arrays.sort(components, new Comparator<Component>() {
            public int compare(Component component, Component component1) {
                return component.getY() - component1.getY();
            }
        });

        int dy = components[0].getY();
        for (Component c : components) {

            Graphics2D g2d = (Graphics2D) g.create();
            g2d.translate(0, dy);

            if (c instanceof TrackPanelScrollPane) {

                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;

                //Skip if panel has no tracks
                if(tsp.getTrackPanel().getTracks().size() == 0) {
                    continue;
                }

                int maxPanelHeight = SnapshotUtilities.MAX_PANEL_HEIGHT;
                int panelHeight = Math.min(maxPanelHeight, Math.max(tsp.getVisibleRect().height, tsp.getDataPanel().getHeight()));

                Rectangle tspRect = new Rectangle(tsp.getBounds());
                tspRect.height = panelHeight;

                g2d.setClip(new Rectangle(0, 0, tsp.getWidth(), tspRect.height));
                tsp.paintOffscreen(g2d, tspRect);

                dy += tspRect.height;

            } else {
                g2d.setClip(new Rectangle(0, 0, c.getWidth(), c.getHeight()));
                c.paint(g2d);
                dy += c.getHeight();
            }

            g2d.dispose();

        }

        super.paintBorder(g);

    }

    public int getOffscreenImageHeight() {
        int height = centerSplitPane.getBounds().y;
        for (Component c : centerSplitPane.getComponents()) {

            if (c instanceof TrackPanelScrollPane) {

                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;

                int maxPanelHeight = SnapshotUtilities.MAX_PANEL_HEIGHT;
                int panelHeight = Math.min(maxPanelHeight, Math.max(tsp.getVisibleRect().height, tsp.getDataPanel().getHeight()));

                Rectangle tspRect = new Rectangle(tsp.getBounds());
                tspRect.height = panelHeight;


                height += tspRect.height;
            } else {
                height += c.getHeight();
            }

        }
        // TODO Not sure why this is neccessary
        height += 35;
        return height;

    }


}
