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

package org.broad.igv.ui.panel;

import com.jidesoft.swing.JideSplitPane;
import org.apache.log4j.Logger;
import org.broad.igv.PreferenceManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.SnapshotUtilities;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 * @date Sep 10, 2010
 */
public class MainPanel extends JPanel implements Paintable {

    private static Logger log = Logger.getLogger(MainPanel.class);

    IGV igv;

    // private static final int DEFAULT_NAME_PANEL_WIDTH = 160;

    private int namePanelX;
    private int namePanelWidth = PreferenceManager.getInstance().getAsInt(PreferenceManager.NAME_PANEL_WIDTH);
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
                revalidate();
                repaint();
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
        revalidate();
    }

    public void expandNamePanel() {
        namePanelWidth = PreferenceManager.getInstance().getAsInt(PreferenceManager.NAME_PANEL_WIDTH);
        revalidate();
    }

    public void setNamePanelWidth(int width) {
        this.namePanelWidth = width;
        revalidate();
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
        layoutFrames();
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
        nameHeaderPanel.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
        nameHeaderPanel.setMinimumSize(new java.awt.Dimension(0, 0));
        nameHeaderPanel.setPreferredSize(new java.awt.Dimension(0, 0));
        nameHeaderPanel.setLayout(new BorderLayout());

        attributeHeaderPanel = new AttributeHeaderPanel();
        attributeHeaderPanel.setBackground(new java.awt.Color(255, 255, 255));
        attributeHeaderPanel.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
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

        if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY)) {
            featureTrackScrollPane = new TrackPanelScrollPane();
            featureTrackScrollPane.setPreferredSize(new java.awt.Dimension(1021, 50));
            featureTrackScrollPane.setViewportView(new TrackPanel(IGV.FEATURE_PANEL_NAME, this));
            add(featureTrackScrollPane, java.awt.BorderLayout.SOUTH);
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

        add(centerSplitPane, BorderLayout.CENTER);

        setBackground(PreferenceManager.getInstance().getAsColor(PreferenceManager.BACKGROUND_COLOR));


    }


    /**
     * Called from IGV.resetSession()
     */
    public void resetPanels() {
        // Remove user added panels
        for (TrackPanel tp : getTrackPanels()) {
            tp.clearTracks();
            final TrackPanelScrollPane tsp = tp.getScrollPane();
            if (tsp == dataTrackScrollPane || tsp == featureTrackScrollPane) {
                continue;
            }
            centerSplitPane.remove(tsp);
            TrackNamePanel.removeDropListenerFor(tsp.getNamePanel());
        }

        igv.reset();
    }

    /**
     * Add a new data panel set
     */
    public synchronized TrackPanelScrollPane addDataPanel(String name) {

        TrackPanel trackPanel = new TrackPanel(name, this);
        final TrackPanelScrollPane sp = new TrackPanelScrollPane();
        sp.setViewportView(trackPanel);
        //sp.setPreferredSize(new Dimension(700, 300));

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

        if (!PreferenceManager.getInstance().getAsBoolean(PreferenceManager.SHOW_SINGLE_TRACK_PANE_KEY)) {
            if (sp.getTrackPanel().getTracks().size() == 0) {
                //If the igv window is too small the divider won't exist and this causes an exception
                //We solved by setting a minimum size
                centerSplitPane.setDividerLocation(0, 3);
            }
        }

        return sp;
    }

    /**
     * Return an ordered list of TrackPanels.  This method is provided primarily for storing sessions, where
     * TrackPanels need to be stored in proper order
     *
     * @return
     */
    public java.util.List<TrackPanel> getTrackPanels() {
        ArrayList panels = new ArrayList();
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

        TrackPanelScrollPane sp = null;
        for (TrackPanel tp : getTrackPanels()) {
            if (name.equals(tp.getName())) {
                sp = tp.getScrollPane();
                break;
            }
        }
        // Don't remove the "special" panes
        if (sp == dataTrackScrollPane || sp == featureTrackScrollPane) {
            return;
        }
        if (sp != null) {
            centerSplitPane.remove(sp);
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
        @Override
        public void doLayout() {
            if (log.isDebugEnabled()) {
                log.debug("Layout");
            }
            super.doLayout();    //To change body of overridden methods use File | Settings | File Templates.
        }
    }


    public void paintOffscreen(Graphics2D g, Rectangle rect) {

//        g.setColor(Color.lightGray);
//        g.fill(rect);


        // Header
        int width = applicationHeaderPanel.getWidth();
        int height = applicationHeaderPanel.getHeight();

        Graphics2D headerGraphics = (Graphics2D) g.create();
        headerGraphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        Rectangle headerRect = new Rectangle(0, 0, width, height);
        applicationHeaderPanel.paintOffscreen(headerGraphics, headerRect);
        headerGraphics.dispose();

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
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2d.translate(0, dy);

            if (c instanceof TrackPanelScrollPane) {

                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;

                //Skip if panel has no tracks
                if (tsp.getTrackPanel().getTracks().size() == 0) {
                    continue;
                }

                int panelHeight = getOffscreenImagePanelHeight(tsp);

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

        //super.paintBorder(g);

    }

    /**
     * Return the image height required to paint this component with current options.  This is used to size bitmap
     * images for offscreen drawing.
     *
     * @return
     */
    public int getOffscreenImageHeight() {
        int height = centerSplitPane.getBounds().y;
        for (Component c : centerSplitPane.getComponents()) {

            if (c instanceof TrackPanelScrollPane) {

                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;

                int panelHeight = getOffscreenImagePanelHeight(tsp);

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

    private int getOffscreenImagePanelHeight(TrackPanelScrollPane tsp) {
        int panelHeight;
        int maxPanelHeight = SnapshotUtilities.getMaxPanelHeight();
        final int visibleHeight = tsp.getVisibleRect().height;
        if (maxPanelHeight < 0) {
            panelHeight = visibleHeight;
        } else {
            panelHeight = Math.min(maxPanelHeight, Math.max(visibleHeight, tsp.getDataPanel().getHeight()));
        }
        return panelHeight;
    }


}
