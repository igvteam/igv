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

import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
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
    private TrackPanelScrollPane featureTrackScrollPane;
    private JPanel outerScrollContainer;
    private NameHeaderPanel nameHeaderPanel;
    private AttributeHeaderPanel attributeHeaderPanel;

    private int hgap = 5;
    private JScrollPane headerScrollPane;


    public MainPanel(IGV igv) {

        this.igv = igv;

        initComponents();

        addComponentListener(new ComponentAdapter() {
            @Override
            public void componentResized(ComponentEvent componentEvent) {
                revalidateTrackPanels();
                igv.repaint();
            }
        });
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

        applicationHeaderPanel = new IGVPanel(this);
        applicationHeaderPanel.add(nameHeaderPanel);
        applicationHeaderPanel.add(attributeHeaderPanel);
        applicationHeaderPanel.add(headerPanelContainer);
        applicationHeaderPanel.setPreferredSize(new java.awt.Dimension(1021, 130));

        // The header panel is wrapped in a scroll pane, never used, so that it will align with the data panels.
        JScrollPane headerPanelScrollpane = new JScrollPane(applicationHeaderPanel);
        headerPanelScrollpane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        headerPanelScrollpane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);

        outerScrollContainer = new JPanel();
        outerScrollContainer.setLayout(new MainPanelLayout(outerScrollContainer));

        JScrollPane outerScrollPane = new JScrollPane(outerScrollContainer);
        outerScrollPane.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        outerScrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        outerScrollPane.setColumnHeaderView(headerPanelScrollpane);

        add(outerScrollPane, BorderLayout.CENTER);

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
            if (tsp == featureTrackScrollPane) {
                continue;
            }
            outerScrollContainer.remove(tsp);
            TrackNamePanel.removeDropListenerFor(tsp.getNamePanel());
        }
    }

    public void addTrack(Track t) {

        Runnable runnable = () -> {
            final TrackPanel trackPanel = new TrackPanel(t.getName(), this);
            trackPanel.addTrack(t);
            t.setTrackPanel(trackPanel);

            final TrackPanelScrollPane sp = new TrackPanelScrollPane();
            sp.setViewportView(trackPanel);
            trackPanel.setScrollPane(sp);

            Insets insets2 = sp.getInsets();
            int dy = insets2.top + insets2.bottom;

            // trackPanel.setPreferredSize(new Dimension(1000, t.getHeight() + dy));
            // trackPanel.setMinimumSize(new Dimension(0, t.getMinimumHeight() + dy));
            final int height = t.getHeight() + dy;
            sp.setPreferredSize(new Dimension(1000, height));
            sp.setMinimumSize(new Dimension(0, t.getMinimumHeight() + dy));
            outerScrollContainer.add(sp);

            trackPanel.setVisible(t.isVisible());

        };
        UIUtilities.invokeAndWaitOnEventThread(runnable);
    }

    /**
     * Add a new data panel set
     */
    public synchronized TrackPanelScrollPane addDataPanel(String name) {

        final TrackPanel trackPanel = new TrackPanel(name, this);
        final TrackPanelScrollPane sp = new TrackPanelScrollPane();
        Runnable runnable = () -> {
            sp.setViewportView(trackPanel);
            trackPanel.setPreferredSize(new Dimension(1000, 200));
            sp.setPreferredSize(new Dimension(1000, 200));
            outerScrollContainer.add(sp);
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
        for (Component c : outerScrollContainer.getComponents()) {
            if (c instanceof TrackPanelScrollPane) {
                panels.add(((TrackPanelScrollPane) c).getTrackPanel());
            }
        }

        return panels;
    }

    public void reorderPanels(java.util.List<String> names) {

        Map<String, TrackPanelScrollPane> panes = new LinkedHashMap<>();
        for (Component c : outerScrollContainer.getComponents()) {
            if (c instanceof TrackPanelScrollPane) {
                TrackPanelScrollPane tsp = (TrackPanelScrollPane) c;
                panes.put(tsp.getTrackPanelName(), tsp);
            }
        }

        outerScrollContainer.removeAll();
        for (String name : names) {
            outerScrollContainer.add(panes.get(name));
        }

        outerScrollContainer.revalidate();
    }


    public void removeEmptyPanels() {
        List<TrackPanelScrollPane> emptyPanels = new ArrayList();
        for (TrackPanel tp : getTrackPanels()) {
            if (tp.getTracks().isEmpty()) {
                emptyPanels.add(tp.getScrollPane());
            }
        }
        for (TrackPanelScrollPane panel : emptyPanels) {
            if (panel != null) {
                outerScrollContainer.remove(panel);
                TrackNamePanel.removeDropListenerFor(panel.getNamePanel());
            }

        }

    }

    public void removeDataPanel(Component panel) {

        outerScrollContainer.remove(panel);
        // TrackNamePanel.removeDropListenerFor(sp.getNamePanel());

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

    public JPanel getOuterScrollContainer() {
        return outerScrollContainer;
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
        Rectangle r = outerScrollContainer.getBounds();
        g.translate(0, r.y);

        // Get the components of the center pane and sort by Y position.
        Component[] components = outerScrollContainer.getComponents();
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

            for (Component c : outerScrollContainer.getComponents()) {

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
