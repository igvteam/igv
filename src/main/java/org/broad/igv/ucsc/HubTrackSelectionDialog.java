package org.broad.igv.ucsc;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.util.HyperlinkFactory;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Dialog to enable selection of tracks defined by track hubs.  Modifies the "visible" property of
 * supplied track configurations in place.
 */
public class HubTrackSelectionDialog extends JDialog {

    private static final Dimension TRACK_SPACER = new Dimension(10, 5);
    private Map<JCheckBox, TrackConfig> configMap;

    public HubTrackSelectionDialog(List<TrackConfigGroup> groupedTrackConfigurations, Frame owner) {
        super(owner);
        setModal(true);
        init(groupedTrackConfigurations);
        setLocationRelativeTo(owner);
    }


    void init(List<TrackConfigGroup> trackConfigurations) {

        configMap = new HashMap<>();

        setTitle("Select tracks to load");
//        JLabel headerMessage = new JLabel("Select tracks to load");
//        headerMessage.setFont(FontManager.getFont(Font.BOLD, 14));
//        headerMessage.setHorizontalAlignment(SwingConstants.CENTER);
//        headerMessage.setPreferredSize(new Dimension(300, 50));
//        add(headerMessage, BorderLayout.NORTH);

        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BorderLayout());
        JScrollPane scrollPane = new JScrollPane(mainPanel);
        add(scrollPane, BorderLayout.CENTER);

        JPanel categoryContainer = new JPanel();
        categoryContainer.setLayout(new BoxLayout(categoryContainer, BoxLayout.PAGE_AXIS));
        mainPanel.add(categoryContainer);

        // Panel for select all/none
        JPanel checkAllPanel = new JPanel();
        ((FlowLayout) checkAllPanel.getLayout()).setAlignment(FlowLayout.LEFT);
        JButton selectAllButton = new JButton("Select All");
        selectAllButton.setFocusPainted(false);
        selectAllButton.addActionListener(e -> configMap.keySet().forEach(cb -> cb.setSelected(true)));
        checkAllPanel.add(selectAllButton);
        JButton selectNoneButton = new JButton("Select None");
        selectNoneButton.addActionListener(e -> configMap.keySet().forEach(cb -> cb.setSelected(false)));
        checkAllPanel.add(selectNoneButton);
        categoryContainer.add(checkAllPanel);
        //mainPanel.add(checkAllPanel, BorderLayout.NORTH);


        List<JPanel> cpl = new ArrayList<>();
        for (TrackConfigGroup configGroup : trackConfigurations) {
            categoryContainer.add(Box.createVerticalStrut(10));
            JPanel categoryPanel = categoryPanel(configGroup);
            categoryContainer.add(categoryPanel);
            cpl.add(categoryPanel);
        }

        JPanel buttonPanel = new JPanel();
        ((FlowLayout) buttonPanel.getLayout()).setAlignment(FlowLayout.RIGHT);

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> setVisible(false));

        JButton okButton = new JButton("OK");
        okButton.addActionListener(e -> okAction());

        if (Globals.IS_MAC) {
            buttonPanel.add(cancelButton);
            buttonPanel.add(okButton);
        } else {
            buttonPanel.add(okButton);
            buttonPanel.add(cancelButton);
        }

        getRootPane().setDefaultButton(okButton);

        mainPanel.add(buttonPanel, BorderLayout.SOUTH);

        // Try to adjust height to "just enough", BoxLayout will potentially leave empty space otherwise.  There might
        // be a better way to achieve this.
        int h = 75;
        for (JPanel p : cpl) h += p.getMinimumSize().height;
        setSize(new Dimension(800, Math.min(800, h)));
        pack();
        revalidate();

    }

    private void okAction() {
        for (Map.Entry<JCheckBox, TrackConfig> entry : configMap.entrySet()) {
            final TrackConfig trackConfig = entry.getValue();
            if (entry.getKey().isSelected()) {
                trackConfig.visible = true;
            } else {
                trackConfig.visible = false;
            }
        }
        setVisible(false);
    }

    /**
     * Create a JPanel for a particular track category, including a label and checkboxes for contained tracks.
     *
     * @param configGroup
     * @return
     */
    private JPanel categoryPanel(TrackConfigGroup configGroup) {

        JPanel container = new JPanel();
        container.setLayout(new BorderLayout());
        Border border = BorderFactory.createLineBorder(Color.lightGray);//   BorderFactory.createLoweredBevelBorder();
        container.setBorder(BorderFactory.createTitledBorder(border, configGroup.label));

        JPanel trackContainer = new JPanel();
        final WrapLayout wrapLayout = new WrapLayout();
        wrapLayout.setAlignment(FlowLayout.LEFT);
        trackContainer.setLayout(wrapLayout);

        container.add(trackContainer);
        for (TrackConfig trackConfig : configGroup.tracks) {

            JPanel p = new JPanel();
            final JCheckBox checkBox = new JCheckBox();
            configMap.put(checkBox, trackConfig);
            checkBox.setSelected(trackConfig.visible);

            JLabel label = trackConfig.html == null ?
                    new JLabel(trackConfig.name) :
                    HyperlinkFactory.createLink(trackConfig.name, trackConfig.html);
            label.setLabelFor(checkBox);

            p.add(checkBox);
            p.add(label);
            p.add(Box.createRigidArea(TRACK_SPACER));
            trackContainer.add(p);
        }

        return container;
    }

    /**
     * Convenience method to extract and return selected track configurations.
     *
     * @return
     */
    public List<TrackConfig> getSelectedConfigs() {
        List<TrackConfig> selected = configMap.values().stream().filter(trackConfig -> trackConfig.visible).collect(Collectors.toList());
        selected.sort((o1, o2) -> o1.order - o2.order);
        return selected;
    }

    /**
     * A FlowLayout extension that wraps components as needed and updates the preferre size accordingly.
     * FlowLayout itself does not update the preferred size, so the component width grows without bounds.
     * <p>
     * Code courtesy Rob Camick.  See https://tips4java.wordpress.com/2008/11/06/wrap-layout/
     * <p>
     * MIT License
     * <p>
     * Copyright (c) 2023 Rob Camick
     * <p>
     * Permission is hereby granted, free of charge, to any person obtaining a copy
     * of this software and associated documentation files (the "Software"), to deal
     * in the Software without restriction, including without limitation the rights
     * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
     * copies of the Software, and to permit persons to whom the Software is
     * furnished to do so, subject to the following conditions:
     * <p>
     * The above copyright notice and this permission notice shall be included in all
     * copies or substantial portions of the Software.
     * <p>
     * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
     * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
     * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
     * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
     * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
     * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
     * SOFTWARE.
     */
    static class WrapLayout extends FlowLayout {
        private Dimension preferredLayoutSize;

        /**
         * Constructs a new <code>WrapLayout</code> with a left
         * alignment and a default 5-unit horizontal and vertical gap.
         */
        public WrapLayout() {
            super();
        }

        /**
         * Constructs a new <code>FlowLayout</code> with the specified
         * alignment and a default 5-unit horizontal and vertical gap.
         * The value of the alignment argument must be one of
         * <code>WrapLayout</code>, <code>WrapLayout</code>,
         * or <code>WrapLayout</code>.
         *
         * @param align the alignment value
         */
        public WrapLayout(int align) {
            super(align);
        }

        /**
         * Creates a new flow layout manager with the indicated alignment
         * and the indicated horizontal and vertical gaps.
         * <p>
         * The value of the alignment argument must be one of
         * <code>WrapLayout</code>, <code>WrapLayout</code>,
         * or <code>WrapLayout</code>.
         *
         * @param align the alignment value
         * @param hgap  the horizontal gap between components
         * @param vgap  the vertical gap between components
         */
        public WrapLayout(int align, int hgap, int vgap) {
            super(align, hgap, vgap);
        }

        /**
         * Returns the preferred dimensions for this layout given the
         * <i>visible</i> components in the specified target container.
         *
         * @param target the component which needs to be laid out
         * @return the preferred dimensions to lay out the
         * subcomponents of the specified container
         */
        @Override
        public Dimension preferredLayoutSize(Container target) {
            return layoutSize(target, true);
        }

        /**
         * Returns the minimum dimensions needed to layout the <i>visible</i>
         * components contained in the specified target container.
         *
         * @param target the component which needs to be laid out
         * @return the minimum dimensions to lay out the
         * subcomponents of the specified container
         */
        @Override
        public Dimension minimumLayoutSize(Container target) {
            Dimension minimum = layoutSize(target, false);
            minimum.width -= (getHgap() + 1);
            return minimum;
        }

        /**
         * Returns the minimum or preferred dimension needed to layout the target
         * container.
         *
         * @param target    target to get layout size for
         * @param preferred should preferred size be calculated
         * @return the dimension to layout the target container
         */
        private Dimension layoutSize(Container target, boolean preferred) {
            synchronized (target.getTreeLock()) {
                //  Each row must fit with the width allocated to the containter.
                //  When the container width = 0, the preferred width of the container
                //  has not yet been calculated so lets ask for the maximum.

                int targetWidth = target.getSize().width;
                Container container = target;

                while (container.getSize().width == 0 && container.getParent() != null) {
                    container = container.getParent();
                }

                targetWidth = container.getSize().width;

                if (targetWidth == 0)
                    targetWidth = Integer.MAX_VALUE;

                int hgap = getHgap();
                int vgap = getVgap();
                Insets insets = target.getInsets();
                int horizontalInsetsAndGap = insets.left + insets.right + (hgap * 2);
                int maxWidth = targetWidth - horizontalInsetsAndGap;

                //  Fit components into the allowed width

                Dimension dim = new Dimension(0, 0);
                int rowWidth = 0;
                int rowHeight = 0;

                int nmembers = target.getComponentCount();

                for (int i = 0; i < nmembers; i++) {
                    Component m = target.getComponent(i);

                    if (m.isVisible()) {
                        Dimension d = preferred ? m.getPreferredSize() : m.getMinimumSize();

                        //  Can't add the component to current row. Start a new row.

                        if (rowWidth + d.width > maxWidth) {
                            addRow(dim, rowWidth, rowHeight);
                            rowWidth = 0;
                            rowHeight = 0;
                        }

                        //  Add a horizontal gap for all components after the first

                        if (rowWidth != 0) {
                            rowWidth += hgap;
                        }

                        rowWidth += d.width;
                        rowHeight = Math.max(rowHeight, d.height);
                    }
                }

                addRow(dim, rowWidth, rowHeight);

                dim.width += horizontalInsetsAndGap;
                dim.height += insets.top + insets.bottom + vgap * 2;

                //	When using a scroll pane or the DecoratedLookAndFeel we need to
                //  make sure the preferred size is less than the size of the
                //  target containter so shrinking the container size works
                //  correctly. Removing the horizontal gap is an easy way to do this.

                Container scrollPane = SwingUtilities.getAncestorOfClass(JScrollPane.class, target);

                if (scrollPane != null && target.isValid()) {
                    dim.width -= (hgap + 1);
                }

                return dim;
            }
        }

        /*
         *  A new row has been completed. Use the dimensions of this row
         *  to update the preferred size for the container.
         *
         *  @param dim update the width and height when appropriate
         *  @param rowWidth the width of the row to add
         *  @param rowHeight the height of the row to add
         */
        private void addRow(Dimension dim, int rowWidth, int rowHeight) {
            dim.width = Math.max(dim.width, rowWidth);

            if (dim.height > 0) {
                dim.height += getVgap();
            }

            dim.height += rowHeight;
        }
    }


    /**
     * main for testing and development
     *
     * @param args
     * @throws InterruptedException
     * @throws InvocationTargetException
     * @throws IOException
     */
    public static void main(String[] args) throws InterruptedException, InvocationTargetException, IOException {

        String hubFile = "test/data/hubs/hub.txt";
        Hub hub = Hub.loadHub(hubFile);
        List<TrackConfigGroup> groupedTrackConfigurations = hub.getGroupedTrackConfigurations();

        final HubTrackSelectionDialog dlf = new HubTrackSelectionDialog(groupedTrackConfigurations, null);

        dlf.setVisible(true);

        for (TrackConfig config : dlf.getSelectedConfigs()) {
            System.out.println(config.name);
        }
    }

}
