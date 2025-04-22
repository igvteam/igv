package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.FontManager;
import org.broad.igv.ui.util.IconFactory;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.net.URI;
import java.util.*;
import java.util.List;
import java.util.function.Function;

public class CollapsiblePanel extends JPanel {

    private static Logger log = LogManager.getLogger(CollapsiblePanel.class);
    public static final Color HEADER_BG = new Color(180, 204, 226);
    public static final Color HEADER_BG2 = new Color(204, 204, 204);
    private final JLabel jlabel;
    final List<SelectionBox> selectionBoxes;
    private final boolean autoselectDefaults;
    private final TrackConfigContainer configContainer;
    private JButton collapseButton;
    private JComponent content;
    private JPanel header;
    private ImageIcon openIcon;
    private ImageIcon closeIcon;

    private static Set<String> autoselectTracks = Set.of("ncbiRefSeq", "augustus", "ensGene");

    public CollapsiblePanel(TrackConfigContainer configContainer, boolean autoselectDefaults, boolean addSearchButton) {

        Color backgroundColor = HEADER_BG;

        this.configContainer = configContainer;
        this.autoselectDefaults = autoselectDefaults;

        setLayout(new BorderLayout());

        JPanel trackPanel = new JPanel();
        trackPanel.setLayout(new BoxLayout(trackPanel, BoxLayout.Y_AXIS));

        // There is a (so far) intractable bug if a large # of JCheckboxes are created for this widget.
        int totalTrackCount = configContainer.countTracks();
        SelectionBox.CheckboxType checkboxType = totalTrackCount < 1000 ? SelectionBox.CheckboxType.SWING : SelectionBox.CheckboxType.CUSTOM;

        selectionBoxes = addSelectionBoxes(null, configContainer, trackPanel, checkboxType);

        boolean isSelected = false;
        int maxWidth = 0;
        int selectionCount = 0;
        for (SelectionBox selectionBox : selectionBoxes) {
            if (selectionBox.isSelected()) {
                selectionCount++;
                isSelected = true;
            }
            maxWidth = Math.max(maxWidth, selectionBox.getPreferredSize().width);
            //  allSelectionBoxes.add(selectionBox);
        }

        for (SelectionBox selectionBox : selectionBoxes) {
            selectionBox.setPreferredWidth(maxWidth);
        }

        String label = configContainer.label + "   (" + selectionBoxes.size() + " tracks, " + selectionCount + " selected)";

        boolean isOpen = isSelected || configContainer.defaultOpen;

        this.content = trackPanel;
        this.openIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.MINUS);
        this.closeIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.PLUS);

        content.setVisible(isOpen);
        this.add(content, BorderLayout.CENTER);

        header = new JPanel();
        header.setLayout(new BorderLayout());
        header.setBackground(backgroundColor);

        this.collapseButton = new JButton();
        this.collapseButton.setBackground(backgroundColor);
        this.collapseButton.setBorderPainted(false);
        this.collapseButton.setFocusPainted(false);
        this.collapseButton.setContentAreaFilled(false);
        collapseButton.setIcon(isOpen ? openIcon : closeIcon);
        //collapseButton.setBorder(new EmptyBorder(10, 5, 10, 0));
        collapseButton.setHorizontalAlignment(SwingConstants.LEFT);

        collapseButton.addActionListener(e -> {
            collapseButton.setIcon(content.isVisible() ? closeIcon : openIcon);
            content.setVisible(!content.isVisible());
            //this.getParent().revalidate();
        });
        header.add(collapseButton, BorderLayout.WEST);

        jlabel = new JLabel(label);
        jlabel.setFont(FontManager.getFont(14));
        jlabel.setHorizontalAlignment(SwingConstants.CENTER);
        header.add(jlabel, BorderLayout.CENTER);

        this.add(header, BorderLayout.NORTH);

        if (addSearchButton) {
            final JButton searchButton = TrackHubSelectionDialog.createSearchButton("Search " + configContainer.label, Arrays.asList(this),
                    (selectedCount) -> {
                        this.updateLabel();
                        return null;
                    });

            this.addSearchButton(searchButton);
        }

        for (SelectionBox selectionBox : selectionBoxes) {
            selectionBox.setCallback(b -> {
                this.updateLabel();
                return null;
            });
        }
    }

    public String containerLabel() {
        return configContainer.label;
    }

    public void resetSelectionBoxes(Set<String> loadedTrackPaths) {

        for (CollapsiblePanel.SelectionBox box : selectionBoxes) {
            final boolean isLoaded = loadedTrackPaths != null && loadedTrackPaths.contains(box.trackConfig.url);
            box.setSelected(
                    isLoaded ||
                            ((loadedTrackPaths == null || loadedTrackPaths.isEmpty()) && autoselectDefaults &&
                                    (autoselectTracks.contains(box.trackConfig.id) || selectionBoxes.size() == 1)));
            box.setEnabled(this.autoselectDefaults || !isLoaded);
        }
        updateLabel();
    }

    public void updateLabel() {
        int count = 0;
        for (SelectionBox selectionBox : selectionBoxes) {
            if (selectionBox.isEnabled() && selectionBox.isSelected()) {
                count++;
            }
        }
        String label = configContainer.label + "   (" + selectionBoxes.size() + " tracks, " + count + " selected)";
        jlabel.setText(label);
    }

    public void clearSelections() {
        for (CollapsiblePanel.SelectionBox box : selectionBoxes) {
            if (box.isEnabled()) {
                box.setSelected(false);
            }
        }
        String label = configContainer.label + "   (" + selectionBoxes.size() + " tracks, 0 selected)";
        jlabel.setText(label);
    }

    public List<TrackConfig> getSelectedTracks() {
        return selectionBoxes.stream()
                .filter(b -> b.isEnabled() && b.isSelected())
                .map(CollapsiblePanel.SelectionBox::getTrackConfig).toList();
    }


    public void addSearchButton(JComponent searchButton) {
        header.add(searchButton, BorderLayout.EAST);
        revalidate();
    }

    public void collapse() {
        collapseButton.setIcon(closeIcon);
        content.setVisible(false);
    }

    public void expand() {
        collapseButton.setIcon(openIcon);
        content.setVisible(true);
    }

    /**
     * Constrain the maximum height to prevent BoxLayout from needlessly resizing the panel to fill space.  This is
     * rather hardcoded for the TrackHubSelectionDialog.
     *
     * @return
     */

    @Override
    public Dimension getMaximumSize() {
        Dimension d4 = header.getMinimumSize();
        if (!content.isVisible()) {
            return new Dimension(Integer.MAX_VALUE, d4.height);
        } else {
            Dimension d5 = content.getMinimumSize();
            return new Dimension(Integer.MAX_VALUE, d4.height + (int) (1.2 * d5.height));
        }
    }

    /**
     * Add selection boxes for the container (a group, superTrack, compositeTrack, or view).  Return true if any
     * tracks are selected.
     *
     * @param labelPrefix
     * @param container
     * @param panel
     * @param checkboxType
     * @return
     */
    private java.util.List<SelectionBox> addSelectionBoxes(String labelPrefix,
                                                           TrackConfigContainer container,
                                                           JPanel panel,
                                                           SelectionBox.CheckboxType checkboxType) {

        String title = labelPrefix == null ? "" :
                labelPrefix + (labelPrefix.length() > 0 ? " - " : "") + container.label;

        java.util.List<SelectionBox> selectionBoxes = new ArrayList<>();

        if (container.tracks.size() > 0) {

            JPanel trackPanel = new JPanel();
            if (labelPrefix != null) {
                trackPanel.setBorder(BorderFactory.createTitledBorder(title));
            }
            final WrapLayout wrapLayout = new WrapLayout();
            wrapLayout.setAlignment(FlowLayout.LEFT);
            trackPanel.setLayout(wrapLayout);

            for (TrackConfig trackConfig : container.tracks) {

                SelectionBox selectionBox = new SelectionBox(trackConfig, checkboxType);
                trackPanel.add(selectionBox);
                selectionBoxes.add(selectionBox);
            }

            panel.add(Box.createVerticalStrut(5));
            panel.add(trackPanel);

        }

        for (TrackConfigContainer childChild : container.children) {
            selectionBoxes.addAll(addSelectionBoxes(title, childChild, panel, checkboxType));
        }
        return selectionBoxes;
    }


    static class SelectionBox extends JPanel {

        enum CheckboxType {SWING, CUSTOM}

        private final JLabel label;
        private TrackConfig trackConfig;
        private TrackHubSelectionDialog.CheckBoxWrapper checkbox;
        private int preferredWidth = -1;
        private int minWidth;
        private Function<Integer, Void> callback;

        public SelectionBox(TrackConfig trackConfig, CheckboxType checkboxType) {

            this.setLayout(new BorderLayout(5, 0));
            this.trackConfig = trackConfig;

            String longLabel = trackConfig.longLabel;
            if (longLabel != null) {
                this.setToolTipText(longLabel);
            }

            this.checkbox = new TrackHubSelectionDialog.CheckBoxWrapper(checkboxType);

            checkbox.setActionListener(e -> {
                trackConfig.visible = (checkbox.isSelected());
                if (callback != null) {
                    callback.apply(checkbox.isSelected() ? 1 : 0);
                }
            });

            label = new JLabel(trackConfig.name);
            label.setLabelFor(checkbox.getComponent());
            add(checkbox.getComponent(), BorderLayout.WEST);

            String infoLink = trackConfig.html;

            if (infoLink == null || "".equals(infoLink.trim())) {
                add(label, BorderLayout.CENTER);
            } else {
                ImageIcon icon = IconFactory.getInstance().getIcon(IconFactory.IconID.INFO);
                JLabel iconLabel = new JLabel(icon);
                iconLabel.setToolTipText(infoLink);
                iconLabel.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
                add(iconLabel, BorderLayout.CENTER);
                iconLabel.addMouseListener(new MouseAdapter() {
                    @Override
                    public void mouseClicked(MouseEvent e) {
                        try {
                            Desktop.getDesktop().browse(new URI(infoLink));
                        } catch (Exception ex) {
                            log.error("Error following hyperlink: " + infoLink, ex);
                        }
                    }
                });

                JPanel panel = new JPanel();
                panel.setLayout(new FlowLayout(FlowLayout.LEFT));
                panel.add(label);
                panel.add(iconLabel);
                add(panel, BorderLayout.CENTER);
            }
        }

        public void setPreferredWidth(int width) {
            minWidth = 200;
            this.preferredWidth = Math.max(minWidth, width);
        }

        @Override
        public Dimension getPreferredSize() {
            return preferredWidth > 0 ? new Dimension(preferredWidth, 20) : super.getPreferredSize();
        }

        @Override
        public Dimension getMinimumSize() {
            return getPreferredSize();
        }

        @Override
        public Dimension getMaximumSize() {
            return getPreferredSize();
        }

        public void setSelected(boolean selected) {
            checkbox.setSelected(selected);
            trackConfig.visible = (selected);
        }

        @Override
        public void setEnabled(boolean enabled) {
            super.setEnabled(enabled);
            checkbox.setEnabled(enabled);
            label.setEnabled(enabled);
        }


        public boolean isSelected() {
            return checkbox.isSelected();
        }

        public TrackConfig getTrackConfig() {
            return trackConfig;
        }

        public void setCallback(Function<Integer, Void> callback) {
            this.callback = callback;
        }
    }
}


