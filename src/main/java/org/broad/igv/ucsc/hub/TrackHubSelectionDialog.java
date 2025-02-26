package org.broad.igv.ucsc.hub;

import org.broad.igv.encode.FileRecord;
import org.broad.igv.encode.TrackChooser;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.HyperlinkFactory;
import org.broad.igv.ui.util.IconFactory;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.net.URI;
import java.util.*;
import java.util.List;

import java.util.function.Function;


/**
 * Dialog to enable selection of tracks defined by track hubs.  Modifies the "visible" property of
 * supplied track configurations.
 */
public class TrackHubSelectionDialog extends JDialog {

    private static Logger log = LogManager.getLogger(TrackHubSelectionDialog.class);

    private static Map<Hub, TrackHubSelectionDialog> hubSelectionDialogs = new HashMap<Hub, TrackHubSelectionDialog>();

    private final List<TrackConfigContainer> trackConfigContainers;
    Hub hub;
    private ArrayList<CollapsiblePanel> categoryPanels;
    List<SelectionBox> allSelectionBoxes;
    boolean canceled = false;

    public static TrackHubSelectionDialog getTrackHubSelectionDialog(Hub hub, Set<String> loadedTrackPaths) {

        if (hubSelectionDialogs.containsKey(hub)) {
            TrackHubSelectionDialog dialog = hubSelectionDialogs.get(hub);
            dialog.resetSelectionBoxes(loadedTrackPaths);
            return dialog;
        } else {
            Frame owner = IGV.getInstance().getMainFrame();
            List<TrackConfigContainer> groups = hub.getGroupedTrackConfigurations();

            // Overide visibility
            if (loadedTrackPaths != null) {
                for (TrackConfigContainer g : groups) {
                    g.setTrackVisibility(loadedTrackPaths);
                }
            }

            TrackHubSelectionDialog dialog = new TrackHubSelectionDialog(hub, groups, owner);
            hubSelectionDialogs.put(hub, dialog);
            return dialog;
        }
    }

    public TrackHubSelectionDialog(Hub hub, List<TrackConfigContainer> trackConfigContainers, Frame owner) {
        super(owner);
        setModal(true);
        this.hub = hub;
        this.trackConfigContainers = trackConfigContainers;
        init(trackConfigContainers);
        setLocationRelativeTo(owner);
    }

    private void resetSelectionBoxes(Set<String> loadedTrackPaths) {
        if (loadedTrackPaths != null) {
            for (SelectionBox box : allSelectionBoxes) {
                box.setSelected(loadedTrackPaths.contains(box.trackConfig.getUrl()));
            }
        }
    }

    void init(List<TrackConfigContainer> trackConfigContainers) {

        setTitle(this.hub.getLongLabel());

        Rectangle ownerBounds = getOwner().getBounds();
        setSize(new Dimension(Math.min(ownerBounds.width, 1200), Math.min(ownerBounds.height, 1000)));

        categoryPanels = new ArrayList<>();
        allSelectionBoxes = new ArrayList<>();

        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BorderLayout());
        add(mainPanel);

        JPanel topPanel = new JPanel();
        topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));
        topPanel.add(getLabeledHyperlink("Hub URL: ", hub.getUrl()));
        String descriptionURL = hub.getDescriptionURL();

        if (descriptionURL != null) {
            topPanel.add(getLabeledHyperlink("Description: ", descriptionURL));
        }

        // Panel for select all/none
        JPanel topButtonPanel = new JPanel();
        topButtonPanel.setLayout(new BorderLayout());

        JPanel expandButtonPanel = new JPanel();
        expandButtonPanel.setLayout(new FlowLayout(FlowLayout.LEFT));
        JButton expandAllButton = new JButton("Expand All");
        expandAllButton.addActionListener(e -> {
            categoryPanels.forEach(cp -> cp.expand());
            this.revalidate();
        });
        expandButtonPanel.add(expandAllButton);

        JButton collapseAllButton = new JButton("Collapse All");
        collapseAllButton.addActionListener(e -> {
            categoryPanels.forEach(cp -> cp.collapse());
            this.revalidate();
        });
        expandButtonPanel.add(collapseAllButton);
        topButtonPanel.add(expandButtonPanel, BorderLayout.WEST);

        topPanel.add(topButtonPanel);
        mainPanel.add(topPanel, BorderLayout.NORTH);

        // Panel for category boxes
        JPanel categoryContainer = new JPanel();
        categoryContainer.setLayout(new BoxLayout(categoryContainer, BoxLayout.Y_AXIS));
        JScrollPane scrollPane = new JScrollPane(categoryContainer);
        mainPanel.add(scrollPane, BorderLayout.CENTER);

        // Loop through track groups
        for (TrackConfigContainer configGroup : trackConfigContainers) {
            categoryContainer.add(Box.createVerticalStrut(10));
            CollapsiblePanel categoryPanel = createCategoryPanel(configGroup);
            categoryContainer.add(categoryPanel);
            categoryPanels.add(categoryPanel);
        }

        // Search button.
        //JButton searchButton = createSearchButton("Search " + hub.getShortLabel(), allSelectionBoxes);
        //topButtonPanel.add(searchButton, BorderLayout.EAST);

        JPanel buttonPanel = new JPanel();
        ((FlowLayout) buttonPanel.getLayout()).setAlignment(FlowLayout.RIGHT);

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> {
            canceled = true;
            setVisible(false);
        });

        JButton okButton = new JButton("OK");
        okButton.addActionListener(e -> {
            canceled = false;
            setVisible(false);
        });

//        if (Globals.IS_MAC) {
//            buttonPanel.add(cancelButton);
//            buttonPanel.add(okButton);
//        } else {
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);
        //}

        getRootPane().setDefaultButton(okButton);

        mainPanel.add(buttonPanel, BorderLayout.SOUTH);

        // mainPanel.validate();

    }

    private JPanel getLabeledHyperlink(String label, String url) {
        JPanel hubURLPanel = new JPanel();
        ((FlowLayout) hubURLPanel.getLayout()).setAlignment(FlowLayout.LEFT);
        hubURLPanel.setBorder(BorderFactory.createEmptyBorder(2, 6, 0, 0));
        hubURLPanel.add(new JLabel(label));
        hubURLPanel.add(HyperlinkFactory.createLink(url, url));
        return hubURLPanel;
    }


    public boolean isCanceled() {
        return canceled;
    }

    /**
     * Create a JPanel for a particular track category, including a label and checkboxes for contained tracks.
     *
     * @param configGroup
     * @return
     */
    private CollapsiblePanel createCategoryPanel(TrackConfigContainer configGroup) {

        JPanel trackContainer = new JPanel();
        trackContainer.setLayout(new BoxLayout(trackContainer, BoxLayout.Y_AXIS));

        // There is a (so far) intractable bug if a large # of JCheckboxes are created for this widget.
        int totalTrackCount = configGroup.countTracks();
        SelectionBox.CheckboxType checkboxType = totalTrackCount < 1000 ? SelectionBox.CheckboxType.SWING : SelectionBox.CheckboxType.CUSTOM;

        List<SelectionBox> selectionBoxes = addSelectionBoxes(null, configGroup, trackContainer, checkboxType);

        boolean isSelected = false;
        int maxWidth = 0;
        int selectionCount = 0;
        for (SelectionBox selectionBox : selectionBoxes) {
            if (selectionBox.checkbox.isSelected()) {
                selectionCount++;
                isSelected = true;
            }
            maxWidth = Math.max(maxWidth, selectionBox.getPreferredSize().width);
            allSelectionBoxes.add(selectionBox);
        }

        for (SelectionBox selectionBox : selectionBoxes) {
            selectionBox.setPreferredWidth(maxWidth);
        }

        String label = configGroup.label + "   (" + selectionBoxes.size() + " tracks, " + selectionCount + " selected)";

        final CollapsiblePanel collapsiblePanel = new CollapsiblePanel(label, trackContainer, isSelected || configGroup.defaultOpen);

        // Add a search button for categories with large numbers of records
        final JButton searchButton = createSearchButton("Search " + configGroup.label, selectionBoxes,
                (selectedCount) -> {
                    String l = configGroup.label + "   (" + selectionBoxes.size() + " tracks, " + selectedCount + " selected)";
                    collapsiblePanel.updateLabel(l);
                    return null;
                });

        collapsiblePanel.addSearchButton(searchButton);

        for (SelectionBox selectionBox : selectionBoxes) {
            selectionBox.setCallback(b -> {
                int selected = configGroup.countSelectedTracks();
                String l = configGroup.label + "   (" + selectionBoxes.size() + " tracks, " + selected + " selected)";
                collapsiblePanel.updateLabel(l);
                return null;

            });
        }

        return collapsiblePanel;
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
    private List<SelectionBox> addSelectionBoxes(String labelPrefix, TrackConfigContainer container, JPanel panel, SelectionBox.CheckboxType checkboxType) {

        String title = labelPrefix == null ? "" :
                labelPrefix + (labelPrefix.length() > 0 ? " - " : "") + container.label;

        List<SelectionBox> selectionBoxes = new ArrayList<>();

        if (container.tracks.size() > 0) {

            JPanel trackPanel = new JPanel();
            if (labelPrefix != null) {
                trackPanel.setBorder(BorderFactory.createTitledBorder(title));
            }
            final WrapLayout wrapLayout = new WrapLayout();
            wrapLayout.setAlignment(FlowLayout.LEFT);
            trackPanel.setLayout(wrapLayout);

            for (TrackConfig trackConfig : container.tracks) {
                SelectionBox p = new SelectionBox(trackConfig, checkboxType);
                trackPanel.add(p);
                selectionBoxes.add(p);
            }

            panel.add(Box.createVerticalStrut(5));
            panel.add(trackPanel);

            // final CollapsiblePanel collapsiblePanel = new CollapsiblePanel(title, trackPanel, false, CollapsiblePanel.HEADER_BG2);
            // trackContainer.add(collapsiblePanel);
        }

        for (TrackConfigContainer childChild : container.children) {
            selectionBoxes.addAll(addSelectionBoxes(title, childChild, panel, checkboxType));
        }
        return selectionBoxes;
    }

    private JButton createSearchButton(String label, List<SelectionBox> selectionBoxes, Function<Integer, Void> callback) {

        JButton searchButton = new JButton("Search");

        searchButton.addActionListener(e -> {

            Set<String> attributeNames = new LinkedHashSet<>();
            attributeNames.add("Name");
            attributeNames.add("Description");
            attributeNames.add("Format");


            Map<FileRecord, SelectionBox> recordSelectionBoxMap = new HashMap<>();

            List<FileRecord> records = new ArrayList<>();

            for (SelectionBox selectionBox : selectionBoxes) {

                TrackConfig trackConfig = selectionBox.getTrackConfig();
                final Map<String, String> trackConfigAttributes = trackConfig.getAttributes();
                Map<String, String> attributes = trackConfigAttributes;
                if (attributes == null) {
                    attributes = new LinkedHashMap<>();
                }

                attributes.put("Name", trackConfig.getName());
                attributes.put("Description", trackConfig.getDescription());
                attributes.put("Format", trackConfig.getFormat());

                if (trackConfigAttributes != null) {
                    attributes.putAll(trackConfigAttributes);
                    attributeNames.addAll(trackConfigAttributes.keySet());
                }

                final FileRecord record = new FileRecord(trackConfig.getUrl(), attributes);
                record.setSelected(trackConfig.getVisible());
                records.add(record);
                recordSelectionBoxMap.put(record, selectionBox);
            }


            List<String> headings = new ArrayList<>(attributeNames);
            // Limit # of columns
            if (headings.size() > 15) {
                headings = headings.subList(0, 15);
            }

            Frame owner = IGV.hasInstance() ? IGV.getInstance().getMainFrame() : null;

            TrackChooser chooser = new TrackChooser(
                    owner,
                    headings,
                    records,
                    label);
            chooser.setSize(this.getSize());
            chooser.setLocationRelativeTo(getOwner());
            chooser.setVisible(true);

            if (!chooser.isCanceled()) {
                Set<FileRecord> selectedRecords = new HashSet<>(chooser.getSelectedRecords());
                for (Map.Entry<FileRecord, SelectionBox> entry : recordSelectionBoxMap.entrySet()) {
                    entry.getValue().setSelected(selectedRecords.contains(entry.getKey()));
                }
                callback.apply(selectedRecords.size());
            }
        });
        return searchButton;
    }

    /**
     * Convenience method to extract and return selected track configurations.
     *
     * @return
     */
    public List<TrackConfig> getSelectedConfigs() {

        List<TrackConfig> selectedConfigs = new ArrayList<>();
        for (TrackConfigContainer container : trackConfigContainers) {
            container.findSelectedConfigs(selectedConfigs);
        }
        return selectedConfigs;
    }

    static class SelectionBox extends JPanel {

        enum CheckboxType {SWING, CUSTOM}

        TrackConfig trackConfig;
        private CheckBoxWrapper checkbox;
        int preferredWidth = -1;
        private int minWidth;
        Function<Integer, Void> callback;

        public SelectionBox(TrackConfig trackConfig, CheckboxType checkboxType) {

            this.setLayout(new BorderLayout(5, 0));
            this.trackConfig = trackConfig;

            String longLabel = trackConfig.getLongLabel();
            if (longLabel != null) {
                this.setToolTipText(longLabel);
            }

            this.checkbox = new CheckBoxWrapper(checkboxType);
            checkbox.setSelected(trackConfig.getVisible());
            checkbox.setActionListener(e -> {
                trackConfig.setVisible(checkbox.isSelected());
                if (callback != null) {
                    callback.apply(checkbox.isSelected() ? 1 : 0);
                }
            });

            JLabel label = new JLabel(trackConfig.getName());
            label.setLabelFor(checkbox.getComponent());
            add(checkbox.getComponent(), BorderLayout.WEST);

            String infoLink = trackConfig.getHtml();

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
            trackConfig.setVisible(selected);
        }

        public TrackConfig getTrackConfig() {
            return trackConfig;
        }

        public void setCallback(Function<Integer, Void> callback) {
            this.callback = callback;
        }
    }

    static class CheckBoxWrapper {

        CheckBox checkBox;
        JCheckBox jCheckBox;

        public CheckBoxWrapper(SelectionBox.CheckboxType checkboxType) {
            if(checkboxType == SelectionBox.CheckboxType.SWING) {
                jCheckBox = new JCheckBox();
            } else {
                checkBox = new CheckBox();
            }
        }

        public JComponent getComponent() {
            return checkBox != null ? checkBox : jCheckBox;
        }

        public void setSelected(boolean selected) {
            if (checkBox != null) {
                checkBox.setSelected(selected);
            } else {
                jCheckBox.setSelected(selected);
            }
        }

        public boolean isSelected() {
            return checkBox != null ? checkBox.isSelected() :  jCheckBox.isSelected();
        }

        public void setActionListener(ActionListener l) {
            if(checkBox != null) {
                checkBox.setActionListener(l);
            } else {
                jCheckBox.addActionListener(l);
            }
        }
    }

    static class CheckBox extends JLabel {

        boolean selected = false;
        Icon checkedIcon;
        Icon uncheckedIcon;
        private ActionListener actionListener;

        public CheckBox() {

            this.checkedIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.CHECKBOX);
            this.uncheckedIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.CHECKBOX_UNCHECKED);
            setIcon(selected ? checkedIcon : uncheckedIcon);
            this.setVerticalAlignment(SwingConstants.BOTTOM);
            setSize(new Dimension(16, 16));

            addMouseListener(new MouseAdapter() {
                @Override
                public void mouseReleased(MouseEvent e) {
                    setSelected(!selected);
                    if (actionListener != null) {
                        actionListener.actionPerformed(new ActionEvent(this, 0, ""));
                    }
                }
            });
        }

        public void setSelected(boolean selected) {
            this.selected = selected;
            setIcon(selected ? checkedIcon : uncheckedIcon);
        }

        public boolean isSelected() {
            return selected;
        }

        public void setActionListener(ActionListener l) {
            this.actionListener = l;
        }

        @Override
        public Dimension getPreferredSize() {
            return new Dimension(16, 16);
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

        String hubFile = "https://hgdownload.soe.ucsc.edu/gbdb/hs1/hubs/public/hub.txt";

        Hub hub = HubParser.loadHub(hubFile, "hs1");

        List<TrackConfigContainer> groupedTrackConfigurations = hub.getGroupedTrackConfigurations();

        final TrackHubSelectionDialog dlf = new TrackHubSelectionDialog(hub, groupedTrackConfigurations, null);
        dlf.setSize(new Dimension(800, 600));
        dlf.setVisible(true);

        for (TrackConfig config : dlf.getSelectedConfigs()) {
            System.out.println(config.getName());
        }
    }

}
