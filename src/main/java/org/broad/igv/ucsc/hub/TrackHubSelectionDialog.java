package org.broad.igv.ucsc.hub;

import org.broad.igv.Globals;
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
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.net.URI;
import java.util.*;
import java.util.List;

import java.util.stream.Collectors;


/**
 * Dialog to enable selection of tracks defined by track hubs.  Modifies the "visible" property of
 * supplied track configurations in place.
 */
public class TrackHubSelectionDialog extends JDialog {

    private static Logger log = LogManager.getLogger(TrackHubSelectionDialog.class);

    private final List<TrackConfigGroup> trackConfigGroups;
    Hub hub;
    private ArrayList<CollapsiblePanel> categoryPanels;
    List<SelectionBox> allSelectionBoxes;
    boolean canceled = false;


    public TrackHubSelectionDialog(Hub hub, List<TrackConfigGroup> trackConfigGroups, Frame owner) {
        super(owner);
        setModal(true);
        this.hub = hub;
        this.trackConfigGroups = trackConfigGroups;
        init(trackConfigGroups);
        setLocationRelativeTo(owner);
    }


    void init(List<TrackConfigGroup> trackConfigGroups) {

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
        //scrollPane.setBorder(BorderFactory.createLineBorder(Color.gray));
        mainPanel.add(scrollPane, BorderLayout.CENTER);

        // Loop through track groups
        for (TrackConfigGroup configGroup : trackConfigGroups) {
            categoryContainer.add(Box.createVerticalStrut(10));
            CollapsiblePanel categoryPanel = createCategoryPanel(configGroup);
            categoryContainer.add(categoryPanel);
            categoryPanels.add(categoryPanel);
        }

        // If only a single category expand it


        // If total # of tracks is small expand all
        if (allSelectionBoxes.size() < 50) {
            for (CollapsiblePanel panel : categoryPanels) {
                panel.expand();
            }
        } else {
            if (categoryPanels.size() == 1) {
                categoryPanels.get(0).expand();
            }
        }

        // Search button.
        JButton searchButton = createSearchButton("Search " + hub.getShortLabel() , allSelectionBoxes);
        topButtonPanel.add(searchButton, BorderLayout.EAST);

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

        mainPanel.validate();

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
    private CollapsiblePanel createCategoryPanel(TrackConfigGroup configGroup) {

        JPanel trackContainer = new JPanel();
        final WrapLayout wrapLayout = new WrapLayout();
        wrapLayout.setAlignment(FlowLayout.LEFT);
        trackContainer.setLayout(wrapLayout);

        int maxWidth = 0;
        List<SelectionBox> selectionBoxes = new ArrayList<>();
        boolean isSelected = false;
        for (TrackConfig trackConfig : configGroup.tracks) {

            isSelected = isSelected || trackConfig.getVisible();

            SelectionBox p = new SelectionBox(trackConfig);
            selectionBoxes.add(p);
            allSelectionBoxes.add(p);

            maxWidth = Math.max(maxWidth, p.getPreferredSize().width);

            String longLabel = trackConfig.getLongLabel();
            if (longLabel != null) {
                p.setToolTipText(longLabel);
            }

            trackContainer.add(p);
        }

        for (SelectionBox selectionBox : selectionBoxes) {
            selectionBox.setPreferredWidth(maxWidth);
        }

        final CollapsiblePanel collapsiblePanel = new CollapsiblePanel(configGroup.label, trackContainer, isSelected || configGroup.defaultOpen);

        // Add a search button for categories with large numbers of records

        if (configGroup.tracks.size() > 2) {
            final JButton searchButton = createSearchButton("Search " + configGroup.label, selectionBoxes);
            searchButton.addActionListener(e -> collapsiblePanel.expand());
            collapsiblePanel.addSearchButton(searchButton);
        }


        return collapsiblePanel;
    }

    private JButton createSearchButton(String label, List<SelectionBox> selectionBoxes) {

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
        return trackConfigGroups.stream()
                .flatMap(group -> group.tracks.stream())
                .filter(trackConfig -> trackConfig.getVisible())
                .collect(Collectors.toList());
    }

    static class SelectionBox extends JPanel {

        TrackConfig trackConfig;
        private JCheckBox checkbox;
        int preferredWidth = -1;
        private int minWidth;


        public SelectionBox(TrackConfig trackConfig) {

            this.setLayout(new BorderLayout());
            this.trackConfig = trackConfig;

            this.checkbox = new JCheckBox();
            checkbox.setSelected(trackConfig.getVisible());
            checkbox.addActionListener(e -> {
                trackConfig.setVisible(checkbox.isSelected());
            });

            JLabel label = new JLabel(trackConfig.getName());
            label.setLabelFor(checkbox);
            add(checkbox, BorderLayout.WEST);

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

        public boolean isSelected() {
            return checkbox.isSelected();
        }

        public TrackConfig getTrackConfig() {
            return trackConfig;
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

        List<TrackConfigGroup> groupedTrackConfigurations = hub.getGroupedTrackConfigurations();

        final TrackHubSelectionDialog dlf = new TrackHubSelectionDialog(hub, groupedTrackConfigurations, null);
        dlf.setVisible(true);

        for (TrackConfig config : dlf.getSelectedConfigs()) {
            System.out.println(config.getName());
        }
    }

}
