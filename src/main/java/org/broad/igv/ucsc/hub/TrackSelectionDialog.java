package org.broad.igv.ucsc.hub;

import org.broad.igv.encode.FileRecord;
import org.broad.igv.encode.TrackChooser;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.FontManager;
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
import java.util.*;
import java.util.List;
import java.util.function.Function;


/**
 * Dialog to enable selection of tracks defined by track hubs.  Modifies the "visible" property of
 * supplied track configurations.
 */
public class TrackSelectionDialog extends JDialog {

    private static Logger log = LogManager.getLogger(TrackSelectionDialog.class);

    private static Map<Hub, TrackSelectionDialog> hubSelectionDialogs = new HashMap<>();

    private Hub hub;
    private boolean autoselectDefaults;
    private ArrayList<CollapsiblePanel> categoryPanels;
    private boolean canceled = true;

    public static TrackSelectionDialog getTrackHubSelectionDialog(
            Hub hub,
            String genomeId,
            Set<String> loadedTrackPaths,
            boolean autoselectDefaults,
            String message) {

        TrackSelectionDialog dialog;
        if (hubSelectionDialogs.containsKey(hub) && !autoselectDefaults) {
            dialog = hubSelectionDialogs.get(hub);
            dialog.autoselectDefaults = autoselectDefaults;
        } else {
            Frame owner = IGV.getInstance().getMainFrame();
            List<TrackConfigContainer> groups = hub.getGroupedTrackConfigurations(genomeId);
            dialog = new TrackSelectionDialog(hub, groups, autoselectDefaults, message, owner);
            hubSelectionDialogs.put(hub, dialog);
        }
        dialog.resetSelectionBoxes(loadedTrackPaths);

        return dialog;
    }

    private TrackSelectionDialog(Hub hub,
                                 List<TrackConfigContainer> trackConfigContainers,
                                 boolean autoselectDefaults,
                                 String message,
                                 Frame owner) {
        super(owner);
        setModal(true);
        this.autoselectDefaults = autoselectDefaults;
        this.hub = hub;
        init(trackConfigContainers, message);
        setLocationRelativeTo(owner);
    }

    /**
     * Called when hub dialog is reused.  Update the selection box state to reflect currently loaded tracks, which
     * could have changed since last invocation.
     *
     * @param loadedTrackPaths
     */
    private void resetSelectionBoxes(Set<String> loadedTrackPaths) {
        for (CollapsiblePanel collapsiblePanel : categoryPanels) {
            collapsiblePanel.resetSelectionBoxes(loadedTrackPaths);
        }
    }

    void init(List<TrackConfigContainer> trackConfigContainers, String message) {

        setTitle(this.hub.getLongLabel());

        Rectangle ownerBounds = getOwner().getBounds();
        setSize(new Dimension(Math.min(ownerBounds.width, 1200), Math.min(ownerBounds.height, 1000)));

        categoryPanels = new ArrayList<>();

        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BorderLayout());
        add(mainPanel);

        JPanel topPanel = new JPanel();
        topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));

        if (message != null) {
            JTextPane messagePane = new JTextPane();
            messagePane.setBorder(BorderFactory.createEmptyBorder(30, 5, 20, 5));
            messagePane.setFont(FontManager.getFont(Font.BOLD, 16));
            messagePane.setText(message);
            messagePane.setEditable(false);
            topPanel.add(messagePane);
        }

        topPanel.add(getLabeledHyperlink("Hub URL: ", hub.getUrl()));
        String descriptionURL = hub.getDescriptionURL();

        if (descriptionURL != null) {
            topPanel.add(getLabeledHyperlink("Description: ", descriptionURL));
        }

        // Panel for select all/none
        JPanel topButtonPanel = new JPanel();
        topButtonPanel.setLayout(new BorderLayout());

        JPanel expandButtonPanel = new JPanel();
        expandButtonPanel.setLayout(new FlowLayout(FlowLayout.RIGHT));
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

//        JPanel clearAllPanel = new JPanel();
//        clearAllPanel.setLayout(new FlowLayout(FlowLayout.RIGHT));
//        JButton clearAllButton = new JButton("Clear All");
//        clearAllPanel.add(clearAllButton);
//        clearAllButton.addActionListener(e -> {
//            for (CollapsiblePanel collapsiblePanel : categoryPanels) {
//                collapsiblePanel.clearSelections();
//            }
//        });
//        topButtonPanel.add(clearAllPanel, BorderLayout.EAST);

        topPanel.add(topButtonPanel);
        mainPanel.add(topPanel, BorderLayout.NORTH);

        // Panel for category boxes
        JPanel categoryContainer = new JPanel();
        categoryContainer.setLayout(new BoxLayout(categoryContainer, BoxLayout.Y_AXIS));
        JScrollPane scrollPane = new JScrollPane(categoryContainer);
        mainPanel.add(scrollPane, BorderLayout.CENTER);

        // Get total track count
        int count = 0;
        for (TrackConfigContainer cp : trackConfigContainers) {
            count += cp.countTracks();
        }

        // Loop through track groups
        boolean catSearch = count > 10000;
        for (TrackConfigContainer configGroup : trackConfigContainers) {
            categoryContainer.add(Box.createVerticalStrut(10));
            CollapsiblePanel categoryPanel = createCategoryPanel(configGroup, catSearch);
            categoryContainer.add(categoryPanel);
            categoryPanels.add(categoryPanel);
        }

        // Search button.
        if (!catSearch) {
            JButton searchButton = createSearchButton("Search " + hub.getShortLabel(), categoryPanels, (selectedCount) -> {
                categoryPanels.stream().forEach(p -> p.updateLabel());
                return null;
            });
            topButtonPanel.add(searchButton, BorderLayout.EAST);
        }

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
    private CollapsiblePanel createCategoryPanel(TrackConfigContainer configGroup, boolean search) {

        boolean autoselectGenes = this.autoselectDefaults && configGroup.name.toLowerCase().equals("genes");

        return new CollapsiblePanel(configGroup, autoselectGenes, search);
    }


    /**
     * Convenience method to extract and return selected track configurations.
     *
     * @return
     */
    public List<TrackConfig> getSelectedConfigs() {

        List<TrackConfig> selectedConfigs = new ArrayList<>();
        for (CollapsiblePanel collapsiblePanel : categoryPanels) {
            selectedConfigs.addAll(collapsiblePanel.getSelectedTracks());
        }
        return selectedConfigs;
    }

    static class CheckBoxWrapper {

        CheckBox checkBox;
        JCheckBox jCheckBox;

        public CheckBoxWrapper(CollapsiblePanel.SelectionBox.CheckboxType checkboxType) {
            if (checkboxType == CollapsiblePanel.SelectionBox.CheckboxType.SWING) {
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
            return checkBox != null ? checkBox.isSelected() : jCheckBox.isSelected();
        }

        public void setActionListener(ActionListener l) {
            if (checkBox != null) {
                checkBox.setActionListener(l);
            } else {
                jCheckBox.addActionListener(l);
            }
        }

        public void setEnabled(boolean enabled) {
            if (checkBox != null) {
                checkBox.setEnabled(enabled);
            } else {
                jCheckBox.setEnabled(enabled);
            }
        }

        public boolean isEnabled() {
            return checkBox != null ? checkBox.isEnabled() : jCheckBox.isEnabled();
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

    public static JButton createSearchButton(String label, List<CollapsiblePanel> panels, Function<Integer, Void> callback) {

        JButton searchButton = new JButton("Search");

        searchButton.addActionListener(e -> {

            Set<String> attributeNames = new LinkedHashSet<>();
            if (panels.size() > 1) {
                attributeNames.add("Group");
            }
            attributeNames.add("Name");
            attributeNames.add("Description");
            attributeNames.add("Format");

            Map<FileRecord, CollapsiblePanel.SelectionBox> recordSelectionBoxMap = new HashMap<>();

            List<FileRecord> records = new ArrayList<>();

            for (CollapsiblePanel panel : panels) {

                for (CollapsiblePanel.SelectionBox selectionBox : panel.selectionBoxes) {

                    if (selectionBox.isEnabled()) {

                        TrackConfig trackConfig = selectionBox.getTrackConfig();
                        final Map<String, String> trackConfigAttributes = trackConfig.attributes;
                        Map<String, String> attributes = trackConfigAttributes;
                        if (attributes == null) {
                            attributes = new LinkedHashMap<>();
                        }
                        attributes.put("Group", panel.containerLabel());
                        attributes.put("Name", trackConfig.name);
                        attributes.put("Description", trackConfig.description);
                        attributes.put("Format", trackConfig.format);

                        if (trackConfigAttributes != null) {
                            attributes.putAll(trackConfigAttributes);
                            attributeNames.addAll(trackConfigAttributes.keySet());
                        }

                        final FileRecord record = new FileRecord(trackConfig.url, attributes);
                        record.setSelected(trackConfig.visible);
                        records.add(record);
                        recordSelectionBoxMap.put(record, selectionBox);
                    }
                }
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

            if (owner != null) {
                Rectangle ownerBounds = owner.getBounds();
                chooser.setSize(new Dimension(Math.min(ownerBounds.width, 1200), Math.min(ownerBounds.height, 800)));
                chooser.setLocationRelativeTo(owner);
            }
            chooser.setVisible(true);

            if (!chooser.isCanceled()) {
                Set<FileRecord> selectedRecords = new HashSet<>(chooser.getSelectedRecords());
                for (Map.Entry<FileRecord, CollapsiblePanel.SelectionBox> entry : recordSelectionBoxMap.entrySet()) {
                    entry.getValue().setSelected(selectedRecords.contains(entry.getKey()));
                }
                callback.apply(selectedRecords.size());
            }
        });
        return searchButton;
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

        Hub hub = HubParser.loadHub(hubFile);

        List<TrackConfigContainer> groupedTrackConfigurations = hub.getGroupedTrackConfigurations("hs1");

        final TrackSelectionDialog dlf = new TrackSelectionDialog(hub, groupedTrackConfigurations, true, null, null);
        dlf.setSize(new Dimension(800, 600));
        dlf.setVisible(true);

        for (TrackConfig config : dlf.getSelectedConfigs()) {
            System.out.println(config.name);
        }
    }

}
