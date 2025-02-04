package org.broad.igv.ucsc.hub;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.ui.panel.CollapsiblePanel;
import org.broad.igv.ui.util.HyperlinkFactory;

import javax.swing.*;
import java.awt.*;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;


/**
 * Dialog to enable selection of tracks defined by track hubs.  Modifies the "visible" property of
 * supplied track configurations in place.
 */
public class TrackHubSelectionDialog extends JDialog {

    Hub hub;
    private Map<JCheckBox, TrackConfig> configMap;
    private ArrayList<CollapsiblePanel> categoryPanels;

    public TrackHubSelectionDialog(Hub hub, List<TrackConfigGroup> groupedTrackConfigurations, Frame owner) {
        super(owner);
        setModal(true);
        this.hub = hub;
        init(groupedTrackConfigurations);
        setLocationRelativeTo(owner);
    }


    void init(List<TrackConfigGroup> trackConfigurations) {

        setTitle(this.hub.getLongLabel());

        setSize(new Dimension(1000, 800));

        configMap = new HashMap<>();
        categoryPanels = new ArrayList<>();

        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BorderLayout());
        add(mainPanel);

        JPanel topPanel = new JPanel();
        topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));
        topPanel.add(getLabeledHyperlink("Hub URL: ", hub.getUrl()));
        String descriptionURL = hub.getDescriptionURL();
        if(descriptionURL != null) {
            topPanel.add(getLabeledHyperlink("Description: ", descriptionURL));
        }

        // Panel for select all/none
        JPanel topButtonPanel = new JPanel();
        ((FlowLayout) topButtonPanel.getLayout()).setAlignment(FlowLayout.LEFT);
        JButton expandAllButton = new JButton("Expand All");
        expandAllButton.addActionListener(e -> {
            categoryPanels.forEach(cp -> cp.expand());
            this.revalidate();
        });
        topButtonPanel.add(expandAllButton);

        JButton collapseAllButton = new JButton("Collapse All");
        collapseAllButton.addActionListener(e -> {
            categoryPanels.forEach(cp -> cp.collapse());
            this.revalidate();
        });
        topButtonPanel.add(collapseAllButton);

        topPanel.add(topButtonPanel);
        mainPanel.add(topPanel, BorderLayout.NORTH);

        // Panel for category boxes
        JPanel categoryContainer = new JPanel();
        categoryContainer.setLayout(new BoxLayout(categoryContainer, BoxLayout.PAGE_AXIS));
        JScrollPane scrollPane = new JScrollPane(categoryContainer);
        scrollPane.setBorder(BorderFactory.createLineBorder(Color.gray));
        mainPanel.add(scrollPane, BorderLayout.CENTER);

        // Loop through track groups
        for (TrackConfigGroup configGroup : trackConfigurations) {
            categoryContainer.add(Box.createVerticalStrut(10));
            CollapsiblePanel categoryPanel = createCategoryPanel(configGroup);
            categoryContainer.add(categoryPanel);
            categoryPanels.add(categoryPanel);
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

    }

    private JPanel getLabeledHyperlink(String label, String url) {
        JPanel hubURLPanel = new JPanel();
        ((FlowLayout) hubURLPanel.getLayout()).setAlignment(FlowLayout.LEFT);
        hubURLPanel.setBorder(BorderFactory.createEmptyBorder(2, 6, 0, 0));
        hubURLPanel.add(new JLabel(label));
        hubURLPanel.add(HyperlinkFactory.createLink(url, url));
        return hubURLPanel;
    }

    private void okAction() {
        for (Map.Entry<JCheckBox, TrackConfig> entry : configMap.entrySet()) {
            final TrackConfig trackConfig = entry.getValue();
            if (entry.getKey().isSelected()) {
                trackConfig.setVisible(true);
            } else {
                trackConfig.setVisible(false);
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
    private CollapsiblePanel createCategoryPanel(TrackConfigGroup configGroup) {

        JPanel trackContainer = new JPanel();
        final WrapLayout wrapLayout = new WrapLayout();
        wrapLayout.setAlignment(FlowLayout.LEFT);
        trackContainer.setLayout(wrapLayout);

        boolean isSelected = false;
        for (TrackConfig trackConfig : configGroup.tracks) {

            final JCheckBox checkBox = new JCheckBox();
            configMap.put(checkBox, trackConfig);
            checkBox.setSelected(trackConfig.getVisible());
            isSelected = isSelected || trackConfig.getVisible();

//            JLabel label = trackConfig.getHtml() == null ?
//                    new JLabel(trackConfig.getName()) :
//                    HyperlinkFactory.createLink(trackConfig.getName(), trackConfig.getHtml());
            JLabel label = new JLabel(trackConfig.getName());


            SelectionBox p = new SelectionBox(checkBox, label);
            trackContainer.add(p);
        }

        return new CollapsiblePanel(configGroup.label, trackContainer, isSelected || configGroup.defaultOpen);
    }

    /**
     * Convenience method to extract and return selected track configurations.
     *
     * @return
     */
    public List<TrackConfig> getSelectedConfigs() {
        List<TrackConfig> selected = configMap.values().stream().filter(trackConfig -> trackConfig.getVisible()).collect(Collectors.toList());
        return selected;
    }

    static class SelectionBox extends JPanel {

        public SelectionBox(JCheckBox checkBox, JLabel label) {
            this.setLayout(new BorderLayout());
            label.setLabelFor(checkBox);
            add(checkBox, BorderLayout.WEST);
            add(label, BorderLayout.CENTER);

        }

        @Override
        public Dimension getPreferredSize() {
            return new Dimension(200, 20);
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
        Hub hub = Hub.loadHub(hubFile);
        List<TrackConfigGroup> groupedTrackConfigurations = hub.getGroupedTrackConfigurations();

        final TrackHubSelectionDialog dlf = new TrackHubSelectionDialog(hub, groupedTrackConfigurations, null);
        dlf.setVisible(true);

        for (TrackConfig config : dlf.getSelectedConfigs()) {
            System.out.println(config.getName());
        }
    }

}
