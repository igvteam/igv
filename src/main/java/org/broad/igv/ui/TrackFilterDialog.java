package org.broad.igv.ui;

import org.broad.igv.track.AttributeManager;
import org.broad.igv.util.Filter;
import org.broad.igv.util.FilterElement;

import javax.swing.*;
import java.awt.*;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class TrackFilterDialog extends JDialog {

    JPanel trackFilterPane;
    JScrollPane scrollPane;

    private JCheckBox showAllTracksFilterCheckBox;
    private JCheckBox matchAllCheckBox;
    private JCheckBox matchAnyCheckBox;

    private boolean cancelled = true;
    private List<String> attributeKeys;
    private final int vgap = 5;

    public TrackFilterDialog(Frame owner, String title, Filter filter) {
        super(owner, title, true);
        setLocationRelativeTo(owner);
        setResizable(true);
        init(filter);
        setSize(700, 600);
    }

    private void init(Filter filter) {

        if (filter == null) {
            filter = new Filter();
        }

        setLayout(new BorderLayout());

        JPanel filterHeaderPanel = new JPanel();
        filterHeaderPanel.setLayout(new GridLayout(0, 1));

        showAllTracksFilterCheckBox = new JCheckBox("Show all tracks");
        showAllTracksFilterCheckBox.setSelected(filter.isShowAll());
        showAllTracksFilterCheckBox.addActionListener(e -> {
            matchAllCheckBox.setEnabled(!showAllTracksFilterCheckBox.isSelected());
            matchAnyCheckBox.setEnabled(!showAllTracksFilterCheckBox.isSelected());
        });
        filterHeaderPanel.add(showAllTracksFilterCheckBox);

        matchAllCheckBox = new JCheckBox("match all of the following");
        matchAnyCheckBox = new JCheckBox("match any of the following");
        matchAllCheckBox.setSelected(filter.isMatchAll());
        matchAnyCheckBox.setSelected(!filter.isMatchAll());
        ButtonGroup booleanButtonGroup = new ButtonGroup();
        booleanButtonGroup.add(matchAllCheckBox);
        booleanButtonGroup.add(matchAnyCheckBox);


        if (filter.isShowAll()) {
            matchAllCheckBox.setEnabled(false);
            matchAnyCheckBox.setEnabled(false);
        } else {
            matchAllCheckBox.setEnabled(true);
            matchAnyCheckBox.setEnabled(true);
        }

        JPanel controls = new JPanel();
        FlowLayout layoutManager = new FlowLayout();
        layoutManager.setAlignment(FlowLayout.LEFT);
        controls.setLayout(layoutManager);
        controls.add(new JLabel("Show tracks that: "));
        controls.add(matchAllCheckBox);
        controls.add(matchAnyCheckBox);

        controls.add(Box.createHorizontalStrut(10));
        JButton addButton = new JButton("+");
        addButton.addActionListener(e -> {
            addComponent();
        });
        controls.add(addButton);

        filterHeaderPanel.add(controls);

        add(filterHeaderPanel, BorderLayout.NORTH);

        List<String> uniqueAttributeKeys = AttributeManager.getInstance().getAttributeNames();
        this.attributeKeys = uniqueAttributeKeys;

        trackFilterPane = new JPanel();
        trackFilterPane.setLayout(new FlowLayout(FlowLayout.LEFT));

        scrollPane = new JScrollPane(trackFilterPane);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        //scrollPane.setPreferredSize(new Dimension(700, 500));

        if(filter.getFilterElements().hasNext()) {
            Iterator iterator = filter.getFilterElements();
            while (iterator.hasNext()) {
                FilterElement element = (FilterElement) iterator.next();
                trackFilterPane.add(new FilterComponent(this, uniqueAttributeKeys, element));
            }
        } else {
        // If there are no filters, add a blank one
            trackFilterPane.add(new FilterComponent(this, uniqueAttributeKeys, null));
        }

        add(scrollPane, BorderLayout.CENTER);

        // Button bar with OK and Cancel buttons
        JPanel buttonBar = new JPanel(new FlowLayout(FlowLayout.RIGHT));
        JButton okButton = new JButton("OK");
        JButton cancelButton = new JButton("Cancel");

        okButton.addActionListener(e -> {
            cancelled = false;
            setVisible(false);
            dispose();
        });

        cancelButton.addActionListener(e -> {
            cancelled = true;
            setVisible(false);
            dispose();
        });

        buttonBar.add(okButton);
        buttonBar.add(cancelButton);

        add(buttonBar, BorderLayout.SOUTH);

    }

    /**
     * Add another filter
     */
    public void addComponent() {

        trackFilterPane.add(new FilterComponent(this, attributeKeys, null));
        resizeFilterPane();
    }

    public void removeComponent(FilterComponent filterComponent) {

        if (trackFilterPane.getComponents().length < 2) {
            // Can not leave less than one filter element on the screen
            return;
        }

        trackFilterPane.remove(filterComponent);
        resizeFilterPane();

    }

    private void resizeFilterPane() {
        int h = 0;
        int w = 0;
        for (Component c : trackFilterPane.getComponents()) {
            h += c.getHeight() + 2*vgap;
            w = c.getWidth();
        }
        trackFilterPane.setPreferredSize(new Dimension(w, h));
        scrollPane.getViewport().revalidate();

    }




    public boolean isCancelled() {
        return cancelled;
    }

    public Filter getFilter() {
        boolean matchAll = matchAllCheckBox.isSelected();
        boolean showAllTracks = showAllTracksFilterCheckBox.isSelected();
        List<FilterElement> filterElements = new ArrayList<>();
        for (Component component : trackFilterPane.getComponents()) {
            if (component instanceof FilterComponent) {
                if (((FilterComponent) component).isComplete()) {
                    filterElements.add(((FilterComponent) component).getFilterElement(matchAll));
                }
            }
        }
        return new Filter(showAllTracks, matchAll, filterElements);
    }

}
