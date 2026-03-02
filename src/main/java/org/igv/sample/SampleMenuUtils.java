package org.igv.sample;

import org.igv.track.AbstractTrack;
import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.ui.AttributeSelectionDialog;
import org.igv.ui.IGV;
import org.igv.ui.SampleFilterDialog;
import org.igv.ui.util.MessageUtils;
import org.igv.ui.util.SortDialog;

import javax.swing.*;
import java.util.List;

public class SampleMenuUtils {

    public static JMenuItem getSortByAttributeItem(AbstractTrack track) {

        JMenuItem item = new JMenuItem("Sort Samples By Attribute...");
        item.addActionListener(evt -> {

            List<String> keys = AttributeManager.getInstance().getAttributeNames();
            Object availableSortKeys[] = keys.toArray();
            SortDialog dialog = new SortDialog(IGV.getInstance().getMainFrame(), true, availableSortKeys);
            dialog.setVisible(true);

            if (dialog.isCanceled()) {
                return;
            }

            final String[] attributeNames = dialog.getSelectedSortKeys();
            if (attributeNames != null) {
                final boolean[] ascending = dialog.isAscending();
                AttributeComparator.SampleAttributeComparator comparator = new AttributeComparator.SampleAttributeComparator(attributeNames, ascending);
                track.sortSamples(comparator);
                track.repaint();
            }
        });

        return item;
    }


    public static JMenuItem getGroupByAttributeItem(AbstractTrack track) {

        JMenuItem item = new JMenuItem("Group Samples By Attribute...");

        item.addActionListener(evt -> {

            final AttributeSelectionDialog dlg = new AttributeSelectionDialog(
                    IGV.getInstance().getMainFrame(),
                    "Group");


            String currentSelection = IGV.getInstance().getGroupByAttribute();
            if (currentSelection == null) {
                dlg.setSelectedIndex(0);
            } else {
                dlg.setSelectedItem(currentSelection);
            }

            dlg.setVisible(true);

            if (!dlg.isCanceled()) {
                String selectedAttribute = dlg.getSelected();
                track.setSampleGroupBy(selectedAttribute);
            }
        });

        return item;
    }


    public static JMenuItem getFilterByAttributeItem(AbstractTrack track) {

        JMenuItem item = new JMenuItem("Filter Samples By Attribute...");

        item.addActionListener(evt -> {

            List<String> uniqueAttributeKeys = AttributeManager.getInstance().getAttributeNames();

            // Sort the attribute keys if we have any
            if (uniqueAttributeKeys != null) {
                //Collections.sort(uniqueAttributeKeys, AttributeManager.getInstance().getAttributeComparator());
            } else // If we have no attribute we can't display the
                // track filter dialog so say so and return
                if (uniqueAttributeKeys == null || uniqueAttributeKeys.isEmpty()) {

                    MessageUtils.showMessage("No attributes found to use in a filter");
                    return;
                }

            SampleFilter sampleFilter = IGV.getInstance().getSession().getFilter();
            SampleFilterDialog dialog = new SampleFilterDialog(IGV.getInstance().getMainFrame(), "Filter Tracks", sampleFilter);
            dialog.setVisible(true);

            if (!dialog.isCancelled()) {
                sampleFilter = dialog.getFilter();
                track.setSampleFilter(sampleFilter);

            }
        });

        return item;
    }

}
