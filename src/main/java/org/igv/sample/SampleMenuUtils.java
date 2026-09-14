package org.igv.sample;

import org.igv.track.AbstractTrack;
import org.igv.track.AttributeManager;
import org.igv.ui.AttributeSelectionDialog;
import org.igv.ui.IGV;
import org.igv.ui.SampleFilterDialog;
import org.igv.ui.SampleSelectionDialog;
import org.igv.ui.util.MessageUtils;
import org.igv.ui.util.SortDialog;

import javax.swing.*;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

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
                track.sortSamplesByAttributes(attributeNames, ascending);
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


            String currentSelection = track.getGroupBy();
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

            SampleFilter sampleFilter = track.getSampleFilter();
            SampleFilterDialog dialog = new SampleFilterDialog(IGV.getInstance().getMainFrame(), "Filter Tracks", sampleFilter);
            dialog.setVisible(true);

            if (!dialog.isCancelled()) {
                sampleFilter = dialog.getFilter();
                track.setSampleFilter(sampleFilter);

            }
        });

        return item;
    }


    public static JMenuItem getFilterByIdItem(AbstractTrack track) {

        JMenuItem item = new JMenuItem("Filter Samples By ID...");

        item.addActionListener(evt -> {

            // Start with the current ID filter, or all samples.  The attribute filter is separate and not reflected here.
            List<String> currentSamples = track.getSelectedSamples() == null ? track.getSampleNames() : track.getSelectedSamples();
            SampleSelectionDialog dialog = new SampleSelectionDialog(IGV.getInstance().getMainFrame(), currentSamples, track.getSampleNames());
            dialog.setVisible(true);

            if (dialog.isCanceled()) {
                return;
            }

            List<String> ids = dialog.getSampleIds();

            // OK without editing the list is the same as Cancel.  Reordering a partial list is an edit, as it sets the
            // display order of an unsorted track.
            if ((ids == null ? List.of() : ids).equals(currentSamples)) {
                return;
            }

            if (ids == null) {
                track.setSelectedSamples(null);
                return;
            }

            Set<String> trackSamples = new HashSet<>(track.getSampleNames());
            List<String> found = ids.stream().filter(trackSamples::contains).collect(Collectors.toList());
            List<String> notFound = ids.stream().filter(id -> !trackSamples.contains(id)).collect(Collectors.toList());

            if (found.isEmpty()) {
                MessageUtils.showMessage("None of the entered sample IDs were found in this track.");
                return;
            }

            track.setSelectedSamples(found);

            if (!notFound.isEmpty()) {
                MessageUtils.showMessage(notFound.size() + " of " + ids.size() + " sample IDs were not found: " +
                        String.join(", ", notFound.subList(0, Math.min(20, notFound.size()))) +
                        (notFound.size() > 20 ? ", ..." : ""));
            }
        });

        return item;
    }

}
