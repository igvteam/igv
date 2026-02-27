package org.igv.sample;

import org.igv.track.AttributeManager;
import org.igv.track.Track;
import org.igv.ui.AttributeSelectionDialog;
import org.igv.ui.IGV;
import org.igv.ui.util.SortDialog;

import javax.swing.*;
import java.util.List;

public class SampleMenuUtils {

    public static JMenuItem getSortByAttributeItem(Track track) {

        JMenuItem item = new JMenuItem("Sort By Attribute...");
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


    public static JMenuItem getGroupByAttributeItem(Track track) {

        JMenuItem item = new JMenuItem("Group By Attribute...");

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
                track.groupSamplesByAttribute(selectedAttribute);
            }
        });

        return item;
    }
}
