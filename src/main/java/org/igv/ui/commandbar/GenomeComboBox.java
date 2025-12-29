package org.igv.ui.commandbar;

import org.igv.feature.genome.GenomeManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.ui.UIConstants;
import org.igv.ui.genome.GenomeListItem;
import org.igv.ui.genome.GenomeListManager;
import org.igv.ui.util.MessageUtils;
import org.igv.util.LongRunningTask;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.Vector;

/**
 * Created by jrobinso on 7/6/17.
 */
public class GenomeComboBox extends JComboBox<GenomeListItem> {

    private static Logger log = LogManager.getLogger(GenomeComboBox.class);

    public GenomeComboBox() {

        setRenderer(new ComboBoxRenderer());
        addActionListener(new GenomeBoxActionListener());
    }


    public void refreshGenomeListComboBox() {

        setModel(buildModel());

        String curId = GenomeManager.getInstance().getGenomeId();
        if (curId == null) return;

        int c = this.getItemCount();
        for (int i = 0; i < c; i++) {
            final GenomeListItem item = this.getItemAt(i);
            if (curId.equals(item.getId())) {
                setSelectedItem(item);
                break;
            }
        }
    }

    /**
     * Build a model for the genome combo box
     *
     * @return
     */
    private DefaultComboBoxModel buildModel() {
        Collection<GenomeListItem> genomes;
        try {
            genomes = GenomeListManager.getInstance().getGenomeItemMap().values();
        } catch (IOException e) {
            log.error("Error reading genome list ", e);
            genomes = new ArrayList<>();
            MessageUtils.showErrorMessage("Error reading genome list ", e);
        }

        Vector<GenomeListItem> vector = new Vector<>(genomes);
        vector.sort(Comparator.comparing(GenomeListItem::getDisplayableName));
        return new DefaultComboBoxModel(vector);
    }

    class GenomeBoxActionListener implements ActionListener {

        public void actionPerformed(ActionEvent actionEvent) {
            Object selItem = getSelectedItem();
            if (!(selItem instanceof GenomeListItem)) {
                return;
            }
            GenomeListItem genomeListItem = (GenomeListItem) selItem;

            // If we haven't changed genomes do nothing
            if (genomeListItem.getId().equalsIgnoreCase(GenomeManager.getInstance().getGenomeId())) {
                return;
            }

            final Runnable runnable = () -> {

                if (genomeListItem != null && genomeListItem.getPath() != null) {
                    try {
                        GenomeManager.getInstance().loadGenome(genomeListItem.getPath());
                    } catch (Exception e) {
                        log.error(e);
                        MessageUtils.showErrorMessage("The genome '" + genomeListItem.getDisplayableName() +
                                "' could not be read.", e);
                    }
                }
            };

            // If we're on the dispatch thread spawn a worker, otherwise just execute.
            if (SwingUtilities.isEventDispatchThread()) {
                LongRunningTask.submit(runnable);
            } else {
                runnable.run();
            }
        }
    }


    static class ComboBoxRenderer implements ListCellRenderer {

        JSeparator separator;

        /**
         * Constructs ...
         */
        public ComboBoxRenderer() {
            separator = new JSeparator(JSeparator.HORIZONTAL);
        }

        /**
         * Method description
         *
         * @param list
         * @param value
         * @param index
         * @param isSelected
         * @param cellHasFocus
         * @return
         */
        public Component getListCellRendererComponent(JList list, Object value, int index,
                                                      boolean isSelected, boolean cellHasFocus) {
            String text = (value == null) ? "" : value.toString();

            Component renderer = null;

            if (UIConstants.GENOME_LIST_SEPARATOR.equals(text)) {
                return separator;
            }

            if (text.equals(UIConstants.REMOVE_GENOME_LIST_MENU_ITEM)) {
                JLabel label = new JLabel(text);

                label.setOpaque(true);
                label.setBorder(new EmptyBorder(1, 1, 1, 1));
                renderer = label;
            } else {

                JLabel label = new JLabel(text);

                label.setOpaque(true);
                label.setBorder(new EmptyBorder(1, 1, 1, 1));
                label.setSize(label.getWidth() + 10, label.getHeight());
                renderer = label;
            }

            //We call with a null list when setting width
            if (list != null) {
                if (isSelected) {
                    renderer.setBackground(list.getSelectionBackground());
                    renderer.setForeground(list.getSelectionForeground());
                } else {
                    renderer.setBackground(list.getBackground());
                    renderer.setForeground(list.getForeground());
                }
                renderer.setFont(list.getFont());
            }


            return renderer;
        }
    }


}
