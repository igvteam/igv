package org.broad.igv.ui.commandbar;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by jrobinso on 7/6/17.
 */
public class ChromosomeComboBox extends JComboBox {


    public ChromosomeComboBox() {
        addActionListener(evt -> chromosomeComboBoxActionPerformed(evt));
    }

    private void chromosomeComboBoxActionPerformed(java.awt.event.ActionEvent evt) {
        JComboBox combobox = (JComboBox) evt.getSource();
        final String chrName = (String) combobox.getSelectedItem();
        if (chrName != null & !chrName.equals(FrameManager.getDefaultFrame().getChrName())) {
            FrameManager.getDefaultFrame().changeChromosome(chrName, true);
        }
    }


    public void updateChromosFromGenome(Genome genome) {

        if (genome == null) return;

        UIUtilities.invokeAndWaitOnEventThread(() -> {

            List<String> tmp = new ArrayList<String>(genome.getAllChromosomeNames().size());
            tmp.addAll(genome.getAllChromosomeNames());
            if (tmp.size() > 1) {
                String homeChr = genome.getHomeChromosome();
                if (homeChr.equals(Globals.CHR_ALL)) {
                    tmp.add(0, Globals.CHR_ALL);
                }
            }

            Graphics2D graphics2D = (Graphics2D) getGraphics();
            Font font = getFont();
            FontMetrics fontMetrics = getFontMetrics(font);

            int w = IGVCommandBar.DEFAULT_CHROMOSOME_DROPDOWN_WIDTH;
            for (String chromosomeName : tmp) {
                Rectangle2D textBounds = fontMetrics.getStringBounds(chromosomeName, graphics2D);
                if (textBounds != null) {
                    int width = textBounds.getBounds().width + 50;

                    // int width = chromosomeName.length()*fontSize-(fontSize*4);  // TODO Hack figure out whats's wrong with previous line
                    if (width > w) {
                        w = width;
                    }
                }
            }

            Object[] chomosomeNames = tmp.toArray();
            final DefaultComboBoxModel defaultModel = new DefaultComboBoxModel(chomosomeNames);
            final int dropdownWidth = w;

            setModel(defaultModel);
            setSelectedItem(genome.getHomeChromosome());
            adjustChromosomeDropdownWidth(dropdownWidth);
        });
    }


    private void adjustChromosomeDropdownWidth(int width) {

        int newWidth = (width > IGVCommandBar.DEFAULT_CHROMOSOME_DROPDOWN_WIDTH)
                ? width : IGVCommandBar.DEFAULT_CHROMOSOME_DROPDOWN_WIDTH;

        setMaximumSize(new java.awt.Dimension(newWidth, 35));
        setMinimumSize(new java.awt.Dimension(newWidth, 27));
        setPreferredSize(new java.awt.Dimension(newWidth, 16));
        revalidate();
    }

}
