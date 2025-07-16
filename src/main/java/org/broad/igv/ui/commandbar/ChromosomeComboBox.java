package org.broad.igv.ui.commandbar;

import org.broad.igv.Globals;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by jrobinso on 7/6/17.
 */
public class ChromosomeComboBox extends JComboBox {


    private static int MAX_CHROMOSOME_COUNT = 10000;
    private static String MAX_EXCEEDED = "Max exceeded";

    public ChromosomeComboBox() {
        addActionListener(evt -> chromosomeComboBoxActionPerformed(evt));
    }

    private void chromosomeComboBoxActionPerformed(java.awt.event.ActionEvent evt) {
        JComboBox combobox = (JComboBox) evt.getSource();
        final String chrDisplayName = (String) combobox.getSelectedItem();
        final String chrName = GenomeManager.getInstance().getCurrentGenome().getCanonicalChrName(chrDisplayName);
        if (chrName != null & !chrName.equals(FrameManager.getDefaultFrame().getChrName()) && !chrName.equals(MAX_EXCEEDED)) {
            FrameManager.getDefaultFrame().changeChromosome(chrName, true);
        }
    }


    public void updateChromosFromGenome(Genome genome) {

        if (genome == null) return;

        UIUtilities.invokeAndWaitOnEventThread(() -> {

            final List<String> chromosomeNames = genome.getChromosomeNames();
            if (chromosomeNames == null) {
                this.setVisible(false);
            } else {
                List<String> allChromosomeNames = chromosomeNames.stream().map(chr -> genome.getChromosomeDisplayName(chr)).collect(Collectors.toList());

                if (allChromosomeNames.size() > 1) {
                    this.setVisible(true);

                    List<String> tmp = new ArrayList<>();
                    tmp.addAll(allChromosomeNames.size() > MAX_CHROMOSOME_COUNT ?
                            allChromosomeNames.subList(0, MAX_CHROMOSOME_COUNT) :
                            allChromosomeNames);

                    String homeChr = genome.getHomeChromosome();
                    if (homeChr.equals(Globals.CHR_ALL)) {
                        tmp.add(0, Globals.CHR_ALL);
                    }

                    if (allChromosomeNames.size() > MAX_CHROMOSOME_COUNT) {
                        tmp.add(MAX_EXCEEDED);
                    }

                    Graphics2D graphics2D = (Graphics2D) getGraphics();
                    Font font = getFont();
                    FontMetrics fontMetrics = getFontMetrics(font);

                    int w = IGVCommandBar.DEFAULT_CHROMOSOME_DROPDOWN_WIDTH;
                    for (String chromosomeName : tmp) {
                        Rectangle2D textBounds = fontMetrics.getStringBounds(chromosomeName, graphics2D);
                        if (textBounds != null) {
                            int width = textBounds.getBounds().width + 50;
                            if (width > w) {
                                w = width;
                            }
                        }
                    }

                    final DefaultComboBoxModel defaultModel = new DefaultComboBoxModel(tmp.toArray());
                    final int dropdownWidth = w;

                    setModel(defaultModel);
                    setSelectedItem(genome.getHomeChromosome());
                    adjustChromosomeDropdownWidth(dropdownWidth);
                } else {
                    this.setVisible(false);
                }
            }
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
