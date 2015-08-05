/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/*
 * Created by JFormDesigner on Wed Apr 27 19:46:06 EDT 2011
 */

package org.broad.igv.lists;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.util.FileDialogUtils;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.NamedRunnable;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;

/**
 * @author Jim Robinson
 */
public class VariantListNavigator extends JDialog {

    public VariantListNavigator(Frame owner) {
        super(owner);
        initComponents();
        setAlwaysOnTop(true);
        this.setLocation(owner.getBounds().width - 200, 50);
        refreshList();
    }

    private void refreshList() {
        VariantListEntry[] variantListArray = VariantListManager.variantList.toArray(new VariantListEntry[]{});
        variantList.setListData(variantListArray);

    }

    private void variantListValueChanged(ListSelectionEvent e) {
        processVariant((VariantListEntry) variantList.getSelectedValue());

    }

    private void loadVariantsButtonActionPerformed(ActionEvent e) {

        File variantListFile = FileDialogUtils.chooseFile("Import variant list");
        if (variantListFile != null) {
            try {
                VariantListManager.loadVariants(new ResourceLocator(variantListFile.getAbsolutePath()));
            } catch (IOException e1) {
                MessageUtils.showMessage("Error loading variant file: " + e.toString());
            }
        }
        refreshList();

    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        scrollPane1 = new JScrollPane();
        variantList = new JList();
        panel1 = new JPanel();
        loadVariantsButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== scrollPane1 ========
        {

            //---- variantList ----
            variantList.addListSelectionListener(new ListSelectionListener() {
                public void valueChanged(ListSelectionEvent e) {
                    variantListValueChanged(e);
                }
            });
            scrollPane1.setViewportView(variantList);
        }
        contentPane.add(scrollPane1, BorderLayout.CENTER);

        //======== panel1 ========
        {
            panel1.setLayout(null);

            //---- loadVariantsButton ----
            loadVariantsButton.setText("Load variants...");
            loadVariantsButton.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    loadVariantsButtonActionPerformed(e);
                }
            });
            panel1.add(loadVariantsButton);
            loadVariantsButton.setBounds(new Rectangle(new Point(0, 0), loadVariantsButton.getPreferredSize()));

            { // compute preferred size
                Dimension preferredSize = new Dimension();
                for (int i = 0; i < panel1.getComponentCount(); i++) {
                    Rectangle bounds = panel1.getComponent(i).getBounds();
                    preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                    preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                }
                Insets insets = panel1.getInsets();
                preferredSize.width += insets.right;
                preferredSize.height += insets.bottom;
                panel1.setMinimumSize(preferredSize);
                panel1.setPreferredSize(preferredSize);
            }
        }
        contentPane.add(panel1, BorderLayout.NORTH);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }


    public static void main(String[] args) {
        (new VariantListNavigator(null)).setVisible(true);

    }

    void processVariant(final VariantListEntry variant) {

        VariantListNavigator.this.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        VariantListNavigator.this.variantList.setEnabled(false);

        NamedRunnable runnable = new NamedRunnable() {

            public String getName() {
                return "Process variant: " + variant.toString();
            }

            public void run() {

                try {

                    // Lookup files from samples  TODO accept paths here
                    java.util.List<ResourceLocator> locators = new ArrayList();
                    for (String sample : variant.samples) {
                        String path = VariantListManager.samplePathMap.get(sample);
                        if (path == null) {
                            MessageUtils.showMessage("Sample not found: " + sample);
                        } else {
                            locators.add(new ResourceLocator(path));
                        }
                    }

                    // New session
                    IGV.getInstance().resetSession(null);
                    IGV.getInstance().loadTracks(locators);

                }
                finally {
                    Runnable runnable = new Runnable() {
                        public void run() {
                            Genome genome = GenomeManager.getInstance().getCurrentGenome();
                            String chr = genome == null ? variant.chr : genome.getChromosomeAlias(variant.chr);
                            String locus = chr + ":" + variant.position;
                            IGV.getInstance().goToLocus(locus);
                            VariantListNavigator.this.setCursor(Cursor.getDefaultCursor());
                            VariantListNavigator.this.variantList.setEnabled(true);
                        }
                    };
                    SwingUtilities.invokeLater(runnable);
                }
            }
        };

        LongRunningTask.submit(runnable);

    }


    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JScrollPane scrollPane1;
    private JList variantList;
    private JPanel panel1;
    private JButton loadVariantsButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
