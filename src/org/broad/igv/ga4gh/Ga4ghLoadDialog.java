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

package org.broad.igv.ga4gh;

import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ui.IGV;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ResourceLocator;

import java.awt.*;
import java.awt.event.*;
import java.io.IOException;
import java.util.Arrays;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;

/**
 * @author James Robinson
 */
public class Ga4ghLoadDialog extends JDialog {

    private final DefaultTreeModel treeModel;
    Ga4ghProvider[] providers;
    String selectedId;

    public Ga4ghLoadDialog(Frame owner, Ga4ghProvider[] providers) {

        super(owner);

        initComponents();

        this.providers = providers;

        treeModel = new DefaultTreeModel(createNodes(providers));

        this.selectionTree.setModel(treeModel);
    }

    private DefaultMutableTreeNode createNodes(Ga4ghProvider[] providers) {

        DefaultMutableTreeNode top = new DefaultMutableTreeNode("Ga4gh");

        for (Ga4ghProvider provider : providers) {

            DefaultMutableTreeNode providerNode = new DefaultMutableTreeNode(provider.getName());
            top.add(providerNode);

            for (Ga4ghDataset dataset : provider.getDatasets()) {

                DefaultMutableTreeNode datasetNode = new DefaultMutableTreeNode(dataset.getName());
                providerNode.add(datasetNode);

                for (Ga4ghReadset readset : dataset.getReadsets()) {

                    DefaultMutableTreeNode readsetNode = new DefaultMutableTreeNode(new LeafNode(provider, readset) );
                    datasetNode.add(readsetNode);

                }

            }
        }
        return top;
    }

    private void loadButtonActionPerformed(ActionEvent e) {

        setVisible(false);

        LongRunningTask.submit(new Runnable() {
            public void run() {
                TreePath[] paths = selectionTree.getSelectionPaths();

                for (TreePath path : paths) {
                    DefaultMutableTreeNode obj = (DefaultMutableTreeNode) path.getLastPathComponent();
                    Object userObject = obj.getUserObject();
                    if (userObject instanceof LeafNode) {
                            Ga4ghProvider provider = ((LeafNode) userObject).provider;
                            Ga4ghReadset readSet = ((LeafNode) userObject).readset;
                            setGenome(readSet.getGenomeId());
                            loadTrack(readSet.getId(), provider, readSet.getName());

                    }
                }
            }
        });
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        selectedId = null;
        setVisible(false);
    }

    class LeafNode {
        Ga4ghProvider provider;
        Ga4ghReadset readset;

        LeafNode(Ga4ghProvider provider, Ga4ghReadset readset) {
            this.provider = provider;
            this.readset = readset;
        }

        public String toString() {return readset.getName();}
    }

    private void loadTrack(String readsetId, Ga4ghProvider provider, String name) {

        ResourceLocator locator = new ResourceLocator(readsetId);
        locator.setName(name);
        locator.setType(Ga4ghAPIHelper.RESOURCE_TYPE);
        locator.setAttribute("provider", provider);
        IGV.getInstance().loadTracks(Arrays.asList(locator));

    }

    private void setGenome(String genomeId) {

        if (genomeId != null && !genomeId.equals(GenomeManager.getInstance().getGenomeId())) {
            try {
                GenomeListItem item = GenomeManager.getInstance().findGenomeListItemById(genomeId);
                if (item != null) {
                    IGV.getInstance().loadGenomeById(genomeId);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        scrollPane1 = new JScrollPane();
        selectionTree = new JTree();
        buttonBar = new JPanel();
        loadButton = new JButton();
        cancelButton = new JButton();
        label1 = new JLabel();

        //======== this ========
        setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("GA4GH Prototype");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.X_AXIS));

                //======== scrollPane1 ========
                {
                    scrollPane1.setViewportView(selectionTree);
                }
                contentPanel.add(scrollPane1);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

                //---- loadButton ----
                loadButton.setText("Load");
                loadButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        loadButtonActionPerformed(e);
                    }
                });
                buttonBar.add(loadButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 5), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);

            //---- label1 ----
            label1.setText("Select a readset to load");
            dialogPane.add(label1, BorderLayout.NORTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        setSize(795, 690);
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JScrollPane scrollPane1;
    private JTree selectionTree;
    private JPanel buttonBar;
    private JButton loadButton;
    private JButton cancelButton;
    private JLabel label1;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
