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

package org.broad.igv.aws;

import htsjdk.samtools.util.Tuple;
import org.apache.log4j.Logger;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.ga4gh.OAuthUtils;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.commandbar.GenomeListManager;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.AmazonUtils;
import org.broad.igv.util.LongRunningTask;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.TreeExpansionEvent;
import javax.swing.event.TreeWillExpandListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.ExpandVetoException;
import javax.swing.tree.TreePath;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;


public class S3LoadDialog extends JDialog {

    private static Logger log = Logger.getLogger(S3LoadDialog.class);

    private final DefaultTreeModel treeModel;
    String selectedId;

    public S3LoadDialog(Frame owner) {
        super(owner);
        ArrayList<String> datasets = AmazonUtils.ListBucketsForUser();
        initComponents();

        S3TreeNode root = new S3TreeNode(new S3Object("S3", true), true);

        treeModel = new DefaultTreeModel(root);
        this.selectionTree.setModel(treeModel);

        // List toplevel buckets
        ArrayList<String> buckets = AmazonUtils.ListBucketsForUser();
        for (String bucket: buckets) {
            S3Object bucket_obj = new S3Object(bucket, true);
            root.add(new S3TreeNode(bucket_obj, true));
        }

        log.debug("Populated S3 load dialog with S3 buckets: "+ buckets.toString());

        // Propagate changes to UI
        treeModel.reload();
    }

    private void loadButtonActionPerformed(ActionEvent e) {

        setVisible(false);

        LongRunningTask.submit(() -> {
            TreePath[] paths = selectionTree.getSelectionPaths();
            ArrayList<Tuple<String, String>> preLocatorPaths = new ArrayList<>();
            ArrayList<ResourceLocator> finalLocators = new ArrayList<>();

            for (TreePath path : paths) {
                DefaultMutableTreeNode obj = (DefaultMutableTreeNode) path.getLastPathComponent();
                Object[] selectedObjects = path.getPath();

                String bucketName = ((S3TreeNode) selectedObjects[1]).getUserObject().getName();
                String s3Key = "";

                for (int i = 2; i < selectedObjects.length; i++) {
                    S3TreeNode selectedObject = (S3TreeNode) selectedObjects[i];
                    S3Object selectedS3Object = selectedObject.getUserObject();
                    s3Key += selectedS3Object.getName() + "/";
                }

                s3Key = s3Key.substring(0, s3Key.length() - 1);
                log.debug("Loading S3 object key: " + s3Key + " from bucket " + bucketName);

                preLocatorPaths.add(new Tuple<>(bucketName, s3Key));
            }

            for (Tuple<String, String> preLocator: preLocatorPaths) {
                String bucketName = preLocator.a;
                String s3objPath = preLocator.b;

                // XXX: Trace back path where this happens on other flows: both main file and index as presigned urls
                URL s3_presigned_url = AmazonUtils.translateAmazonCloudURL( bucketName,
                                                                            s3objPath,
                                                                            new java.util.Date(OAuthUtils.getExpirationTime()));

                // XXX: Make it general for all non-bam.bai.
                URL s3_presigned_url_idx = AmazonUtils.translateAmazonCloudURL( bucketName,
                                                                                s3objPath.replace(".bam", ".bam.bai"),
                                                                                new java.util.Date(OAuthUtils.getExpirationTime()));

                ResourceLocator locator = new ResourceLocator(s3_presigned_url.toString());
                locator.setName("");
                locator.setType(".bam"); // XXX: Trace back where this is auto-detected
                locator.setIndexPath(s3_presigned_url_idx.toString());

                //locator.setType(Ga4ghAPIHelper.RESOURCE_TYPE);
                //locator.setAttribute("provider", Ga4ghAPIHelper.GA4GH_GOOGLE_PROVIDER);

                finalLocators.add(locator);

            }

            IGV.getInstance().loadTracks(finalLocators);
        });
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        selectedId = null;
        setVisible(false);
    }

    class LeafNode {
        ResourceLocator locator;
        ArrayList<String> datasets;
        String dataset;

        public LeafNode(ResourceLocator locator) {
            this.locator = locator;
        }

        public LeafNode(ArrayList<String> datasets) {
            this.datasets = datasets;
        }

        public LeafNode(String dataset){
            this.dataset = dataset;
        }

        @Override
        public String toString() {
            return dataset;
        }
    }

    private void loadTrack(String url, String name) {
        ResourceLocator locator = new ResourceLocator(url);
        locator.setName(name);
        IGV.getInstance().loadTracks(Arrays.asList(locator));

    }

    private void setGenome(String genomeId) {

        if (genomeId != null && !genomeId.equals(GenomeManager.getInstance().getGenomeId())) {
            GenomeListItem item = GenomeListManager.getInstance().getGenomeListItem(genomeId);
            if (item != null) {
                try {
                    GenomeManager.getInstance().loadGenomeById(genomeId);
                } catch (IOException e) {
                    MessageUtils.showErrorMessage("Error loading genome: " + genomeId, e);
                    log.error("Error loading genome: " + genomeId, e);
                }
            }
        }
    }

    private void updateModel(DefaultMutableTreeNode parent) {
        DefaultTreeModel model = treeModel;
        if (parent != null) {
            model.reload(parent);
        } else {
            model.reload();
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
        setTitle("Amazon S3 datasets");
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== selectionTree ========
        selectionTree.addTreeWillExpandListener(new TreeWillExpandListener() {
              public void treeWillExpand(TreeExpansionEvent event) throws ExpandVetoException {
                  log.debug("TreeWillExpand on S3LoadDialog panel");

                  Object parent = event.getPath().getLastPathComponent();
                  S3TreeNode parentNode = (S3TreeNode) parent;
                  S3Object s3Object = parentNode.getUserObject();
                  //String s3Bucket = s3Object.toString();

                  if (s3Object.isDir()) {
                      Object[] path = parentNode.getUserObjectPath(); // fullpath to S3 object
                      String currentBucket = path[1].toString();
                      String prefix = "";

                      for (int i = 2; i < path.length; i++) {
                          prefix += path[i] + "/";
                      }

                      log.debug("S3 bucket prefix is: "+prefix);

                      // List contents of bucket with path-prefix passed
                      // XXX: Determine reliable source for bucket root instead of hardcoding it here
                      ArrayList<S3Object> s3Objects = AmazonUtils.ListBucketObjects(currentBucket, prefix);

                      // For each item in the bucket:
                      //  1) create an S3Object
                      //    1.1) Name of object.
                      //    1.2) Dir or file: .getPrefix or null, according to S3 API
                      for (S3Object s3Obj: s3Objects) {
                          // Add it to the corresponding POJO...
                          // XXX: Cannot if (s3Obj.getPrefix) here... find string alternatives for PRE detection?:
                          // https://docs.aws.amazon.com/AmazonS3/latest/dev/UsingMetadata.html
                          // """
                          // The Amazon S3 data model is a flat structure: you create a bucket, and the bucket stores
                          // objects. There is no hierarchy of subbuckets or subfolders; however, you can infer logical
                          // hierarchy using key name prefixes and delimiters as the Amazon S3 console does. The Amazon
                          // S3 console supports a concept of folders.
                          // """
                          parentNode.add(new S3TreeNode(s3Obj));
                      }

                      // ... and update the model
                      updateModel(parentNode);
                  }
              }

              public void treeWillCollapse(TreeExpansionEvent event) throws ExpandVetoException {
                  log.debug("TreeWillCollapse on S3LoadDialog panel");
              }
        });

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
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

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
            label1.setText("Select a objects to load");
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
