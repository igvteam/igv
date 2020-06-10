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

import org.apache.log4j.Logger;
import org.apache.commons.lang3.tuple.Triple;
import org.broad.igv.ui.IGV;
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
import javax.swing.tree.TreePath;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import software.amazon.awssdk.services.s3.model.S3Exception;

import static org.broad.igv.util.AmazonUtils.isObjectAccessible;


public class S3LoadDialog extends JDialog {

    private static Logger log = Logger.getLogger(S3LoadDialog.class);

    private final DefaultTreeModel treeModel;
    String selectedId;

    public S3LoadDialog(Frame owner) {
        super(owner);
        initComponents();

        S3TreeNode root = new S3TreeNode(new IGVS3Object("S3", true, "STANDARD"), true);

        treeModel = new DefaultTreeModel(root);
        this.selectionTree.setModel(treeModel);

        // List toplevel buckets
        List<String> buckets = AmazonUtils.ListBucketsForUser();
        for (String bucket: buckets) {
            IGVS3Object bucket_obj = new IGVS3Object(bucket, true, "STANDARD");
            root.add(new S3TreeNode(bucket_obj, true));
        }

        log.debug("Populated S3 load dialog with S3 buckets: "+ buckets.toString());

        // Propagate changes to UI
        treeModel.reload();
    }

    /**
     * Loads the selected resources when the Load button is pressed.
     * This one only one way of loading files. The other one is via mouse double click (see initComponents())
     * The main differences are that this method can handle multiple selected files, which are defined by the selection
     * in the tree at the time of loading, vs the mouse coordinates in when double clicking.
     * @param e
     */
    private void loadButtonActionPerformed(ActionEvent e) {

        setVisible(false);

        LongRunningTask.submit(() -> {
            TreePath[] paths = selectionTree.getSelectionPaths();
            // Reminder: ArrayList<Triple<Bucket, Key, StorageClass>>
            ArrayList<Triple<String, String, String>> preLocatorPaths = new ArrayList<>();
            ArrayList<ResourceLocator> finalLocators = new ArrayList<>();

            for (TreePath path : paths) {
                if (isFilePath(path)) {
                    Triple<String, String, String> bucketKeyTier = getBucketKeyTierFromTreePath(path);

                    AmazonUtils.s3ObjectAccessResult res = isObjectAccessible(bucketKeyTier.getLeft(), bucketKeyTier.getMiddle());
                    if(!res.getObjAvailable()) { MessageUtils.showErrorMessage(res.getErrorReason(), null); return; }

                    preLocatorPaths.add(bucketKeyTier);
                }
            }

            for (Triple<String, String, String> preLocator: preLocatorPaths) {
                ResourceLocator locator = getResourceLocatorFromBucketKey(preLocator);
                finalLocators.add(locator);
            }

            IGV.getInstance().loadTracks(finalLocators);
        });
    }

    private boolean isFilePath(TreePath path) {
        return ((S3TreeNode) path.getLastPathComponent()).isLeaf();
    }

    private ResourceLocator getResourceLocatorFromBucketKey(Triple<String, String, String> preLocator) {
        String bucketName = preLocator.getLeft();
        String s3objPath = preLocator.getMiddle();

        String s3Path = "s3://"+bucketName+"/"+s3objPath;

        return new ResourceLocator(s3Path);
    }

    private Triple<String, String, String> getBucketKeyTierFromTreePath(TreePath path) {
        Object[] selectedObjects = path.getPath();

        String bucketName = ((S3TreeNode) selectedObjects[1]).getUserObject().getName();
        String s3Key = "";
        String storageClass = "STANDARD";

        for (int i = 2; i < selectedObjects.length; i++) {
            S3TreeNode selectedObject = (S3TreeNode) selectedObjects[i];
            IGVS3Object selectedIGVS3Object = selectedObject.getUserObject();
            storageClass = selectedIGVS3Object.getStorageClass();
            s3Key += selectedIGVS3Object.getName() + "/";
        }

        s3Key = s3Key.substring(0, s3Key.length() - 1);
        log.debug("Loading S3 object key: " + s3Key + " from bucket " + bucketName);

        return Triple.of(bucketName, s3Key, storageClass);
    }

    private void treeWillExpandActionPerformed(TreeExpansionEvent event) {
        log.debug("TreeWillExpand on S3LoadDialog panel");

        Object parent = event.getPath().getLastPathComponent();
        S3TreeNode parentNode = (S3TreeNode) parent;
        IGVS3Object IGVS3Object = parentNode.getUserObject();

        // only load child folder/files if they haven't been loaded already
        if (parentNode.getChildCount() == 0) {
            if (IGVS3Object.isDir()) {
                Object[] path = parentNode.getUserObjectPath(); // fullpath to S3 object
                String currentBucket = path[1].toString();
                String prefix = "";

                for (int i = 2; i < path.length; i++) {
                    prefix += path[i] + "/";
                }

                log.debug("S3 bucket prefix is: "+prefix);

                try {
                    // List contents of bucket with path-prefix passed
                    ArrayList<IGVS3Object> IGVS3Objects = AmazonUtils.ListBucketObjects(currentBucket, prefix);
                    parentNode.addS3Children(IGVS3Objects);
                } catch (S3Exception e) {
                    MessageUtils.showErrorMessage("Amazon S3: Access denied to bucket: "+currentBucket, e);
                    log.error("Permission denied on S3 bucket ListObjects: ");
                }

                // ... and update the model
                updateModel(parentNode);
            }
        }
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        selectedId = null;
        setVisible(false);
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
            public void treeWillExpand(TreeExpansionEvent event) {
                treeWillExpandActionPerformed(event);
            }

            public void treeWillCollapse(TreeExpansionEvent event) {
                log.debug("Tree collapsing");
            }
        });
        selectionTree.addMouseListener(new MouseAdapter() {
            public void mousePressed(MouseEvent e) {
                int selRow = selectionTree.getRowForLocation(e.getX(), e.getY());
                TreePath selPath = selectionTree.getPathForLocation(e.getX(), e.getY());
                if(selRow != -1) {
                    if(e.getClickCount() == 2) {
                        // similar behaviour to loadButtonActionPerformed, see docs there for details
                        if (isFilePath(selPath)) {
                            Triple<String, String, String> bucketKeyTier = getBucketKeyTierFromTreePath(selPath);

                            AmazonUtils.s3ObjectAccessResult res = isObjectAccessible(bucketKeyTier.getLeft(), bucketKeyTier.getMiddle());
                            if(!res.getObjAvailable()) { MessageUtils.showErrorMessage(res.getErrorReason(), null); return;}

                            ResourceLocator loc = getResourceLocatorFromBucketKey(bucketKeyTier);
                            IGV.getInstance().loadTracks(Collections.singletonList(loc));
                            setVisible(false);
                        }
                    }
                }
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
                loadButton.addActionListener(e -> loadButtonActionPerformed(e));
                buttonBar.add(loadButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 5), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(e -> cancelButtonActionPerformed(e));
                buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                        GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                        new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);

            //---- label1 ----
            label1.setText("Select objects to load");
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
