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

package org.broad.igv.ui;

import org.apache.log4j.Logger;
import org.broad.igv.track.AttributeManager;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.ui.util.LinkCheckBox;
import org.broad.igv.util.ResourceLocator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.swing.*;
import javax.swing.tree.*;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.*;
import java.util.List;

import static org.broad.igv.util.ResourceLocator.AttributeType.*;

/**
 * Parses XML file of IGV resources, and displays them in tree format.
 *
 * @author eflakes
 */
public class ResourceTree {

    private StringBuffer buffer = new StringBuffer();
    private static String FAILED_TO_CREATE_RESOURCE_TREE_DIALOG = "Failure while creating the resource tree dialog";
    private static Logger log = Logger.getLogger(ResourceTree.class);
    private static String XML_ROOT = "Global";
    private List<CheckableResource> leafResources = new ArrayList();
    private HashMap<String, TreeNode> leafNodeMap = new HashMap();
    private JTree tree;
    private Set<String> selectedLeafNodePaths = new LinkedHashSet();

    private static enum TreeExpansionFlag {
        EXPAND_ALL,
        EXPAND_ROOT_ONLY,
        EXPAND_SELECTED_ONLY,
    }

    private ResourceTree() {
    }

    static class CancelableOptionPane extends JOptionPane {
        private boolean canceled = false;

        CancelableOptionPane(Object o, int i, int i1) {
            super(o, i, i1);
        }

        public boolean isCanceled() {
            return canceled;
        }

        public void setCanceled(boolean canceled) {
            this.canceled = canceled;
        }
    }


    /**
     * Shows a tree of selectable resources.
     *
     * @param document    The document that represents an XML resource tree.
     * @param dialogTitle
     * @return the resources selected by user.
     */
    public static LinkedHashSet<ResourceLocator> showResourceTreeDialog(Component parent, Document document, String dialogTitle) {

        JDialog dialog = null;

        final LinkedHashSet<ResourceLocator> locators = new LinkedHashSet();
        try {

            final ResourceTree resourceTree = new ResourceTree();
            final JTree dialogTree = resourceTree.createTreeFromDOM(document);

            int optionType = JOptionPane.OK_CANCEL_OPTION;
            int messageType = JOptionPane.PLAIN_MESSAGE;
            final CancelableOptionPane optionPane =
                    new CancelableOptionPane(new JScrollPane(dialogTree), messageType, optionType);

            optionPane.setPreferredSize(new Dimension(650, 500));
            optionPane.setOpaque(true);
            optionPane.setBackground(Color.WHITE);
            optionPane.addPropertyChangeListener(JOptionPane.VALUE_PROPERTY,
                    new PropertyChangeListener() {

                        public void propertyChange(PropertyChangeEvent e) {

                            Object value = e.getNewValue();
                            if (value instanceof Integer) {
                                int option = (Integer) value;
                                if (option == JOptionPane.CANCEL_OPTION) {
                                    optionPane.setCanceled(true);

                                } else {

                                    LinkedHashSet<ResourceLocator> selectedLocators =
                                            resourceTree.getSelectedResourceLocators();
                                    for (ResourceLocator locator : selectedLocators) {
                                        locators.add(locator);
                                    }
                                }
                            }
                        }
                    });

            dialog = optionPane.createDialog(parent, (dialogTitle == null ? "Resource Tree" : dialogTitle));
            dialog.setBackground(Color.WHITE);
            dialog.getContentPane().setBackground(Color.WHITE);

            Component[] children = optionPane.getComponents();
            if (children != null) {
                for (Component child : children) {
                    child.setBackground(Color.WHITE);
                }
            }

            dialog.setResizable(true);
            dialog.pack();
            dialog.setLocationRelativeTo(parent);
            dialog.setVisible(true);

            return optionPane.isCanceled() ? null : locators;
        } catch (Exception e) {
            log.error(FAILED_TO_CREATE_RESOURCE_TREE_DIALOG, e);
            return null;
        }
    }

    private void initTree(DefaultMutableTreeNode rootNode) {
        tree = new JTree(rootNode);
        tree.setExpandsSelectedPaths(true);
        tree.setCellRenderer(new NodeRenderer());
        tree.setCellEditor(new ResourceEditor(tree));
        tree.setEditable(true);
    }

    private JTree createTreeFromDOM(Document document) {

        Element rootElement =
                (Element) document.getElementsByTagName(XML_ROOT).item(0);

        if (rootElement == null) {
            return new JTree(new DefaultMutableTreeNode(""));
        }

        String nodeName = rootElement.getNodeName();
        if (!nodeName.equalsIgnoreCase(XML_ROOT)) {
            throw new RuntimeException(rootElement +
                    " is not the root of the xml document!");
        }

        String rootLabel = getAttribute(rootElement, "name");
        DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode(rootLabel);

        initTree(rootNode);

        Set<ResourceLocator> loadedResources = IGV.hasInstance() ?
                IGV.getInstance().getDataResourceLocators() : Collections.<ResourceLocator>emptySet();
        loadedResources.addAll(AttributeManager.getInstance().getLoadedResources());

        // Build and attach descendants of the root node to the tree
        buildLocatorTree(rootNode, rootElement, loadedResources, this);

        // Force proper checks on startup
        recheckTree();

        expandTree(TreeExpansionFlag.EXPAND_SELECTED_ONLY, true);
        tree.updateUI(); // Tree state not always correct without this call
        return tree;
    }

    /**
     * Build a tree of all resources, placed under {@code treeNode}, starting
     * from {@code xmlNode}.
     *
     * @param treeNode
     * @param xmlNode
     * @param alreadyLoaded Resources which have already been loaded. These will be checked
     * @param resourceTree  ResourceTree instance. Can be null if not intending UI interaction
     */
    public static void buildLocatorTree(DefaultMutableTreeNode treeNode, Element xmlNode,
                                        Set<ResourceLocator> alreadyLoaded, ResourceTree resourceTree) {

        String name = getAttribute(xmlNode, NAME.getText());

        ResourceLocator locator = new ResourceLocator(
                getAttribute(xmlNode, PATH.getText())
        );
        locator.setName(name);

        String infoLink = getAttribute(xmlNode, HYPERLINK.getText());
        if (infoLink == null) {
            infoLink = getAttribute(xmlNode, INFOLINK.getText());
        }
        locator.setTrackInforURL(infoLink);

        if (xmlNode.getTagName().equalsIgnoreCase("Resource")) {

            String resourceType = getAttribute(xmlNode, RESOURCE_TYPE.getText());
            locator.setType(resourceType);


            String sampleId = getAttribute(xmlNode, SAMPLE_ID.getText());
            if (sampleId == null) {
                // legacy option
                sampleId = getAttribute(xmlNode, ID.getText());
            }
            locator.setSampleId(sampleId);
            locator.setFeatureInfoURL(getAttribute(xmlNode, URL.getText()));
            locator.setDescription(getAttribute(xmlNode, DESCRIPTION.getText()));
            locator.setTrackLine(getAttribute(xmlNode, TRACK_LINE.getText()));
            locator.setName(name);
            // Special element for alignment tracks
            locator.setCoverage(getAttribute(xmlNode, COVERAGE.getText()));

            String colorString = getAttribute(xmlNode, COLOR.getText());
            if (colorString != null) {
                try {
                    Color c = ColorUtilities.stringToColor(colorString);
                    locator.setColor(c);
                } catch (Exception e) {
                    log.error("Error setting color: ", e);
                }
            }
        }

        NodeList nodeList = xmlNode.getChildNodes();
        Node xmlChildNode;

        // If we have children treat it as a category not a leaf
        for (int i = 0; i < nodeList.getLength(); i++) {

            xmlChildNode = nodeList.item(i);
            String nodeName = xmlChildNode.getNodeName();
            if (nodeName.equalsIgnoreCase("#text")) {
                continue;
            }

            // Need to check class of child node, its not necessarily an
            // element (could be a comment for example).
            if (xmlChildNode instanceof Element) {
                String categoryLabel = getAttribute((Element) xmlChildNode, NAME.getText());
                DefaultMutableTreeNode treeChildNode = new DefaultMutableTreeNode(categoryLabel);
                treeNode.add(treeChildNode);
                buildLocatorTree(treeChildNode, (Element) xmlChildNode, alreadyLoaded, resourceTree);
            }
        }

        CheckableResource resource = new CheckableResource(name, false, locator);
        treeNode.setUserObject(resource);

        if (resourceTree != null) {
            resource.setEnabled(resourceTree.tree.isEnabled());

            // If it's a leaf set the checkbox to represent the resource
            if (treeNode.isLeaf()) {

                resourceTree.expandPath(new TreePath(treeNode.getPath()));
                treeNode.setAllowsChildren(false);

                // If data already loaded disable the check box
                if (alreadyLoaded.contains(locator)) {
                    resource.setEnabled(false);
                    resource.setSelected(true);
                    resourceTree.checkParentNode(treeNode, true);
                }
                resourceTree.leafResources.add(resource);
            } else {

                treeNode.setAllowsChildren(true);

                boolean hasSelectedChildren = resourceTree.hasSelectedChildren(treeNode);
                resource.setSelected(hasSelectedChildren);

                if (true || hasSelectedChildren) {
                    ResourceEditor.checkOrUncheckParentNodesRecursively(treeNode, true);
                }
            }

            // Store the paths to all the leaf nodes for easy access
            if (treeNode.isLeaf()) {
                resourceTree.leafNodeMap.put(resourceTree.getPath(treeNode), treeNode);
            }
        }
    }

    public TreeNode checkParentNode(TreeNode childNode, boolean isSelected) {

        TreeNode parentNode = childNode.getParent();
        Object parentsUserObject =
                ((DefaultMutableTreeNode) parentNode).getUserObject();

        if (parentsUserObject instanceof CheckableResource) {
            CheckableResource parentResource =
                    ((CheckableResource) parentsUserObject);

            if (parentResource.isEnabled()) {
                parentResource.setSelected(isSelected);
            }
        }

        return parentNode;
    }

    public List<CheckableResource> getLeafResources() {
        return leafResources;
    }

    public LinkedHashSet<ResourceLocator> getSelectedResourceLocators() {

        LinkedHashSet<ResourceLocator> resourceLocators = new LinkedHashSet();

        for (CheckableResource resource : leafResources) {

            if (resource.isSelected()) {
                resourceLocators.add(resource.getResourceLocator());
            }
        }
        return resourceLocators;
    }

    private static String getAttribute(Element element, String key) {

        String value = element.getAttribute(key);
        if (value != null) {
            if (value.trim().equals("")) {
                value = null;
            }
        }
        return value;
    }


    private LinkedHashSet<DefaultMutableTreeNode> findSelectedLeafNodes(
            TreeModel model, TreeNode parentNode, LinkedHashSet<DefaultMutableTreeNode> allCheckedLeafNodes) {

        if (allCheckedLeafNodes == null) {
            allCheckedLeafNodes = new LinkedHashSet();
        }

        int count = model.getChildCount(parentNode);
        for (int i = 0; i < count; i++) {
            TreeNode childNode = (TreeNode) model.getChild(parentNode, i);
            if (!childNode.isLeaf()) {
                findSelectedLeafNodes(model, childNode, allCheckedLeafNodes);
            } else {
                DefaultMutableTreeNode node = (DefaultMutableTreeNode) childNode;
                CheckableResource resource = (CheckableResource) node.getUserObject();

                if (resource.isSelected()) {
                    allCheckedLeafNodes.add(node);
                }
            }
        }

        return allCheckedLeafNodes;
    }

    /**
     * Node's Renderer
     */
    static class NodeRenderer implements TreeCellRenderer {

        private LinkCheckBox renderer = new LinkCheckBox();
        private Color selectionForeground;
        private Color selectionBackground;
        private Color textForeground;
        private Color textBackground;

        public NodeRenderer() {

            Font fontValue;
            fontValue = UIManager.getFont("Tree.font");
            if (fontValue != null) {
                renderer.setFont(fontValue);
            }
            Boolean booleanValue =
                    (Boolean) UIManager.get("Tree.drawsFocusBorderAroundIcon");
            renderer.setFocusPainted(
                    (booleanValue != null) && (booleanValue.booleanValue()));

            selectionForeground = UIManager.getColor("Tree.selectionForeground");
            selectionBackground = UIManager.getColor("Tree.selectionBackground");
            textForeground = UIManager.getColor("Tree.textForeground");
            textBackground = UIManager.getColor("Tree.textBackground");
            renderer.setSelected(false);
        }

        public Component getTreeCellRendererComponent(JTree tree, Object value,
                                                      boolean isNodeSelected, boolean isNodeExpanded, boolean isLeaf,
                                                      int row, boolean hasFocus) {

            // Convert value into a usable string
            String stringValue = "";
            if (value != null) {
                String toStringValue = value.toString();
                if (toStringValue != null) {
                    stringValue = toStringValue;
                }
            }

            // Initialize checkbox state and selection
            renderer.setSelected(false);
            renderer.setText(stringValue);
            renderer.setEnabled(tree.isEnabled());

            // Tell renderer how to highlight nodes on selection
            if (isNodeSelected) {
                renderer.setForeground(selectionForeground);
                renderer.setBackground(selectionBackground);
            } else {
                renderer.setForeground(textForeground);
                renderer.setBackground(textBackground);
            }

            if (value != null) {
                if (value instanceof DefaultMutableTreeNode) {

                    DefaultMutableTreeNode node =
                            (DefaultMutableTreeNode) value;

                    Object userObject = node.getUserObject();
                    if (userObject instanceof CheckableResource) {

                        CheckableResource resource = (CheckableResource) userObject;
                        renderer.setText(resource.getText());
                        renderer.setSelected(resource.isSelected());
                        renderer.setEnabled(resource.isEnabled());

                        String hyperLink = resource.getResourceLocator().getTrackInfoURL();
                        if (hyperLink == null) {
                            renderer.showHyperLink(false);
                        } else {
                            renderer.setHyperLink(hyperLink);
                            renderer.showHyperLink(true);
                        }
                    }
                }
            }

            return renderer;
        }

        protected LinkCheckBox getRendereringComponent() {
            return renderer;
        }
    }

    /**
     * Node's Resource Editor
     */
    static class ResourceEditor extends AbstractCellEditor
            implements TreeCellEditor {

        NodeRenderer renderer = new NodeRenderer();
        JTree tree;

        public ResourceEditor(JTree tree) {
            this.tree = tree;
        }

        public Object getCellEditorValue() {

            DataResource resource = null;
            TreePath treePath = tree.getEditingPath();
            if (treePath != null) {

                Object node = treePath.getLastPathComponent();

                if ((node != null) && (node instanceof DefaultMutableTreeNode)) {

                    LinkCheckBox checkbox = renderer.getRendereringComponent();

                    DefaultMutableTreeNode treeNode =
                            (DefaultMutableTreeNode) node;

                    Object userObject = treeNode.getUserObject();

                    resource = (CheckableResource) userObject;

                    // Don't change resource if disabled
                    if (!resource.isEnabled()) {
                        return resource;
                    }

                    boolean isChecked = checkbox.isSelected();

                    // Check/Uncheck the selected node. This code ONLY handles
                    // the clicked node. Not it's ancestors or decendants.
                    if (isChecked) {
                        ((CheckableResource) resource).setSelected(true);
                    } else {

                        // See if we are allowed to unchecking this specific 
                        // node - if not, it won't be done. This does not 
                        // prevent it's children from being unchecked.
                        uncheckCurrentNodeIfAllowed((CheckableResource) resource,
                                treeNode);
                    }


                    /*
                    * Now we have to check or uncheck the descendants and
                    * ancestors depending on what we did above.
                    */

                    boolean checkRelatives = isChecked;

                    // If we found a mix of select leave and selected but
                    // but disabled leave we must be trying to toggle off
                    // the children
                    if (hasSelectedAndLockedDescendants(treeNode)) {
                        checkRelatives = false;
                    }
                    // If we found only locked leave we must be trying to toggle 
                    // on the unlocked children
                    else if (hasLockedDescendants(treeNode)) {
                        checkRelatives = true;
                    }
                    // Otherwise, just use the value of the checkbox


                    if (!treeNode.isLeaf()) { //check up and down the tree

                        // If not a leaf check/uncheck children as requested
                        checkOrUncheckChildNodesRecursively(treeNode, checkRelatives);

                        // If not a leaf check/uncheck ancestors
                        checkOrUncheckParentNodesRecursively(treeNode,
                                ((CheckableResource) resource).isSelected());
                    } else { // it must be a leaf - so check up the tree
                        checkOrUncheckParentNodesRecursively(treeNode, checkRelatives);
                    }
                }
                tree.treeDidChange();
            }
            return resource;
        }

        /*
        * Uncheck a node unless rule prevent this behavior.
        */

        private void uncheckCurrentNodeIfAllowed(CheckableResource resource,
                                                 TreeNode treeNode) {

            // If we are unchecking a parent make sure there are
            // no checked children
            if (!hasSelectedChildren(treeNode)) {
                ((CheckableResource) resource).setSelected(false);
            } else {

                // If node has selected children and has disabled descendants we
                // must not unselect
                if (hasLockedDescendants(treeNode)) {
                    ((CheckableResource) resource).setSelected(true);
                } else {
                    // No disabled descendants so we can uncheck at will
                    ((CheckableResource) resource).setSelected(false);
                }
            }
        }

        /**
         * Call to recursively check or uncheck the parent ancestors of the
         * passed node.
         */
        static public void checkOrUncheckParentNodesRecursively(TreeNode node,
                                                                boolean checkParentNode) {

            if (node == null) {
                return;
            }

            TreeNode parentNode = node.getParent();
            if (parentNode == null) {
                return;
            }

            Object parentUserObject =
                    ((DefaultMutableTreeNode) parentNode).getUserObject();

            CheckableResource parentNodeResource = null;
            if (parentUserObject instanceof CheckableResource) {
                parentNodeResource = ((CheckableResource) parentUserObject);
            }

            if (parentNodeResource != null) {

                // If parent's current check state matchs what we want there
                // is nothing to do so just leave
                if (parentNodeResource.isSelected() == checkParentNode) {
                    return;
                } else if (checkParentNode) {
                    parentNodeResource.setSelected(true);
                } else { // Uncheck Only if their are no selected descendants

                    if (!hasSelectedChildren(parentNode)) {
                        parentNodeResource.setSelected(false);
                    }
                }
            }

            checkOrUncheckParentNodesRecursively(parentNode,
                    checkParentNode);
        }

        /**
         * Can only be called from getCellEditorValue() to recursively check
         * or uncheck the children of the passed parent node.
         */
        private void checkOrUncheckChildNodesRecursively(TreeNode currentNode,
                                                         boolean isCheckingNeeded) {

            Object parentUserObject =
                    ((DefaultMutableTreeNode) currentNode).getUserObject();

            CheckableResource currentTreeNodeResource = null;
            if (parentUserObject instanceof CheckableResource) {
                currentTreeNodeResource = ((CheckableResource) parentUserObject);
            }

            if (currentTreeNodeResource != null) {

                // Set all enabled children to the checked state of their parent
                Enumeration children = currentNode.children();
                while (children.hasMoreElements()) {

                    TreeNode childNode = (TreeNode) children.nextElement();

                    Object childsUserObject =
                            ((DefaultMutableTreeNode) childNode).getUserObject();
                    if (childsUserObject instanceof CheckableResource) {

                        CheckableResource childResource =
                                ((CheckableResource) childsUserObject);

                        if (childResource.isEnabled()) {

                            // Child must be checked if it has selected
                            // selected and disabled descendants
                            if (hasLockedDescendants(childNode)) {
                                childResource.setSelected(true);
                            } else { // else check/uncheck  as requested
                                childResource.setSelected(isCheckingNeeded);
                            }
                        }
                    }
                    checkOrUncheckChildNodesRecursively(childNode,
                            isCheckingNeeded);
                }
            }
        }

        public boolean hasLockedDescendants(TreeNode treeNode) {

            Enumeration children = treeNode.children();
            while (children.hasMoreElements()) {

                TreeNode childNode = (TreeNode) children.nextElement();

                Object childsUserObject =
                        ((DefaultMutableTreeNode) childNode).getUserObject();
                if (childsUserObject instanceof CheckableResource) {
                    CheckableResource childResource =
                            ((CheckableResource) childsUserObject);

                    // If disabled say so
                    if (!childResource.isEnabled()) {
                        return true;
                    }
                }

                // If a descendant is disabled say so
                if (hasLockedDescendants(childNode)) {
                    return true;
                }
            }
            return false;
        }

        static public boolean hasSelectedDescendants(TreeNode treeNode) {

            Enumeration children = treeNode.children();
            while (children.hasMoreElements()) {

                TreeNode childNode = (TreeNode) children.nextElement();

                Object childsUserObject =
                        ((DefaultMutableTreeNode) childNode).getUserObject();
                if (childsUserObject instanceof CheckableResource) {
                    CheckableResource childResource =
                            ((CheckableResource) childsUserObject);

                    // If has selected say so
                    if (childResource.isSelected()) {
                        return true;
                    }
                }

                // If has selected descendant say so
                if (hasSelectedDescendants(childNode)) {
                    return true;
                }
            }
            return false;
        }

        static public boolean hasSelectedChildren(TreeNode treeNode) {

            Enumeration children = treeNode.children();
            while (children.hasMoreElements()) {

                TreeNode childNode = (TreeNode) children.nextElement();

                Object childsUserObject =
                        ((DefaultMutableTreeNode) childNode).getUserObject();
                if (childsUserObject instanceof CheckableResource) {
                    CheckableResource childResource =
                            ((CheckableResource) childsUserObject);
                    if (childResource.isSelected()) {
                        return true;
                    }
                }
            }
            return false;
        }

        /**
         * Return true if it find nodes that ar both selected and disabled
         *
         * @param treeNode
         * @return true if we are working with preselected nodes
         */
        public boolean hasLockedChildren(TreeNode treeNode) {

            boolean hasSelectedAndDisabled = false;
            Enumeration children = treeNode.children();
            while (children.hasMoreElements()) {

                TreeNode childNode = (TreeNode) children.nextElement();

                Object childsUserObject =
                        ((DefaultMutableTreeNode) childNode).getUserObject();
                if (childsUserObject instanceof CheckableResource) {

                    CheckableResource childResource =
                            ((CheckableResource) childsUserObject);

                    if (!childResource.isEnabled() && childResource.isSelected()) {
                        hasSelectedAndDisabled = true;
                    }

                    if (hasSelectedAndDisabled) {
                        break;
                    }
                }
            }
            return (hasSelectedAndDisabled);
        }

        /**
         * @param treeNode
         * @return true if we are working with preselected nodes
         */
        public boolean hasSelectedAndLockedChildren(TreeNode treeNode) {

            boolean hasSelected = false;
            boolean hasSelectedAndDisabled = false;
            Enumeration children = treeNode.children();
            while (children.hasMoreElements()) {

                TreeNode childNode = (TreeNode) children.nextElement();

                Object childsUserObject =
                        ((DefaultMutableTreeNode) childNode).getUserObject();
                if (childsUserObject instanceof CheckableResource) {

                    CheckableResource childResource =
                            ((CheckableResource) childsUserObject);

                    if (childResource.isSelected() && childResource.isEnabled()) {
                        hasSelected = true;
                    }
                    if (!childResource.isEnabled() && childResource.isSelected()) {
                        hasSelectedAndDisabled = true;
                    }

                    if (hasSelected & hasSelectedAndDisabled) {
                        break;
                    }
                }
            }

            // If we have both we can return true
            return (hasSelected & hasSelectedAndDisabled);
        }

        /**
         * @param treeNode
         * @return true if we are working with preselected nodes
         */
        public boolean hasSelectedAndLockedDescendants(TreeNode treeNode) {

            boolean hasSelected = false;
            boolean hasSelectedAndDisabled = false;
            Enumeration children = treeNode.children();
            while (children.hasMoreElements()) {

                TreeNode childNode = (TreeNode) children.nextElement();

                Object childsUserObject =
                        ((DefaultMutableTreeNode) childNode).getUserObject();
                if (childsUserObject instanceof CheckableResource) {

                    CheckableResource childResource =
                            ((CheckableResource) childsUserObject);

                    if (childResource.isSelected() && childResource.isEnabled()) {
                        hasSelected = true;
                    }
                    if (!childResource.isEnabled() && childResource.isSelected()) {
                        hasSelectedAndDisabled = true;
                    }

                    if (hasSelected & hasSelectedAndDisabled) {
                        break;
                    }
                }

                // If has a mix of selected and checked but disableddescendant
                if (hasSelectedAndLockedDescendants(childNode)) {
                    return true;
                }

            }

            // If we have both we can return true
            return (hasSelected & hasSelectedAndDisabled);
        }

        @Override
        public boolean isCellEditable(EventObject event) {

            boolean returnValue = false;

            if (event instanceof MouseEvent) {

                MouseEvent mouseEvent = (MouseEvent) event;
                TreePath treePath = tree.getPathForLocation(
                        mouseEvent.getX(), mouseEvent.getY());

                if (treePath != null) {

                    Object node = treePath.getLastPathComponent();
                    if ((node != null) &&
                            (node instanceof DefaultMutableTreeNode)) {
                        DefaultMutableTreeNode treeNode =
                                (DefaultMutableTreeNode) node;
                        Object userObject = treeNode.getUserObject();

                        if (userObject instanceof CheckableResource) {
                            returnValue = true;
                        } else if (userObject instanceof CheckableResource) {
                            returnValue = true;
                        }
                    }
                }
            }
            return returnValue;
        }

        public Component getTreeCellEditorComponent(JTree tree, Object value,
                                                    boolean selected, boolean expanded, boolean leaf, int row) {

            Component rendererComponent = renderer.getTreeCellRendererComponent(
                    tree, value, true, expanded, leaf, row, true);


            ItemListener itemListener = new ItemListener() {

                public void itemStateChanged(ItemEvent itemEvent) {
                    if (stopCellEditing()) {
                        fireEditingStopped();
                    }
                }
            };
            if (rendererComponent instanceof LinkCheckBox) {
                ((LinkCheckBox) rendererComponent).addItemListener(itemListener);
            }

            return rendererComponent;
        }
    }

    public static class CheckableResource implements SelectableResource {

        final static protected Color partialSelectionColor =
                new Color(255, 128, 128);
        protected boolean isParentOfPartiallySelectedChildren = false;
        protected String text;
        protected boolean selected;
        protected ResourceLocator dataResourceLocator;
        protected boolean isEnabled = true;

        public CheckableResource() {
        }

        public CheckableResource(String text, boolean selected,
                                 ResourceLocator dataResourceLocator) {

            this.text = text;
            this.selected = selected;
            this.dataResourceLocator = dataResourceLocator;
        }

        public boolean isSelected() {
            return selected;
        }

        public void setSelected(boolean newValue) {
            selected = newValue;
        }

        public boolean isEnabled() {
            return isEnabled;
        }

        public void setEnabled(boolean value) {
            isEnabled = value;
        }

        public String getText() {
            return text;
        }

        public void setText(String newValue) {
            text = newValue;
        }

        public ResourceLocator getResourceLocator() {
            return dataResourceLocator;
        }

        public void setResourceLocator(ResourceLocator dataResourceLocator) {
            this.dataResourceLocator = dataResourceLocator;
        }

        public boolean isParentOfPartiallySelectedChildren() {
            return isParentOfPartiallySelectedChildren;
        }

        public void setIsParentOfPartiallySelectedChildren(boolean value) {
            this.isParentOfPartiallySelectedChildren = value;
        }

        public Color getBackground() {

            if (isParentOfPartiallySelectedChildren()) {
                return partialSelectionColor;
            } else {
                return Color.WHITE;
            }
        }

        @Override
        public String toString() {
            return text + ":" + selected;
        }
    }


    /**
     *
     */
    static interface SelectableResource extends DataResource {

        public boolean isSelected();

        public void setSelected(boolean newValue);
    }

    /**
     *
     */
    static interface DataResource {

        public ResourceLocator getResourceLocator();

        public void setText(String newValue);

        public String getText();

        public void setEnabled(boolean value);

        public boolean isEnabled();
    }

    /**
     * Expands tree.
     */
    private void expandTree(TreeExpansionFlag expansionFlag,
                            boolean skipSelectingEnabledNodes) {

        TreeNode root = (TreeNode) tree.getModel().getRoot();
        if (expansionFlag == TreeExpansionFlag.EXPAND_SELECTED_ONLY) {

            if (selectedLeafNodePaths.isEmpty()) {
                TreePath rootPath = new TreePath(root);
                expandPath(rootPath);
            } else {
                boolean[] expansionVetoed = {true}; // Default to not expanded
                expandSelectedDescendants(new TreePath(root), true,
                        expansionVetoed, skipSelectingEnabledNodes);
            }
        } else if (expansionFlag == TreeExpansionFlag.EXPAND_ALL) {
            expandAllDescendants(new TreePath(root), true,
                    skipSelectingEnabledNodes);
        } else if (expansionFlag == TreeExpansionFlag.EXPAND_ROOT_ONLY) {
            TreePath rootPath = new TreePath(root);
            expandPath(rootPath);
            checkAllSelectedLeafNodes(rootPath, skipSelectingEnabledNodes);
        }
    }

    /**
     * Expands all children of a tree node.
     */
    private void expandAllDescendants(TreePath parentPath, boolean isExpanding,
                                      boolean skipSelectingEnabledNodes) {

        // Traverse children
        TreeNode node = (TreeNode) parentPath.getLastPathComponent();
        for (Enumeration e = node.children(); e.hasMoreElements(); ) {
            TreeNode n = (TreeNode) e.nextElement();
            TreePath childPath = parentPath.pathByAddingChild(n);
            expandAllDescendants(childPath, isExpanding,
                    skipSelectingEnabledNodes);
        }

        // Expansion or collapse must be done bottom-up
        if (isExpanding) {
            expandPath(parentPath);
        } else {
            collapsePath(parentPath);
        }

        // Leaf nodes processing 
        if (node.isLeaf()) {
            String path = getPath((DefaultMutableTreeNode) node);
            if (selectedLeafNodePaths.contains(path)) {
                manuallySelectNode((DefaultMutableTreeNode) node,
                        skipSelectingEnabledNodes);
            }
        }
    }

    /**
     * Expands all children of a tree node.
     */
    private void expandSelectedDescendants(TreePath parentPath,
                                           boolean isExpanding, boolean[] expansionVetoed,
                                           boolean skipSelectingEnabledNodes) {

        boolean[] expansionIsVetoed = {true}; // true: Defaults to not expanded

        // Traverse children
        TreeNode node = (TreeNode) parentPath.getLastPathComponent();
        for (Enumeration e = node.children(); e.hasMoreElements(); ) {

            DefaultMutableTreeNode childNode =
                    (DefaultMutableTreeNode) e.nextElement();

            // Leaf nodes processing 
            if (childNode.isLeaf()) {

                String path = getPath(childNode);
                if (selectedLeafNodePaths.contains(path)) {
                    manuallySelectNode(childNode, skipSelectingEnabledNodes);
                    expansionIsVetoed[0] = false;
                }
            }
            TreePath childPath = parentPath.pathByAddingChild(childNode);
            expandSelectedDescendants(childPath, isExpanding, expansionVetoed,
                    skipSelectingEnabledNodes);
        }

        if (expansionIsVetoed[0])
            return;

        // Expansion or collapse must be done bottom-up
        if (isExpanding) {
            expandPath(parentPath);
        } else {
            collapsePath(parentPath);
        }
    }

    private String getPath(DefaultMutableTreeNode treeNode) {

        buffer.delete(0, buffer.length());

        TreeNode[] nodesInPath = treeNode.getPath();
        for (Object element : nodesInPath) {
            DefaultMutableTreeNode node = (DefaultMutableTreeNode) element;
            Object userObject = node.getUserObject();
            if (userObject instanceof CheckableResource) {
                CheckableResource resource = (CheckableResource) userObject;
                buffer.append(resource.getText());

                if (!node.isLeaf()) {
                    buffer.append("/");
                }
            } else if (userObject instanceof CheckableResource) {
                CheckableResource resource = (CheckableResource) userObject;
                buffer.append(resource.getText());

                if (!node.isLeaf()) {
                    buffer.append("/");
                }
            } else {
                buffer.append(userObject.toString());

                if (!node.isLeaf()) {
                    buffer.append("/");
                }
            }
        }

        return buffer.toString();
    }


    private void manuallySelectNode(DefaultMutableTreeNode childNode,
                                    boolean skipSelectingEnabledNodes) {


        Object userObject = childNode.getUserObject();
        if (userObject instanceof CheckableResource) {
            CheckableResource resource = (CheckableResource) userObject;

            // Leave if node is disable
            if (skipSelectingEnabledNodes && resource.isEnabled()) {
                return;
            }

            // Don't select the node if already selected
            if (resource.isSelected())
                return;

            resource.setSelected(true);

            if (childNode.isLeaf()) {
                TreeNode parentNode = childNode.getParent();
                Object parentsUserObject =
                        ((DefaultMutableTreeNode) parentNode).getUserObject();
                if (parentsUserObject instanceof CheckableResource) {
                    CheckableResource parentResource =
                            ((CheckableResource) parentsUserObject);
                    if (!parentResource.isSelected()) {
                        parentResource.setSelected(true);
                    }
                    tree.treeDidChange();
                }
            }
        }
    }

    private void checkAllSelectedLeafNodes(TreePath parentPath,
                                           boolean skipSelectingEnabledNodes) {

        // Traverse children
        TreeNode node = (TreeNode) parentPath.getLastPathComponent();
        for (Enumeration e = node.children(); e.hasMoreElements(); ) {
            TreeNode n = (TreeNode) e.nextElement();
            TreePath childPath = parentPath.pathByAddingChild(n);
            checkAllSelectedLeafNodes(childPath,
                    skipSelectingEnabledNodes);
        }

        // Leaf nodes processing 
        if (node.isLeaf()) {
            String path = getPath((DefaultMutableTreeNode) node);
            if (selectedLeafNodePaths.contains(path)) {
                manuallySelectNode((DefaultMutableTreeNode) node,
                        skipSelectingEnabledNodes);
            }
        }
    }

    private void collapsePath(TreePath treePath) {
        if (!tree.isCollapsed(treePath)) {
            tree.collapsePath(treePath);
        }
    }

    private void expandPath(TreePath treePath) {
        if (!tree.isExpanded(treePath)) {
            tree.expandPath(treePath);
        }
    }

    public boolean hasSelectedChildren(TreeNode treeNode) {

        Enumeration children = treeNode.children();
        while (children.hasMoreElements()) {

            TreeNode childNode = (TreeNode) children.nextElement();

            Object childsUserObject =
                    ((DefaultMutableTreeNode) childNode).getUserObject();
            if (childsUserObject instanceof CheckableResource) {
                CheckableResource childResource =
                        ((CheckableResource) childsUserObject);
                if (childResource.isSelected()) {
                    return true;
                }
                tree.treeDidChange();
            }
        }
        return false;
    }

    static private Set<ResourceLocator> getLoadedResources() {
        Set<ResourceLocator> loadedResources = IGV.getInstance().getDataResourceLocators();
        loadedResources.addAll(AttributeManager.getInstance().getLoadedResources());
        return loadedResources;
    }

    /**
     * This method will only check tree items.
     */
    private void recheckTree() {

        if (leafNodeMap != null && !leafNodeMap.isEmpty()) {
            Collection<TreeNode> leaves = leafNodeMap.values();
            for (TreeNode leafNode : leaves) {
                TreeNode parent = leafNode.getParent();
                if (parent != null) {

                    Object userObject =
                            ((DefaultMutableTreeNode) leafNode).getUserObject();
                    CheckableResource checkableLeafResource =
                            ((CheckableResource) userObject);
                    if (userObject != null) {

                        if (checkableLeafResource.isSelected()) {
                            checkNode(parent, true);
                            while (true) {
                                parent = parent.getParent();
                                if (parent == null) {
                                    break;
                                }
                                checkNode(parent, true);
                            }
                        }
                    }
                }
            }
        }
    }

    protected void checkNode(TreeNode node, boolean checked) {

        // Only allowed to check
        if (!checked)
            return;

        Object userObject = ((DefaultMutableTreeNode) node).getUserObject();
        CheckableResource nodeResource = null;
        if (userObject instanceof CheckableResource) {
            nodeResource = ((CheckableResource) userObject);
            nodeResource.setSelected(checked);
        }
    }

}
