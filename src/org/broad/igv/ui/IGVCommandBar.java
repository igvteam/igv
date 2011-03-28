/*
 * Copyright (c) 2007-2011 by The Broad Institute, Inc. and the Massachusetts Institute of
 * Technology.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */
/*
 * IGVCommandBar.java
 *
 * Created on April 5, 2008, 10:02 AM
 */
package org.broad.igv.ui;


import com.jidesoft.swing.JideBoxLayout;
import com.jidesoft.swing.JideButton;
import com.jidesoft.swing.JideToggleButton;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.PreferenceManager;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.genome.GenomeManager.GenomeListItem;
import org.broad.igv.feature.genome.GenomeServerException;
import org.broad.igv.session.History;
import org.broad.igv.ui.action.FitDataToWindowMenuAction;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.action.SearchCommand;
import org.broad.igv.ui.panel.ZoomSliderPanel;
import org.broad.igv.ui.util.*;
import org.broad.igv.ui.util.ProgressMonitor;
import org.broad.igv.util.NamedRunnable;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.LineBorder;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.io.FileNotFoundException;
import java.net.NoRouteToHostException;
import java.util.*;
import java.util.List;

/**
 * @author jrobinso
 */
public class IGVCommandBar extends javax.swing.JPanel {

    private static Logger log = Logger.getLogger(IGVCommandBar.class);
    private IGVMainFrame owner;
    private LinkedHashSet<GenomeListItem> userDefinedGenomeItemList;
    private LinkedHashSet<GenomeListItem> serverGenomeItemList;
    private LinkedHashSet<GenomeListItem> cachedGenomeItemList;
    private JComboBox chromosomeComboBox;
    private JComboBox genomeComboBox;
    //private JPanel geneListPanel;
    // private JideButton geneListLabel;
    private JideButton goButton;
    private JideButton homeButton;
    private JPanel locationPanel;
    private JideButton refreshButton;
    private JideToggleButton roiToggleButton;
    private JTextField searchTextField;
    private JPanel toolPanel;
    private JPanel zoomControl;
    private ItemListener genomeComboBoxItemListener = null;
    final private int DEFAULT_CHROMOSOME_DROPDOWN_WIDTH = 120;
    private JideButton backButton;
    private JideButton forwardButton;
    private JideButton fitToWindowButton;
    private static final Font GENE_LIST_FONT = FontManager.getScalableFont(Font.BOLD, 14);

    /**
     * Creates new form IGVCommandBar
     *
     * @param owner
     */
    public IGVCommandBar(IGVMainFrame owner) {
        this.owner = owner;
        initComponents();

        // Initialize controls

        String currentChr = getDefaultReferenceFrame().getChrName();
        boolean isWholeGenome = currentChr.equals(Globals.CHR_ALL);

        chromosomeComboBox.setSelectedItem(currentChr);
        roiToggleButton.setEnabled(!isWholeGenome);
        zoomControl.setEnabled(!isWholeGenome);
    }

    /**
     * Method description
     *
     * @param genome
     * @return
     */
    public boolean isGenomeCached(String genome) {
        boolean isCached = false;

        if ((cachedGenomeItemList != null) && !cachedGenomeItemList.isEmpty()) {
            for (GenomeListItem item : cachedGenomeItemList) {
                if (item.getId().equalsIgnoreCase(genome)) {
                    isCached = true;
                }
            }
        }
        return isCached;
    }

    /**
     * This method is called once on startup
     *
     * @param monitor
     * @throws FileNotFoundException
     * @throws NoRouteToHostException
     */
    public void initializeGenomeList(final ProgressMonitor monitor)
            throws FileNotFoundException, NoRouteToHostException {

        if (log.isDebugEnabled()) {
            log.debug("Enter initializeGenomeList");
        }

        if (monitor != null) {
            monitor.fireProgressChange(1);
        }

        genomeComboBox.removeAllItems();
        genomeComboBox.setRenderer(new ComboBoxRenderer());
        genomeComboBox.setToolTipText(UIConstants.CHANGE_GENOME_TOOLTIP);
        rebuildGenomeItemList(null);

        if (monitor != null) {
            monitor.fireProgressChange(50);
        }

        genomeComboBox.addItemListener(getGenomeComboBoxItemListener());

        // Initialize the genome list
        LinkedHashSet<GenomeListItem> mergedServerAndCacheItemList = new LinkedHashSet();
        if (this.userDefinedGenomeItemList != null) {
            mergedServerAndCacheItemList.addAll(this.userDefinedGenomeItemList);
        }
        if (this.cachedGenomeItemList != null) {
            mergedServerAndCacheItemList.addAll(this.cachedGenomeItemList);
        }
        if (this.serverGenomeItemList != null) {
            mergedServerAndCacheItemList.addAll(this.serverGenomeItemList);
        }

        // Post creation widget setup.
        searchTextField.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                goButtonActionPerformed(actionevent);
            }
        });

        if (log.isDebugEnabled()) {
            log.debug("Exit initializeGenomeList");
        }

    }

    /**
     * Method description
     *
     * @return
     */
    public ItemListener getGenomeComboBoxItemListener() {

        if (genomeComboBoxItemListener == null) {

            genomeComboBoxItemListener = new ItemListener() {

                int previousSelectedItemId = 0;

                private void reselectPreviousGenome() {
                    if (previousSelectedItemId != -1) {

                        final int index = previousSelectedItemId;
                        UIUtilities.invokeOnEventThread(new Runnable() {

                            public void run() {
                                genomeComboBox.setSelectedIndex(index);
                            }
                        });
                    }
                    previousSelectedItemId = -1;
                }

                public void itemStateChanged(final ItemEvent itemEvent) {

                    // Only interested in selections
                    if (ItemEvent.SELECTED != itemEvent.getStateChange()) {
                        return;
                    }

                    SwingWorker worker = new SwingWorker() {

                        @Override
                        protected Object doInBackground() throws Exception {
                            Object item = itemEvent.getItem();

                            if (item instanceof GenomeListItem) {
                                GenomeListItem genomeListItem = (GenomeListItem) item;
                                String requestedGenomeId = genomeListItem.getId();
                                String currentGenomeId = GenomeManager.getInstance().getGenomeId();
                                if ((currentGenomeId != null) && requestedGenomeId.equalsIgnoreCase(currentGenomeId)) {

                                    // Nothing to do if genome already loaded
                                    return null;
                                }

                                final ProgressMonitor monitor = new ProgressMonitor();
                                final ProgressBar bar =
                                        ProgressBar.showProgressDialog(IGVMainFrame.getInstance(),
                                                "Loading Genome...", monitor,
                                                false);


                                if (requestedGenomeId != null) {
                                    try {
                                        monitor.fireProgressChange(50);

                                        GenomeManager.getInstance().loadGenome(genomeListItem, genomeListItem.isUserDefined());

                                        monitor.fireProgressChange(25);

                                        if (!isGenomeCached(genomeListItem.getId())) {
                                            cachedGenomeItemList.add(genomeListItem);
                                        }

                                        IGVMainFrame.getInstance().doChooseGenome(
                                                GenomeManager.getInstance().getGenomeDescriptor(
                                                        requestedGenomeId));

                                        previousSelectedItemId = genomeComboBox.getSelectedIndex();

                                        monitor.fireProgressChange(25);

                                        UIUtilities.invokeOnEventThread(new Runnable() {

                                            public void run() {
                                                if (bar != null) {
                                                    bar.close();
                                                }
                                            }
                                        });
                                    } catch (GenomeServerException e) {

                                        UIUtilities.invokeOnEventThread(new Runnable() {

                                            public void run() {
                                                JOptionPane.showMessageDialog(
                                                        IGVMainFrame.getInstance(),
                                                        UIConstants.CANNOT_ACCESS_SERVER_GENOME_LIST);
                                            }
                                        });
                                        reselectPreviousGenome();
                                        return null;
                                    } catch (FileNotFoundException e) {
                                        UIUtilities.invokeOnEventThread(new Runnable() {

                                            public void run() {
                                                if (bar != null) {
                                                    bar.close();
                                                }
                                            }
                                        });

                                        int choice =
                                                JOptionPane.showConfirmDialog(
                                                        IGVMainFrame.getInstance(), "The genome file [" + e.getMessage() + "] could not be located. " + "Would you like to remove the selected entry?", "", JOptionPane.OK_CANCEL_OPTION);

                                        if (choice == JOptionPane.OK_OPTION) {
                                            Set<String> excludedArchivesUrls = new HashSet();
                                            excludedArchivesUrls.add(genomeListItem.getLocation());
                                            IGVCommandBar.this.rebuildGenomeItemList(excludedArchivesUrls);
                                        }
                                        reselectPreviousGenome();
                                        return null;
                                    } catch (Exception e) {
                                        reselectPreviousGenome();
                                        return null;
                                    } finally {
                                        IGVMainFrame.getInstance().repaint();
                                    }
                                }
                            }

                            return null;
                        }
                    };
                    worker.execute();
                }
            };
        }
        return genomeComboBoxItemListener;
    }

    /**
     * Adds the new user-defined genome to the drop down list.
     *
     * @param newItem
     */
    public void addToUserDefinedGenomeItemList(GenomeListItem newItem) {


        if (userDefinedGenomeItemList == null) {
            userDefinedGenomeItemList = new LinkedHashSet<GenomeListItem>();
            userDefinedGenomeItemList.add(newItem);
        } else {

            LinkedHashSet tempItemList = new LinkedHashSet<GenomeListItem>();
            tempItemList.add(newItem);
            for (GenomeListItem item : userDefinedGenomeItemList) {
                tempItemList.add(item);
            }
            userDefinedGenomeItemList = tempItemList;
        }
        setGenomeItemListModel();
    }

    /**
     * Completely rebuild the genome drop down info from scratch.
     *
     * @param excludedArchivesUrls
     */
    public void rebuildGenomeItemList(Set excludedArchivesUrls) {

        try {

            // Build a single available genome list from both client, server
            // and cached information. This allows us to process
            // everything the same way.
            LinkedHashSet<GenomeListItem> serverSideItemList = null;

            try {

                serverSideItemList =
                        GenomeManager.getInstance().getServerGenomeArchiveList(excludedArchivesUrls);
            } catch (Exception e) {

                UIUtilities.invokeOnEventThread(new Runnable() {

                    public void run() {
                        JOptionPane.showMessageDialog(
                                IGVMainFrame.getInstance(),
                                UIConstants.CANNOT_ACCESS_SERVER_GENOME_LIST);
                    }
                });
            }

            LinkedHashSet<GenomeListItem> cacheGenomeItemList =
                    GenomeManager.getInstance().getCachedGenomeArchiveList();

            LinkedHashSet<GenomeListItem> clientSideItemList =
                    GenomeManager.getInstance().getUserDefinedGenomeArchiveList();

            setGenomeItemList(clientSideItemList, serverSideItemList, cacheGenomeItemList);
            setGenomeItemListModel();

        } catch (Exception e) {
            log.error("Failed to get genome archive list " + "information from the server!", e);
        }
    }


    void updateChromosomeDropdown() {

        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (genome == null) return;

        List<String> tmp = new LinkedList(genome.getChromosomeNames());
        if (tmp.size() > 1) {
            String homeChr = genome.getHomeChromosome();
            if (homeChr.equals(Globals.CHR_ALL)) {
                tmp.add(Globals.CHR_ALL);
            }
        }

        Graphics2D graphics2D = (Graphics2D) chromosomeComboBox.getGraphics();
        Font font = chromosomeComboBox.getFont();
        FontMetrics fontMetrics = chromosomeComboBox.getFontMetrics(font);

        int w = DEFAULT_CHROMOSOME_DROPDOWN_WIDTH;
        for (String chromosomeName : tmp) {
            Rectangle2D textBounds = fontMetrics.getStringBounds(chromosomeName,
                    graphics2D);
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

        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                adjustChromosomeDropdownWidth(dropdownWidth);
                chromosomeComboBox.setModel(defaultModel);
                chromosomeComboBox.setSelectedItem(getDefaultReferenceFrame().getChrName());
            }
        });

    }

    protected void chromosomeChanged() {
        roiToggleButton.setEnabled(!getDefaultReferenceFrame().getChrName().equals(Globals.CHR_ALL));
        zoomControl.setEnabled(!getDefaultReferenceFrame().getChrName().equals(Globals.CHR_ALL));

        if (chromosomeComboBox.getSelectedItem() != null) {
            if (!chromosomeComboBox.getSelectedItem().equals(getDefaultReferenceFrame().getChrName())) {
                chromosomeComboBox.setSelectedItem(getDefaultReferenceFrame().getChrName());
            }
        }
    }

    /**
     * Method description
     */
    public void updateCurrentCoordinates() {
        searchTextField.setText("");

        final String chr = getDefaultReferenceFrame().getChrName();

        if (!chr.equals(chromosomeComboBox.getSelectedItem())) {
            chromosomeChanged();
            chromosomeComboBox.setSelectedItem(chr);
            owner.chromosomeChangeEvent(false);
        }

        String p = "";

        if (!chr.equals(Globals.CHR_ALL)) {
            p = getDefaultReferenceFrame().getFormattedLocusString();
        }
        final String position = p;
        UIUtilities.invokeOnEventThread(new Runnable() {

            public void run() {
                searchTextField.setText(position);
            }
        });

        final History history = owner.getSession().getHistory();
        forwardButton.setEnabled(history.canGoForward());
        backButton.setEnabled(history.canGoBack());


    }

    private ReferenceFrame getDefaultReferenceFrame() {
        return FrameManager.getDefaultFrame();
    }

    public void setGeneListMode(boolean geneListMode) {

        genomeComboBox.setEnabled(!geneListMode);
        locationPanel.setEnabled(!geneListMode);
        chromosomeComboBox.setEnabled(!geneListMode);
        searchTextField.setEnabled(!geneListMode);
        goButton.setEnabled(!geneListMode);
        zoomControl.setEnabled(!geneListMode);
        homeButton.setEnabled(!geneListMode);
        roiToggleButton.setEnabled(!geneListMode);
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
                renderer = label;
            }

            if (isSelected) {
                renderer.setBackground(list.getSelectionBackground());
                renderer.setForeground(list.getSelectionForeground());
            } else {
                renderer.setBackground(list.getBackground());
                renderer.setForeground(list.getForeground());
            }

            renderer.setFont(list.getFont());

            return renderer;
        }
    }

    /**
     * Selects the first genome from the list,
     */
    public void selectFirstGenomeInList() {

        int itemCount = genomeComboBox.getItemCount();

        for (int i = 0; i < itemCount; i++) {
            Object object = genomeComboBox.getItemAt(i);

            if (object instanceof GenomeListItem) {
                GenomeListItem genomeListItem = (GenomeListItem) object;
                genomeComboBox.setSelectedIndex(i);
                String id = genomeListItem.getId();
                IGVMainFrame.getInstance().setGenomeId(id);
                break;
            }
        }
    }

    /**
     * Gets the collection of genome display names currently in use.
     *
     * @return Set of display names.
     */
    public Collection<String> getGenomeDisplayNames() {

        Set<String> displayNames = new HashSet();
        int itemCount = genomeComboBox.getItemCount();

        for (int i = 0; i < itemCount; i++) {
            Object object = genomeComboBox.getItemAt(i);
            if (object instanceof GenomeListItem) {
                GenomeListItem genomeListItem = (GenomeListItem) object;
                displayNames.add(genomeListItem.getDisplayableName());
            }
        }
        return displayNames;
    }

    /**
     * Gets the collection of genome list items ids currently in use.
     *
     * @return Set of ids.
     */
    public Collection<String> getGenomeIds() {

        Set<String> ids = new HashSet();
        int itemCount = genomeComboBox.getItemCount();

        for (int i = 0; i < itemCount; i++) {
            Object object = genomeComboBox.getItemAt(i);
            if (object instanceof GenomeListItem) {
                GenomeListItem genomeListItem = (GenomeListItem) object;
                ids.add(genomeListItem.getId());
            }
        }
        return ids;
    }

    /**
     * Selects the first genome from the list,
     *
     * @param genomeId
     */
    public void selectGenomeFromListWithNoImport(String genomeId) {

        int itemCount = genomeComboBox.getItemCount();

        for (int i = 0; i < itemCount; i++) {
            Object object = genomeComboBox.getItemAt(i);

            if (object instanceof GenomeListItem) {
                GenomeListItem genomeListItem = (GenomeListItem) object;

                // If the list genome matchs the one we are interested in
                // process it
                String id = genomeListItem.getId();
                if ((id != null) && id.trim().equalsIgnoreCase(genomeId)) {
                    genomeComboBox.setSelectedIndex(i);
                    IGVMainFrame.getInstance().setGenomeId(id);
                    break;
                }
            }
        }
    }


    /**
     * Called from session loading,  command line listener, startup code
     */
    public void selectGenomeFromList(final String genomeId)
            throws FileNotFoundException, NoRouteToHostException {


        // See if this genome is already loaded
        String currentGenomeId = GenomeManager.getInstance().getGenomeId();
        if (currentGenomeId != null && genomeId != null && genomeId.equalsIgnoreCase(currentGenomeId)) {
            return;
        }

        NamedRunnable runnable = new NamedRunnable() {

            public void run() {

                log.debug("Run selectGenomeFromList");

                boolean ignoreCache = false;
                boolean wasFound = false;

                // Now select this item in tne comboBox

                if (genomeId != null) {
                    try {
                        int itemCount = genomeComboBox.getItemCount();
                        for (int i = 0; i < itemCount; i++) {
                            Object object = genomeComboBox.getItemAt(i);

                            if (object instanceof GenomeListItem) {

                                GenomeListItem genomeListItem = (GenomeListItem) object;
                                String id = genomeListItem.getId();

                                // If the list genome matchs the one we are interested in
                                // process it
                                if ((id != null) && id.trim().equalsIgnoreCase(genomeId)) {


                                    if (log.isDebugEnabled()) {
                                        log.debug("Call loadGenome");
                                    }

                                    GenomeManager.getInstance().loadGenome(
                                            genomeListItem.getLocation(),
                                            genomeListItem.isUserDefined(),
                                            null);

                                    genomeComboBox.setSelectedIndex(i);
                                    wasFound = true;

                                    IGVMainFrame.getInstance().setGenomeId(genomeId);

                                    break;
                                }
                            }
                        }
                    } catch (FileNotFoundException ex) {
                        log.error("Genome file not found. ", ex);
                    }
                }


                // If genome archive was not found use first item
                // we have in the list
                if (!wasFound) {
                    int count = genomeComboBox.getItemCount();

                    for (int i = 0; i < count; i++) {
                        Object object = genomeComboBox.getItemAt(i);

                        if (object instanceof GenomeListItem) {
                            GenomeListItem item = (GenomeListItem) object;

                            // We found the genome we want moved it to the local cache
                            // if it's not there already
                            try {
                                GenomeManager.getInstance().loadGenome(item, item.isUserDefined());

                                IGVMainFrame.getInstance().setGenomeId(item.getId());
                                genomeComboBox.setSelectedIndex(i);
                                break;
                            } catch (FileNotFoundException ex) {
                                log.error("Genome file not found ", ex);
                            }
                        }
                    }
                }
                log.debug("Finish selectGenomeFromList");

            }

            public String getName() {
                return "selectGenomeFromList: " + genomeId;  //To change body of implemented methods use File | Settings | File Templates.
            }
        };

        WaitCursorManager.CursorToken token = WaitCursorManager.showWaitCursor();
        try {
            runnable.run();
        }
        finally {
            WaitCursorManager.removeWaitCursor(token);
        }
    }

    /**
     * Method description
     *
     * @return
     */
    public DefaultComboBoxModel getModelForGenomeListComboBox() {

        LinkedHashSet<Object> list = new LinkedHashSet();

        if ((userDefinedGenomeItemList != null) && !userDefinedGenomeItemList.isEmpty()) {
            for (GenomeListItem item : userDefinedGenomeItemList) {
                list.add(item);
            }
            list.add(UIConstants.GENOME_LIST_SEPARATOR);
        }

        if ((serverGenomeItemList != null) && !serverGenomeItemList.isEmpty()) {
            for (GenomeListItem item : this.serverGenomeItemList) {
                list.add(item);
            }

            if ((cachedGenomeItemList == null) || cachedGenomeItemList.isEmpty()) {
                list.add(UIConstants.GENOME_LIST_SEPARATOR);
            }
        }

        if ((cachedGenomeItemList != null) && !cachedGenomeItemList.isEmpty()) {
            for (GenomeListItem item : this.cachedGenomeItemList) {
                list.add(item);
            }

            list.add(UIConstants.GENOME_LIST_SEPARATOR);
        }

        if (list.isEmpty()) {
            list.add(GenomeManager.getInstance().getDefaultGenomeListItem());
        }

        return new DefaultComboBoxModel(list.toArray());
    }

    /**
     * Method description
     *
     * @param monitor
     * @throws FileNotFoundException
     * @throws NoRouteToHostException
     */
    public void intializeDefaultGenome(final ProgressMonitor monitor)
            throws FileNotFoundException, NoRouteToHostException {


        LinkedHashSet<GenomeListItem> mergedServerAndCacheItemList = new LinkedHashSet();

        if (this.userDefinedGenomeItemList != null) {
            mergedServerAndCacheItemList.addAll(this.userDefinedGenomeItemList);
        }

        if (this.cachedGenomeItemList != null) {
            mergedServerAndCacheItemList.addAll(this.cachedGenomeItemList);
        }

        if (this.serverGenomeItemList != null) {
            mergedServerAndCacheItemList.addAll(this.serverGenomeItemList);
        }

        // If the genome has not been set user the one at the top of the list
        /*String genome = GenomeManager.getInstance().getGenomeId();
        if (genome == null) {
        for (GenomeListItem item : mergedServerAndCacheItemList) {

        if (!item.isUserDefined()) {
        GenomeManager.getInstance().loadGenome(item, false, item.isUserDefined());

        String idOfFirstServerGenomeInList = item.getId();

        if (idOfFirstServerGenomeInList == null) {
        return;
        }

        IGVMainFrame.getInstance().addInaccessibleGenomeId(idOfFirstServerGenomeInList);
        selectGenomeFromList(idOfFirstServerGenomeInList);
        IGVMainFrame.getInstance().removeInaccessibleGenomeId(
        idOfFirstServerGenomeInList);
        break;
        }
        }

        }*/
        this.grabFocus();
    }

    /**
     * Method description
     *
     * @param clientItemList
     * @param serverItemList
     * @param cachedGenomeItemList
     */
    public void setGenomeItemList(LinkedHashSet<GenomeListItem> clientItemList,
                                  LinkedHashSet<GenomeListItem> serverItemList,
                                  LinkedHashSet<GenomeListItem> cachedGenomeItemList) {

        if (clientItemList == null) {
            clientItemList = new LinkedHashSet<GenomeListItem>();
        }

        if (serverItemList == null) {
            serverItemList = new LinkedHashSet<GenomeListItem>();
        }

        if (cachedGenomeItemList == null) {
            cachedGenomeItemList = new LinkedHashSet<GenomeListItem>();
        }


        this.userDefinedGenomeItemList = clientItemList;
        this.cachedGenomeItemList = cachedGenomeItemList;
        this.serverGenomeItemList = serverItemList;
    }

    /**
     * Method description
     */
    public void setGenomeItemListModel() {
        genomeComboBox.setModel(getModelForGenomeListComboBox());
    }

    /**
     * This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    private void initComponents() {

        setMinimumSize(new Dimension(200, 32));

        // setPreferredSize(new Dimension(800, 32));

        JideBoxLayout layout = new JideBoxLayout(this, JideBoxLayout.X_AXIS);

        setLayout(layout);

        // This controls the vertical height of the command bar

        locationPanel = new javax.swing.JPanel();
        locationPanel.setBorder(new LineBorder(Color.lightGray, 1, true));

        // BorderFactory.createMatteBorder(2, 2, 2, 2, Color.lightGray));
        // new javax.swing.border.SoftBevelBorder(javax.swing.border.BevelBorder.RAISED));
        locationPanel.setPreferredSize(new java.awt.Dimension(150, 20));
        locationPanel.setLayout(new JideBoxLayout(locationPanel, JideBoxLayout.X_AXIS));
        locationPanel.setAlignmentY(CENTER_ALIGNMENT);
        locationPanel.add(Box.createRigidArea(new Dimension(10, 36)), JideBoxLayout.FIX);

        genomeComboBox = new JComboBox() {

            @Override
            public void setSelectedIndex(int index) {

                Object object = getItemAt(index);
                String text = object.toString();

                if (object instanceof GenomeListItem) {

                    // Only allow selection of genomes
                    super.setSelectedIndex(index);
                } else if (text.equals(UIConstants.LOAD_GENOME_LIST_MENU_ITEM)) {
                    IGVMainFrame.getInstance().doLoadGenome(null);
                }
            }
        };
        genomeComboBox.setMinimumSize(new Dimension(180, 27));
        genomeComboBox.setPreferredSize(new Dimension(180, 27));
        locationPanel.add(genomeComboBox, JideBoxLayout.FIX);
        locationPanel.add(Box.createHorizontalStrut(5), JideBoxLayout.FIX);

        chromosomeComboBox = new javax.swing.JComboBox();
        chromosomeComboBox.setToolTipText("Select a chromosome to view");
        chromosomeComboBox.setMaximumSize(new java.awt.Dimension(DEFAULT_CHROMOSOME_DROPDOWN_WIDTH, 35));
        chromosomeComboBox.setMinimumSize(new java.awt.Dimension(DEFAULT_CHROMOSOME_DROPDOWN_WIDTH, 27));
        chromosomeComboBox.setPreferredSize(new java.awt.Dimension(DEFAULT_CHROMOSOME_DROPDOWN_WIDTH, 16));
        chromosomeComboBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                chromosomeComboBoxActionPerformed(evt);
            }
        });
        locationPanel.add(chromosomeComboBox, JideBoxLayout.FIX);
        locationPanel.add(Box.createHorizontalStrut(5), JideBoxLayout.FIX);

        searchTextField = new SearchTextField();
        searchTextField.setToolTipText("Enter a gene of locus, e.f. EGFR,   chr1,   or chr1:100,000-200,000");
        searchTextField.setMaximumSize(new java.awt.Dimension(250, 15));
        searchTextField.setMinimumSize(new java.awt.Dimension(100, 28));
        searchTextField.setPreferredSize(new java.awt.Dimension(230, 28));
        searchTextField.setAlignmentY(CENTER_ALIGNMENT);
        locationPanel.add(searchTextField, JideBoxLayout.FIX);

        goButton = new JideButton("Go");
        // goButton.setButtonStyle(ButtonStyle.TOOLBOX_STYLE);

        // goButton.setPreferredSize(new java.awt.Dimension(30, 30));
        // goButton.setMaximumSize(new java.awt.Dimension(30, 30));
        // goButton.setMinimumSize(new java.awt.Dimension(30, 30));
        // goButton.setText("Go");
        goButton.setToolTipText("Jump to gene or locus");
        goButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                goButtonActionPerformed(evt);
            }
        });
        locationPanel.add(goButton, JideBoxLayout.FIX);

        add(locationPanel, JideBoxLayout.FIX);

        add(Box.createHorizontalStrut(10), JideBoxLayout.FIX);

        toolPanel = new javax.swing.JPanel();
        toolPanel.setAlignmentX(RIGHT_ALIGNMENT);
        toolPanel.setLayout(new JideBoxLayout(toolPanel, JideBoxLayout.X_AXIS));

        homeButton = new com.jidesoft.swing.JideButton();
        homeButton.setAlignmentX(RIGHT_ALIGNMENT);
        homeButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        homeButton.setIcon(new javax.swing.ImageIcon(
                getClass().getResource("/toolbarButtonGraphics/navigation/Home24.gif")));
        homeButton.setMaximumSize(new java.awt.Dimension(32, 32));
        homeButton.setMinimumSize(new java.awt.Dimension(32, 32));
        homeButton.setPreferredSize(new java.awt.Dimension(32, 32));
        homeButton.setToolTipText("Jump to whole genome view");
        homeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                homeButtonActionPerformed(evt);
            }
        });
        toolPanel.add(homeButton, JideBoxLayout.FIX);


        // toolPanel.setBorder(
        // new javax.swing.border.SoftBevelBorder(javax.swing.border.BevelBorder.RAISED));
        backButton = new JideButton();
        backButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        backButton.setIcon(new javax.swing.ImageIcon(
                getClass().getResource("/images/left-arrow.gif")));
        backButton.setToolTipText("Go back");
        backButton.setMaximumSize(new java.awt.Dimension(32, 32));
        backButton.setMinimumSize(new java.awt.Dimension(32, 32));
        backButton.setPreferredSize(new java.awt.Dimension(32, 32));
        backButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                IGVMainFrame.getInstance().getSession().getHistory().back();

            }
        });
        backButton.setEnabled(false);
        toolPanel.add(backButton, JideBoxLayout.FIX);

        forwardButton = new JideButton();
        forwardButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        forwardButton.setIcon(new javax.swing.ImageIcon(
                getClass().getResource("/images/right-arrow.gif")));    // NOI18N
        forwardButton.setToolTipText("Go forward");
        forwardButton.setMaximumSize(new java.awt.Dimension(32, 32));
        forwardButton.setMinimumSize(new java.awt.Dimension(32, 32));
        forwardButton.setPreferredSize(new java.awt.Dimension(32, 32));
        forwardButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                IGVMainFrame.getInstance().getSession().getHistory().forward();
            }
        });
        forwardButton.setEnabled(false);
        toolPanel.add(forwardButton, JideBoxLayout.FIX);

        refreshButton = new com.jidesoft.swing.JideButton();
        refreshButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        refreshButton.setAlignmentX(RIGHT_ALIGNMENT);
        refreshButton.setIcon(new javax.swing.ImageIcon(
                getClass().getResource("/toolbarButtonGraphics/general/Refresh24.gif")));    // NOI18N
        refreshButton.setMaximumSize(new java.awt.Dimension(32, 32));
        refreshButton.setMinimumSize(new java.awt.Dimension(32, 32));
        refreshButton.setPreferredSize(new java.awt.Dimension(32, 32));
        refreshButton.setToolTipText("Refresh the screen");
        refreshButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                refreshButtonActionPerformed(evt);
            }
        });
        toolPanel.add(refreshButton, JideBoxLayout.FIX);


        Icon regionOfInterestIcon =
                IconFactory.getInstance().getIcon(IconFactory.IconID.REGION_OF_INTEREST);

        roiToggleButton = new JideToggleButton(regionOfInterestIcon);
        roiToggleButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        roiToggleButton.setAlignmentX(RIGHT_ALIGNMENT);
        roiToggleButton.setToolTipText("Define a region of interest.");
        roiToggleButton.setMaximumSize(new java.awt.Dimension(32, 32));
        roiToggleButton.setMinimumSize(new java.awt.Dimension(32, 32));
        roiToggleButton.setPreferredSize(new java.awt.Dimension(32, 32));
        roiToggleButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                roiToggleButtonActionPerformed(evt);
            }
        });
        toolPanel.add(roiToggleButton, JideBoxLayout.FIX);

        Icon fitToWindowIcon =
                IconFactory.getInstance().getIcon(IconFactory.IconID.REGION_OF_INTEREST);

        fitToWindowButton = new JideButton();
        fitToWindowButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        fitToWindowButton.setAlignmentX(RIGHT_ALIGNMENT);
        fitToWindowButton.setIcon(new javax.swing.ImageIcon(getClass().getResource("/images/collapseall.gif")));
        fitToWindowButton.setMaximumSize(new java.awt.Dimension(32, 32));
        fitToWindowButton.setMinimumSize(new java.awt.Dimension(32, 32));
        fitToWindowButton.setPreferredSize(new java.awt.Dimension(32, 32));
        fitToWindowButton.setToolTipText("Squish tracks to fit in window.");
        fitToWindowButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                (new FitDataToWindowMenuAction(null, 0, owner)).actionPerformed(evt);
            }
        });
        toolPanel.add(fitToWindowButton, JideBoxLayout.FIX);


        this.add(toolPanel);

        this.add(Box.createHorizontalGlue(), JideBoxLayout.VARY);

        zoomControl = new ZoomSliderPanel();

        // zoomControl.setAlignmentX(RIGHT_ALIGNMENT);
        zoomControl.setPreferredSize(new Dimension(200, 30));
        zoomControl.setMinimumSize(new Dimension(200, 30));
        zoomControl.setMaximumSize(new Dimension(200, 30));
        zoomControl.setToolTipText("Click + to zoom in,  - to zoom out.");
        zoomControl.setOpaque(false);
        this.add(zoomControl, JideBoxLayout.FIX);

        this.add(Box.createHorizontalStrut(20), JideBoxLayout.FIX);
    }

    /**
     * Method description
     *
     * @return
     */
    public GenomeListItem getGenomeSelectedInDropdown() {
        return (GenomeListItem) genomeComboBox.getSelectedItem();
    }

    private void adjustChromosomeDropdownWidth(int width) {

        int newWidth = (width > DEFAULT_CHROMOSOME_DROPDOWN_WIDTH)
                ? width : DEFAULT_CHROMOSOME_DROPDOWN_WIDTH;

        chromosomeComboBox.setMaximumSize(new java.awt.Dimension(newWidth, 35));
        chromosomeComboBox.setMinimumSize(new java.awt.Dimension(newWidth, 27));
        chromosomeComboBox.setPreferredSize(new java.awt.Dimension(newWidth, 16));
        revalidate();
    }

    private void homeButtonActionPerformed(java.awt.event.ActionEvent evt) {
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (FrameManager.isGeneListMode()) {
            owner.setGeneList(null);
        }
        if (genome != null) {
            String chr = genome.getHomeChromosome();
            getDefaultReferenceFrame().setChromosomeName(chr);
            IGVMainFrame.getInstance().getSession().getHistory().push(chr, getDefaultReferenceFrame().getZoom());
            chromosomeComboBox.setSelectedItem(chr);
            updateCurrentCoordinates();
            owner.chromosomeChangeEvent();
            owner.repaint();
        }
    }

    private void refreshButtonActionPerformed(java.awt.event.ActionEvent evt) {
        //LRUCache.clearCaches();
        owner.doRefresh();
    }

    private void chromosomeComboBoxActionPerformed(java.awt.event.ActionEvent evt) {
        JComboBox combobox = (JComboBox) evt.getSource();
        String chromosomeName = (String) combobox.getSelectedItem();
        if (chromosomeName != null) {

            if (!chromosomeName.equals(getDefaultReferenceFrame().getChrName())) {
                getDefaultReferenceFrame().setChromosomeName(chromosomeName);
                getDefaultReferenceFrame().recordHistory();
                updateCurrentCoordinates();
                owner.chromosomeChangeEvent();
                owner.repaint();
                PreferenceManager.getInstance().setLastChromosomeViewed(chromosomeName);
            }
        }
    }

    private void goButtonActionPerformed(java.awt.event.ActionEvent evt) {    // GEN-FIRST:event_goButtonActionPerformed
        String searchText = searchTextField.getText();
        searchByLocus(searchText);
    }


    public void searchByLocus(String searchText) {


        if (log.isDebugEnabled()) {
            log.debug("Enter search by locus: " + searchText);
        }

        if ((searchText != null) && (searchText.length() > 0)) {
            searchTextField.setText(searchText);
            (new SearchCommand(getDefaultReferenceFrame(), searchText)).execute();
            chromosomeComboBox.setSelectedItem(getDefaultReferenceFrame().getChrName());
        }

        if (log.isDebugEnabled()) {
            log.debug("Exit search by locus: " + searchText);
        }
    }


    private void roiToggleButtonActionPerformed(java.awt.event.ActionEvent evt) {    // GEN-FIRST:event_roiToggleButtonActionPerformed
        if (roiToggleButton.isSelected()) {
            owner.beginROI(roiToggleButton);
        } else {
            owner.endROI();
        }
    }

    static class SearchTextField extends javax.swing.JTextField {

        /**
         * Method description
         *
         * @param arg0
         * @param arg1
         * @param arg2
         * @param arg3
         */
        @Override
        public void reshape(int arg0, int arg1, int arg2, int arg3) {
            super.reshape(arg0, arg1, arg2, getPreferredSize().height);
        }

    }


}
