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
 * IGVCommandBar.java
 *
 * Created on April 5, 2008, 10:02 AM
 */
package org.broad.igv.ui.commandbar;

import com.jidesoft.swing.JideBoxLayout;
import com.jidesoft.swing.JideButton;
import com.jidesoft.swing.JideToggleButton;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.event.*;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.GenomeListItem;
import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.prefs.Constants;
import org.broad.igv.prefs.PreferencesManager;
import org.broad.igv.session.History;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.ShowDetailsBehavior;
import org.broad.igv.ui.UIConstants;
import org.broad.igv.ui.action.FitDataToWindowMenuAction;
import org.broad.igv.ui.action.ReloadTracksMenuAction;
import org.broad.igv.ui.panel.FrameManager;
import org.broad.igv.ui.panel.IGVPopupMenu;
import org.broad.igv.ui.panel.ReferenceFrame;
import org.broad.igv.ui.panel.ZoomSliderPanel;
import org.broad.igv.ui.util.IconFactory;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.ui.util.UIUtilities;

import javax.swing.*;
import javax.swing.border.LineBorder;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;

/**
 * @author jrobinso
 */
public class IGVCommandBar extends javax.swing.JPanel implements IGVEventObserver {

    private static Logger log = Logger.getLogger(IGVCommandBar.class);

    final static String MODIFY_DETAILS_TOOLTIP = "Modify popup text behavior in data panels";
    final static int DEFAULT_CHROMOSOME_DROPDOWN_WIDTH = 120;

    private ChromosomeComboBox chromosomeComboBox;
    private GenomeComboBox genomeComboBox;
    private JideButton goButton;
    private JideButton homeButton;
    private JPanel locationPanel;
    private JideButton refreshButton;
    private JideToggleButton roiToggleButton;
    private JideButton detailsBehaviorButton;
    private JideToggleButton rulerLineButton;
    private SearchTextField searchTextField;
    private JPanel toolPanel;
    private JPanel zoomControl;

    private JideButton backButton;
    private JideButton forwardButton;
    private JideButton fitToWindowButton;

    private ShowDetailsBehavior detailsBehavior;


    public IGVCommandBar() {

        initComponents();

        // Post creation widget setup.
        refreshGenomeListComboBox();

        String currentChr = FrameManager.getDefaultFrame().getChrName();
        boolean isWholeGenome = currentChr.equals(Globals.CHR_ALL);

        chromosomeComboBox.setSelectedItem(currentChr);
        roiToggleButton.setEnabled(!isWholeGenome);
        zoomControl.setEnabled(!isWholeGenome);

        detailsBehaviorButton.addMouseListener(new MouseAdapter() {
            public void mousePressed(MouseEvent e) {
                getPopupMenuToolTipBehavior().show(e.getComponent(), e.getX(), e.getY());
            }
        });

        IGVEventBus.getInstance().subscribe(ViewChange.class, this);
        IGVEventBus.getInstance().subscribe(GenomeChangeEvent.class, this);
        IGVEventBus.getInstance().subscribe(GenomeResetEvent.class, this);

    }


    private JPopupMenu getPopupMenuToolTipBehavior() {
        final JPopupMenu popup = new IGVPopupMenu();
        for (final ShowDetailsBehavior behavior : ShowDetailsBehavior.values()) {
            JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem(behavior.getLabel());
            menuItem.setSelected(detailsBehavior == behavior);
            menuItem.addActionListener(new AbstractAction() {
                public void actionPerformed(ActionEvent e) {
                    detailsBehavior = behavior;
                    PreferencesManager.getPreferences().put(Constants.DETAILS_BEHAVIOR_KEY, behavior.name());
                }
            });
            popup.add(menuItem);
        }
        return popup;
    }

    public ShowDetailsBehavior getDetailsBehavior() {
        return detailsBehavior;
    }


    public void setGeneListMode(boolean geneListMode) {

        genomeComboBox.setEnabled(!geneListMode);
//        locationPanel.setEnabled(!geneListMode);
        chromosomeComboBox.setEnabled(!geneListMode);
//        searchTextField.setEnabled(!geneListMode);
//        goButton.setEnabled(!geneListMode);
        zoomControl.setEnabled(!geneListMode);
//        homeButton.setEnabled(true);
//        roiToggleButton.setEnabled(!geneListMode);
    }


    /**
     * Selects the first genome from the list which matches this genomeId.
     * If not found, checks genomes from the server/user-defined list
     *
     * @param genomeId
     */
    public void selectGenome(String genomeId) {

        log.info("Selecting genome " + genomeId);

        GenomeListItem selectedItem = GenomeListManager.getInstance().getGenomeListItem(genomeId);

        if (selectedItem == null || !genomeComboBox.hasItem(selectedItem)) {
            try {
                GenomeManager.getInstance().loadGenomeById(genomeId);
            } catch (IOException e) {
                MessageUtils.showErrorMessage("Error loading genome: " + genomeId, e);
                log.error("Error loading genome: " + genomeId, e);
            }
        }

        if (selectedItem != null) {
            UIUtilities.invokeAndWaitOnEventThread(() -> genomeComboBox.setSelectedItem(selectedItem));
        }
    }

    public void updateCurrentCoordinates() {

        if(IGV.hasInstance()) {
            String p = "";
            ReferenceFrame defaultFrame = FrameManager.getDefaultFrame();
            final String chrName = defaultFrame.getChrName();
            if (!Globals.CHR_ALL.equals(chrName) && !FrameManager.isGeneListMode()) {
                p = defaultFrame.getFormattedLocusString();
            }
            final String position = p;
            final History history = IGV.getInstance().getSession().getHistory();

            UIUtilities.invokeOnEventThread(new Runnable() {
                public void run() {
                    searchTextField.setText(position);
                    forwardButton.setEnabled(history.canGoForward());
                    backButton.setEnabled(history.canGoBack());
                    roiToggleButton.setEnabled(!Globals.CHR_ALL.equals(chrName));
                    zoomControl.setEnabled(!Globals.CHR_ALL.equals(chrName));
                }
            });
        }
    }


    public void refreshGenomeListComboBox() {
        UIUtilities.invokeAndWaitOnEventThread(() -> {
            genomeComboBox.refreshGenomeListComboBox();
        });
    }


    /**
     * Adjust the popup for the combobox to be at least as wide as
     * the widest item.
     *
     * @param box
     */
    private void adjustPopupWidth(JComboBox box) {
        if (box.getItemCount() == 0) return;
        Object comp = box.getUI().getAccessibleChild(box, 0);
        if (!(comp instanceof JPopupMenu)) {
            return;
        }
        JPopupMenu popup = (JPopupMenu) comp;
        JScrollPane scrollPane = null;
        for (Component scomp : popup.getComponents()) {
            if (scomp instanceof JScrollPane) {
                scrollPane = (JScrollPane) scomp;
            }
        }
        if (scrollPane == null) return;

        //Loop through and set width to widest component, plus some padding
        int rendererWidth = box.getWidth();
        for (int index = 0; index < box.getItemCount(); index++) {
            Object value = box.getItemAt(index);
            Component rendererComp = box.getRenderer().getListCellRendererComponent(null, value, index, false, false);
        }

        Dimension size = scrollPane.getPreferredSize();
        size.width = Math.max(size.width, rendererWidth);
        scrollPane.setPreferredSize(size);
        scrollPane.setMaximumSize(size);
        scrollPane.revalidate();
    }

    //<editor-fold desc="Action methods">
    private void homeButtonActionPerformed(java.awt.event.ActionEvent evt) {
        Genome genome = GenomeManager.getInstance().getCurrentGenome();
        if (FrameManager.isGeneListMode()) {
            IGV.getInstance().setGeneList(null);
        }
        if (genome != null) {
            String chrName = genome.getHomeChromosome();
            if (chrName != null && !chrName.equals(chromosomeComboBox.getSelectedItem())) {
                FrameManager.getDefaultFrame().changeChromosome(chrName, false);
            }
        }
    }

    private void refreshButtonActionPerformed(java.awt.event.ActionEvent evt) {

        IGVEventBus.getInstance().post(new org.broad.igv.event.RefreshEvent());
        (new ReloadTracksMenuAction("",-1, IGV.getInstance())).actionPerformed(evt);

    }

    private void goButtonActionPerformed(java.awt.event.ActionEvent evt) {    // GEN-FIRST:event_goButtonActionPerformed
        String searchText = searchTextField.getText();
        searchByLocus(searchText);
    }

    private void roiToggleButtonActionPerformed(java.awt.event.ActionEvent evt) {    // GEN-FIRST:event_roiToggleButtonActionPerformed
        if (roiToggleButton.isSelected()) {
            IGV.getInstance().beginROI(roiToggleButton);
        } else {
            IGV.getInstance().endROI();
        }
    }

    //</editor-fold>

    public void receiveEvent(Object e) {

        if (e instanceof ViewChange) {
            ViewChange event = (ViewChange) e;
            if (event.type == ViewChange.Type.ChromosomeChange || event.type == ViewChange.Type.LocusChange) {
                String chrName = FrameManager.getDefaultFrame().getChrName();
                roiToggleButton.setEnabled(!Globals.CHR_ALL.equals(chrName));
                zoomControl.setEnabled(!Globals.CHR_ALL.equals(chrName));
                if (!chrName.equals(chromosomeComboBox.getSelectedItem())) {
                    chromosomeComboBox.setSelectedItem(chrName);
                }
            }

            updateCurrentCoordinates();
            repaint(); // TODO Is this neccessary?
        } else if (e instanceof GenomeChangeEvent) {
            GenomeChangeEvent event = (GenomeChangeEvent) e;
            Genome genome = event.genome;
            refreshGenomeListComboBox();
            chromosomeComboBox.updateChromosFromGenome(genome);
            String chrName = FrameManager.getDefaultFrame().getChrName();
            zoomControl.setEnabled(!Globals.CHR_ALL.equals(chrName));
        } else if (e instanceof GenomeResetEvent) {
            refreshGenomeListComboBox();
        } else {
            log.info("Unknown event class: " + e.getClass());
        }
    }


    //<editor-fold desc="Search box">
    // Set the focus in the search box

    public void focusSearchBox() {
        searchTextField.requestFocusInWindow();
        searchTextField.selectAll();
    }


    public void searchByLocus(final String searchText) {

        if ((searchText != null) && (searchText.length() > 0)) {
            String homeChr = IGV.getInstance().getGenomeManager().getCurrentGenome().getHomeChromosome();
            if (searchText.equalsIgnoreCase("home") || searchText.equalsIgnoreCase(homeChr)) {
                homeButtonActionPerformed(null);
            } else {
                searchTextField.setText(searchText);
                searchTextField.searchByLocus(searchText);
            }
        }
    }


    /**
     * This method is called from within the constructor
     */
    private void initComponents() {

        setMinimumSize(new Dimension(200, 32));

        JideBoxLayout layout = new JideBoxLayout(this, JideBoxLayout.X_AXIS);

        setLayout(layout);

        final String detailsPreference = PreferencesManager.getPreferences().get(Constants.DETAILS_BEHAVIOR_KEY);
        detailsBehavior = ShowDetailsBehavior.valueOf((detailsPreference.toUpperCase()));

        // This controls the vertical height of the command bar

        locationPanel = new javax.swing.JPanel();
        locationPanel.setBorder(new LineBorder(Color.lightGray, 1, true));

        // BorderFactory.createMatteBorder(2, 2, 2, 2, Color.lightGray));
        // new javax.swing.border.SoftBevelBorder(javax.swing.border.BevelBorder.RAISED));
        locationPanel.setPreferredSize(new java.awt.Dimension(150, 20));
        locationPanel.setLayout(new JideBoxLayout(locationPanel, JideBoxLayout.X_AXIS));
        locationPanel.setAlignmentY(CENTER_ALIGNMENT);
        locationPanel.add(Box.createRigidArea(new Dimension(10, 36)), JideBoxLayout.FIX);

        genomeComboBox = new GenomeComboBox();
        genomeComboBox.setMinimumSize(new Dimension(180, 27));
        genomeComboBox.setPreferredSize(new Dimension(180, 27));
        genomeComboBox.setToolTipText(UIConstants.CHANGE_GENOME_TOOLTIP);

        genomeComboBox.addPopupMenuListener(new PopupMenuListener() {
            @Override
            public void popupMenuWillBecomeVisible(PopupMenuEvent e) {
                try {
                    adjustPopupWidth(genomeComboBox);
                } catch (Exception e1) {
                    log.warn(e1.getMessage(), e1);
                }
            }

            @Override
            public void popupMenuWillBecomeInvisible(PopupMenuEvent e) {
                //TODO
            }

            @Override
            public void popupMenuCanceled(PopupMenuEvent e) {
                //TODO
            }
        });

        locationPanel.add(genomeComboBox, JideBoxLayout.FIX);
        locationPanel.add(Box.createHorizontalStrut(5), JideBoxLayout.FIX);


        chromosomeComboBox = new ChromosomeComboBox();
        chromosomeComboBox.setToolTipText("Select a chromosome to view");
        chromosomeComboBox.setMaximumSize(new java.awt.Dimension(DEFAULT_CHROMOSOME_DROPDOWN_WIDTH, 30));
        chromosomeComboBox.setMinimumSize(new java.awt.Dimension(DEFAULT_CHROMOSOME_DROPDOWN_WIDTH, 30));
        chromosomeComboBox.setPreferredSize(new java.awt.Dimension(DEFAULT_CHROMOSOME_DROPDOWN_WIDTH, 30));

        locationPanel.add(chromosomeComboBox, JideBoxLayout.FIX);
        locationPanel.add(Box.createHorizontalStrut(5), JideBoxLayout.FIX);

        searchTextField = new SearchTextField();
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
        //final Border toolButtonBorder = BorderFactory.createLineBorder(Color.gray, 1);

        homeButton = new com.jidesoft.swing.JideButton();
        homeButton.setAlignmentX(RIGHT_ALIGNMENT);
        //homeButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        // homeButton.setBorder(toolButtonBorder);
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
        //backButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        //backButton.setBorder(toolButtonBorder);
        backButton.setIcon(new javax.swing.ImageIcon(getClass().getResource("/images/left-arrow.gif")));
        backButton.setToolTipText("Go back");
        backButton.setMaximumSize(new java.awt.Dimension(32, 32));
        backButton.setMinimumSize(new java.awt.Dimension(32, 32));
        backButton.setPreferredSize(new java.awt.Dimension(32, 32));
        backButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                IGV.getInstance().getSession().getHistory().back();

            }
        });
        backButton.setEnabled(false);
        toolPanel.add(backButton, JideBoxLayout.FIX);

        forwardButton = new JideButton();
        //forwardButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        //forwardButton.setBorder(toolButtonBorder);
        forwardButton.setIcon(new javax.swing.ImageIcon(getClass().getResource("/images/right-arrow.gif")));
        forwardButton.setToolTipText("Go forward");
        forwardButton.setMaximumSize(new java.awt.Dimension(32, 32));
        forwardButton.setMinimumSize(new java.awt.Dimension(32, 32));
        forwardButton.setPreferredSize(new java.awt.Dimension(32, 32));
        forwardButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                IGV.getInstance().getSession().getHistory().forward();
            }
        });
        forwardButton.setEnabled(false);
        toolPanel.add(forwardButton, JideBoxLayout.FIX);

        refreshButton = new com.jidesoft.swing.JideButton();
        //refreshButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        //refreshButton.setBorder(toolButtonBorder);
        refreshButton.setAlignmentX(RIGHT_ALIGNMENT);
        refreshButton.setIcon(new javax.swing.ImageIcon(
                getClass().getResource("/toolbarButtonGraphics/general/Refresh24.gif")));    // NOI18N
        refreshButton.setMaximumSize(new java.awt.Dimension(32, 32));
        refreshButton.setMinimumSize(new java.awt.Dimension(32, 32));
        refreshButton.setPreferredSize(new java.awt.Dimension(32, 32));
        refreshButton.setToolTipText("Reload tracks and refresh the screen");
        refreshButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                refreshButtonActionPerformed(evt);
            }
        });
        toolPanel.add(refreshButton, JideBoxLayout.FIX);


        Icon regionOfInterestIcon =
                IconFactory.getInstance().getIcon(IconFactory.IconID.REGION_OF_INTEREST);

        roiToggleButton = new JideToggleButton(regionOfInterestIcon);
        //roiToggleButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        //roiToggleButton.setBorder(toolButtonBorder);
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


        fitToWindowButton = new JideButton();
        //fitToWindowButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        //fitToWindowButton.setBorder(toolButtonBorder);
        fitToWindowButton.setAlignmentX(RIGHT_ALIGNMENT);
        fitToWindowButton.setIcon(new javax.swing.ImageIcon(getClass().getResource("/images/collapseall.gif")));
        fitToWindowButton.setMaximumSize(new java.awt.Dimension(32, 32));
        fitToWindowButton.setMinimumSize(new java.awt.Dimension(32, 32));
        fitToWindowButton.setPreferredSize(new java.awt.Dimension(32, 32));
        fitToWindowButton.setToolTipText("Resize tracks to fit in window.");
        fitToWindowButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(java.awt.event.ActionEvent evt) {
                (new FitDataToWindowMenuAction(null, 0, IGV.getInstance())).actionPerformed(evt);
            }
        });
        toolPanel.add(fitToWindowButton, JideBoxLayout.FIX);

        final Icon noTooltipIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.NO_TOOLTIP);
        final Icon tooltipIcon = IconFactory.getInstance().getIcon(IconFactory.IconID.TOOLTIP);
        detailsBehaviorButton = new JideButton(noTooltipIcon);

        //detailsBehaviorButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        //detailsBehaviorButton.setBorder(toolButtonBorder);
        detailsBehaviorButton.setAlignmentX(RIGHT_ALIGNMENT);
        detailsBehaviorButton.setToolTipText(MODIFY_DETAILS_TOOLTIP);
        detailsBehaviorButton.setMaximumSize(new java.awt.Dimension(32, 32));
        detailsBehaviorButton.setMinimumSize(new java.awt.Dimension(32, 32));
        detailsBehaviorButton.setPreferredSize(new java.awt.Dimension(32, 32));
        toolPanel.add(detailsBehaviorButton, JideBoxLayout.FIX);

        rulerLineButton = new JideToggleButton();
        //roiToggleButton.setButtonStyle(JideButton.TOOLBOX_STYLE);
        //roiToggleButton.setBorder(toolButtonBorder);
        rulerLineButton.setIcon(new javax.swing.ImageIcon(getClass().getResource("/images/vertical-line.gif")));
        rulerLineButton.setAlignmentX(RIGHT_ALIGNMENT);
        rulerLineButton.setToolTipText("Enable ruler line in data panels");
        rulerLineButton.setMaximumSize(new java.awt.Dimension(32, 32));
        rulerLineButton.setMinimumSize(new java.awt.Dimension(32, 32));
        rulerLineButton.setPreferredSize(new java.awt.Dimension(32, 32));
        rulerLineButton.addActionListener(evt -> {
            IGV.getInstance().setRulerEnabled(rulerLineButton.isSelected());
            IGV.getInstance().repaintContentPane();
        });
        toolPanel.add(rulerLineButton, JideBoxLayout.FIX);

        this.add(toolPanel);

        this.add(Box.createHorizontalGlue(), JideBoxLayout.VARY);

        zoomControl = new ZoomSliderPanel();

        // zoomControl.setAlignmentX(RIGHT_ALIGNMENT);
        Dimension dimSize = new Dimension(200, 30);
        zoomControl.setPreferredSize(dimSize);
        zoomControl.setMinimumSize(dimSize);
        zoomControl.setMaximumSize(dimSize);
        zoomControl.setToolTipText("Click + to zoom in,  - to zoom out");
        zoomControl.setOpaque(false);
        this.add(zoomControl, JideBoxLayout.FIX);

        this.add(Box.createHorizontalStrut(20), JideBoxLayout.FIX);


    }

}
