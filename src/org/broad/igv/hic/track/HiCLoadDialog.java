/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not
 * responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which is
 * available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.hic.track;

import org.broad.igv.track.Track;
import org.broad.igv.util.ResourceLocator;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;

/**
 * @author Jim Robinson
 * @date 5/8/12
 */
public class HiCLoadDialog extends JDialog {

    private boolean canceled = false;
    private Collection<String> selectedTracks = new HashSet<String>();

    public static void main(String[] args) {
        Map<String, List<ResourceLocator>> locators = HiCTrackManager.getTrackLocators();
        List<HiCTrack> tracks = HiCTrackManager.getLoadedTracks();
        HiCLoadDialog dlg = new HiCLoadDialog(null, locators, tracks);
        dlg.setVisible(true);
        System.exit(1);
    }

    public HiCLoadDialog(Frame parent, Map<String, List<ResourceLocator>> locators, List<HiCTrack> tracks) {
        super(parent);
        initComponents(locators, tracks);
        setModal(true);
        this.setSize(750, 800);
    }

    private void initComponents(Map<String, List<ResourceLocator>> locators, List<HiCTrack> tracks) {

        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());


        //======== dialogPane ========
        JPanel dialogPane = new JPanel();
        dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
        dialogPane.setLayout(new BorderLayout());

        final Box mainPanel = new Box(BoxLayout.Y_AXIS);
        mainPanel.setAlignmentX(LEFT_ALIGNMENT);

        Set<String> loadedTrackNames = new HashSet<String>(tracks.size());
        for (HiCTrack t : tracks) {
            loadedTrackNames.add(t.getName());
        }

        for (Map.Entry<String, List<ResourceLocator>> entry : locators.entrySet()) {
            String catName = entry.getKey();
            List<ResourceLocator> locatorList = entry.getValue();
            mainPanel.add(new CategoryPanel(catName, locatorList, loadedTrackNames));
        }

        JScrollPane sp = new JScrollPane(mainPanel);
        sp.setBackground(mainPanel.getBackground());
        dialogPane.add(sp, BorderLayout.CENTER);
        contentPane.add(dialogPane, BorderLayout.CENTER);


        JPanel buttonBar = new JPanel();
        buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
        buttonBar.setLayout(new GridBagLayout());
        ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
        ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

        //---- okButton ----
        JButton okButton = new JButton("OK");
        okButton.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                canceled = false;
                for (Component c : mainPanel.getComponents()) {
                    if (c instanceof CategoryPanel) {
                        selectedTracks.addAll(((CategoryPanel) c).getSelectedTracks());
                    }
                }
                setVisible(false);
            }
        });


        //---- cancelButton ----
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                canceled = true;
                setVisible(false);
            }
        });

        buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 0), 0, 0));
        buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                new Insets(0, 0, 0, 5), 0, 0));


        contentPane.add(buttonBar, BorderLayout.SOUTH);
        pack();
        setLocationRelativeTo(getOwner());
    }

    public boolean isCanceled() {
        return canceled;
    }

    public Collection<String> getSelectedTracks() {
        return selectedTracks;
    }
}
