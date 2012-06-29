/*
 * Created by JFormDesigner on Thu Jun 14 08:42:46 EDT 2012
 */

package org.broad.igv.ui.panel;

import java.awt.event.*;

import com.jidesoft.swing.*;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.track.Track;
import org.broad.igv.ui.IGV;

import java.awt.*;
import java.util.ArrayList;
import java.util.Collection;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.ListDataListener;

/**
 * @author Stan Diamond
 */


public class FeatureTrackSelectionDialog extends JDialog {

    boolean isCanceled = false;

    public FeatureTrackSelectionDialog(Frame owner) {
        super(owner);
        initComponents();
        setModal(true);

        Collection<Track> tracks = IGV.getInstance().getAllTracks();
        ArrayList<TrackWrapper> wrappers = new ArrayList<TrackWrapper>();
        ArrayList<TrackWrapper> selectedObjects = new ArrayList<TrackWrapper>();
        for (Track t : tracks) {
            if (t instanceof FeatureTrack) {
                TrackWrapper trackWrapper = new TrackWrapper((FeatureTrack) t);
                wrappers.add(trackWrapper);
            }
        }
        featureTrackList.setListData(wrappers.toArray());
        featureTrackList.setSelectedObjects(selectedObjects.toArray());

    }

    public FeatureTrackSelectionDialog(Dialog owner) {
        super(owner);
        initComponents();
    }

    private void okButtonActionPerformed(ActionEvent e) {
        setVisible(false);
        dispose();
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        isCanceled = true;
        setVisible(false);
        dispose();
    }

    public java.util.List<FeatureTrack> getSelectedTracks() {
        if (isCanceled) return null;
        Object[] selections = featureTrackList.getSelectedValues();
        ArrayList<FeatureTrack> tracks = new ArrayList<FeatureTrack>(selections.length);
        for (Object wrapper : selections) {
            tracks.add(((TrackWrapper) wrapper).track);
        }
        return tracks;
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        featureTrackPanel = new JScrollPane();
        featureTrackList = new CheckBoxList();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.X_AXIS));

                //======== featureTrackPanel ========
                {
                    featureTrackPanel.setViewportView(featureTrackList);
                }
                contentPanel.add(featureTrackPanel);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout) buttonBar.getLayout()).columnWidths = new int[]{0, 85, 80};
                ((GridBagLayout) buttonBar.getLayout()).columnWeights = new double[]{1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("OK");
                okButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
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
        }
        contentPane.add(dialogPane, BorderLayout.SOUTH);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JScrollPane featureTrackPanel;
    private CheckBoxList featureTrackList;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


    static class TrackWrapper {
        FeatureTrack track;

        TrackWrapper(FeatureTrack track) {
            this.track = track;
        }

        public String toString() {
            return track.getName();
        }
    }
}
