/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

/*
 * Created by JFormDesigner on Thu Jun 14 08:42:46 EDT 2012
 */

package org.broad.igv.ui.panel;

import org.broad.igv.track.FeatureTrack;
import org.broad.igv.ui.IGV;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;
import java.util.List;

/**
 * @author Stan Diamond
 */


public class FeatureTrackSelectionDialog extends JDialog {

    boolean isCanceled = false;

    public FeatureTrackSelectionDialog(Frame owner) {
        super(owner);
        initComponents();
        setModal(true);

        Collection<FeatureTrack> tracks = IGV.getInstance().getFeatureTracks();
        List<TrackWrapper> wrappers = TrackWrapper.wrapTracks(tracks);
        featureTrackList.setListData(wrappers.toArray());

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

    public FeatureTrack getSelectedTrack() {
        if (isCanceled) return null;
        Object selection = featureTrackList.getCheckBoxListSelectedValue();
        if (selection == null) return null;
        return ((TrackWrapper) selection).getTrack();
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        vSpacer1 = new JPanel(null);
        contentPanel = new JPanel();
        featureTrackPanel = new JScrollPane();
        featureTrackList = new RadioButtonSelectionList();
        buttonBar = new JPanel();
        okButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setMinimumSize(new Dimension(220, 220));
            dialogPane.setLayout(new BorderLayout());

            //---- vSpacer1 ----
            vSpacer1.setPreferredSize(new Dimension(10, 30));
            dialogPane.add(vSpacer1, BorderLayout.NORTH);

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.X_AXIS));

                //======== featureTrackPanel ========
                {

                    //---- featureTrackList ----
                    featureTrackList.setClickInCheckBoxOnly(false);
                    featureTrackList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
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
    private JPanel vSpacer1;
    private JPanel contentPanel;
    private JScrollPane featureTrackPanel;
    private RadioButtonSelectionList featureTrackList;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables


//    /**
//     * Just like a CheckBoxList, but renders with radio buttons
//     * and only allows single selection.
//     * User: jacob
//     * Date: 2012-Jul-17
//     */
//    public static class RadioList extends CheckBoxList {
//
//        public RadioList() {
//            super();
//        }
//
//        @Override
//        protected void init() {
//            super.init();
//            getCheckBoxListSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
//        }
//
//        @Override
//        protected CheckBoxListCellRenderer createCellRenderer() {
//            return new RadioListCellRenderer();
//        }
//
//        @Override
//        protected Handler createHandler() {
//            return new SingleSelectionHandler(this);
//        }
//
//
//        protected static class SingleSelectionHandler extends CheckBoxList.Handler{
//
//            public SingleSelectionHandler(CheckBoxList list) {
//                super(list);
//            }
//
//            protected void toggleSelection(int index) {
//                if(_list.getCheckBoxListSelectionModel().isSelectedIndex(index)){
//                    return;
//                }
//                super.toggleSelection(index);
//            }
//        }
//
//    }
//
//    public static class RadioListCellRenderer extends CheckBoxListCellRenderer {
//
//        protected AbstractButton button = new NullRadioButton();
//
//        public RadioListCellRenderer() {
//            this(null);
//        }
//
//        public RadioListCellRenderer(ListCellRenderer renderer) {
//            super(renderer);
//            //Really have no idea why this is necessary, or how the check box gets added in the first place
//            if (getComponentCount() > 0)
//                remove(0);
//            button.setBorder(BorderFactory.createEmptyBorder(0, 2, 0, 2));
//            button.setOpaque(false);
//            add(button, BorderLayout.BEFORE_LINE_BEGINS);
//            set_checkBox(button);
//        }
//
//        public void set_checkBox(AbstractButton button){
//            //This is here for JFormDesigner
//            //It gets a NoSuchFieldError when doing this.
//            //Doesn't really matter to the compiled program, only makes
//            //editing the dialog easier
//            try{
//                this._checkBox = button;
//            }catch(NoSuchFieldError e){
//
//            }
//
//        }
//    }


}
