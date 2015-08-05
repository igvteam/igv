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
 * Created by JFormDesigner on Wed May 16 16:09:07 EDT 2012
 */

package org.broad.igv.track;

import com.google.common.base.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.cli_plugin.ui.TrackArgument;
import org.broad.igv.data.CombinedDataSource;
import org.broad.igv.ui.IGV;
import org.broad.igv.ui.PanelName;
import org.broad.igv.util.StringUtils;
import org.broad.igv.util.collections.CollUtils;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * @author jacob
 */
public class CombinedDataSourceDialog extends JDialog {

    private static Logger log = Logger.getLogger(CombinedDataSourceDialog.class);


    public CombinedDataSourceDialog(Frame owner) {
        super(owner);
        initComponents();

        ArrayList<CombinedDataSource.Operation> dialogOperations = new ArrayList<CombinedDataSource.Operation>(
                Arrays.asList(CombinedDataSource.Operation.values()));

        operation.setModel(new DefaultComboBoxModel(dialogOperations.toArray()));

        List<DataTrack> visibleTracks = CollUtils.filter(IGV.getInstance().getDataTracks(), new Predicate<DataTrack>() {
            @Override
            public boolean apply(DataTrack input) {
                return input.isVisible();
            }
        });

        trackABox.setModel(new DefaultComboBoxModel(visibleTracks.toArray()));
        trackBBox.setModel(new DefaultComboBoxModel(visibleTracks.toArray()));
        //Show 1st and 2nd track by default
        if(visibleTracks.size() >= 2) trackBBox.setSelectedIndex(1);

        operation.setRenderer(new OperationComboBoxRenderer());

        trackABox.setRenderer(new TrackArgument.TrackComboBoxRenderer());
        trackBBox.setRenderer(new TrackArgument.TrackComboBoxRenderer());

        setOutputTrackName();
        operation.addItemListener(new SetOutputTrackNameListener());
    }

    public CombinedDataSourceDialog(Frame owner, Iterator<Track> tracks) {
        this(owner);

        trackABox.setSelectedItem(tracks.next());
        trackBBox.setSelectedItem(tracks.next());
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        this.setVisible(false);
    }

    private void okButtonActionPerformed(ActionEvent e) {

        DataTrack track0 = (DataTrack) trackABox.getSelectedItem();
        DataTrack track1 = (DataTrack) trackBBox.getSelectedItem();
        CombinedDataSource.Operation op = (CombinedDataSource.Operation) operation.getSelectedItem();
        String text = resultName.getText();

        CombinedDataSource source = new CombinedDataSource(track0, track1, op);

        DataSourceTrack newTrack = new DataSourceTrack(null, track0.getId() + track1.getId() + text, text, source);
        TrackMenuUtils.changeRenderer(Arrays.<Track>asList(newTrack), track0.getRenderer().getClass());
        newTrack.setDataRange(track0.getDataRange());
        newTrack.setColorScale(track0.getColorScale());
        IGV.getInstance().addTracks(Arrays.<Track>asList(newTrack), PanelName.DATA_PANEL);

        this.setVisible(false);

        IGV.getInstance().repaint();
    }

    private void helpButtonActionPerformed(ActionEvent e) {
        String defInfo = "Error retrieving help. See the IGV User Guide.";
        InputStream is = this.getClass().getResourceAsStream("/resources/bedtools_help.txt");
        String info = "";
        String line;
        BufferedReader reader = new BufferedReader(new InputStreamReader(is));
        try {
            while ((line = reader.readLine()) != null) {
                info += line + "\n";
            }
        } catch (IOException exc) {
            log.error(exc);
            info = defInfo;
        }

        JTextArea textArea = new JTextArea(info);
        textArea.setLineWrap(true);
        textArea.setWrapStyleWord(true);

        JScrollPane scrollPane = new JScrollPane(textArea);

        JOptionPane pane = new JOptionPane(scrollPane, JOptionPane.PLAIN_MESSAGE);
        JDialog dialog = pane.createDialog(null, "Help");
        dialog.setAlwaysOnTop(true);
        dialog.setResizable(true);

        dialog.setSize(400, 400);
        dialog.setVisible(true);
        //JOptionPane.showMessageDialog(this, scrollPane, "Help", JOptionPane.PLAIN_MESSAGE);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        contentPanel = new JPanel();
        label2 = new JLabel();
        trackABox = new JComboBox();
        label4 = new JLabel();
        operation = new JComboBox();
        label3 = new JLabel();
        trackBBox = new JComboBox();
        label1 = new JLabel();
        scrollPane1 = new JScrollPane();
        resultName = new JTextArea();
        buttonBar = new JPanel();
        okButton = new JButton();
        helpButton = new JButton();
        cancelButton = new JButton();

        //======== this ========
        setAlwaysOnTop(true);
        Container contentPane = getContentPane();
        contentPane.setLayout(new BorderLayout());

        //======== dialogPane ========
        {
            dialogPane.setBorder(new EmptyBorder(12, 12, 12, 12));
            dialogPane.setLayout(new BorderLayout());

            //======== contentPanel ========
            {
                contentPanel.setLayout(null);

                //---- label2 ----
                label2.setText("Track A:");
                label2.setLabelFor(trackABox);
                contentPanel.add(label2);
                label2.setBounds(new Rectangle(new Point(5, 20), label2.getPreferredSize()));
                contentPanel.add(trackABox);
                trackABox.setBounds(70, 15, 205, trackABox.getPreferredSize().height);

                //---- label4 ----
                label4.setText("Operation:");
                contentPanel.add(label4);
                label4.setBounds(new Rectangle(new Point(5, 65), label4.getPreferredSize()));
                contentPanel.add(operation);
                operation.setBounds(70, 60, 150, operation.getPreferredSize().height);

                //---- label3 ----
                label3.setText("Track B:");
                label3.setLabelFor(trackBBox);
                contentPanel.add(label3);
                label3.setBounds(5, 110, 52, 16);
                contentPanel.add(trackBBox);
                trackBBox.setBounds(70, 105, 205, trackBBox.getPreferredSize().height);

                //---- label1 ----
                label1.setText("Result Track Name");
                contentPanel.add(label1);
                label1.setBounds(50, 150, 128, label1.getPreferredSize().height);

                //======== scrollPane1 ========
                {

                    //---- resultName ----
                    resultName.setLineWrap(true);
                    scrollPane1.setViewportView(resultName);
                }
                contentPanel.add(scrollPane1);
                scrollPane1.setBounds(10, 170, 270, 65);

                { // compute preferred size
                    Dimension preferredSize = new Dimension();
                    for(int i = 0; i < contentPanel.getComponentCount(); i++) {
                        Rectangle bounds = contentPanel.getComponent(i).getBounds();
                        preferredSize.width = Math.max(bounds.x + bounds.width, preferredSize.width);
                        preferredSize.height = Math.max(bounds.y + bounds.height, preferredSize.height);
                    }
                    Insets insets = contentPanel.getInsets();
                    preferredSize.width += insets.right;
                    preferredSize.height += insets.bottom;
                    contentPanel.setMinimumSize(preferredSize);
                    contentPanel.setPreferredSize(preferredSize);
                }
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 80, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

                //---- okButton ----
                okButton.setText("  OK  ");
                okButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        okButtonActionPerformed(e);
                    }
                });
                buttonBar.add(okButton, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.VERTICAL,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- helpButton ----
                helpButton.setText("Help");
                helpButton.setVisible(false);
                helpButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        helpButtonActionPerformed(e);
                    }
                });
                buttonBar.add(helpButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.VERTICAL,
                    new Insets(0, 0, 0, 0), 0, 0));

                //---- cancelButton ----
                cancelButton.setText("Cancel");
                cancelButton.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        cancelButtonActionPerformed(e);
                    }
                });
                buttonBar.add(cancelButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0,
                    GridBagConstraints.CENTER, GridBagConstraints.VERTICAL,
                    new Insets(0, 0, 0, 0), 0, 0));
            }
            dialogPane.add(buttonBar, BorderLayout.SOUTH);
        }
        contentPane.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel contentPanel;
    private JLabel label2;
    private JComboBox trackABox;
    private JLabel label4;
    private JComboBox operation;
    private JLabel label3;
    private JComboBox trackBBox;
    private JLabel label1;
    private JScrollPane scrollPane1;
    private JTextArea resultName;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton helpButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables

    private class OperationComboBoxRenderer extends DefaultListCellRenderer {
        @Override
        public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
            CombinedDataSource.Operation op = (CombinedDataSource.Operation) value;
            String toShow = StringUtils.capWords(op.name());
            return super.getListCellRendererComponent(list, toShow, index, isSelected, cellHasFocus);
        }

    }

    private class SetOutputTrackNameListener implements ItemListener {

        @Override
        public void itemStateChanged(ItemEvent e) {
            setOutputTrackName();
        }
    }

    private void setOutputTrackName() {
        String name = "Sum";
        CombinedDataSource.Operation op = (CombinedDataSource.Operation) operation.getSelectedItem();
        switch (op){
            case ADD:
                break;
            case SUBTRACT:
                name = "Difference";
                break;
        }
        resultName.setText(name);

//        String name = ((Track) trackABox.getSelectedItem()).getName();
//        name += " " + operation.getSelectedItem() + " ";
//        name += ((Track) trackBBox.getSelectedItem()).getName();

    }

}
