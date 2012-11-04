/*
 * Created by JFormDesigner on Tue Jan 03 22:58:44 EST 2012
 */

package org.broad.igv.hic;

import org.broad.igv.renderer.ColorScale;
import slider.RangeSlider;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;

/**
 * @author Jim Robinson
 */
public class ColorRangeDialog extends JDialog {

    RangeSlider colorSlider;

    public ColorRangeDialog(Frame owner, RangeSlider colorSlider) {
        super(owner);
        initComponents();
        this.colorSlider = colorSlider;
        minimumField.setText(String.valueOf(colorSlider.getMinimum()));
        maximumField.setText(String.valueOf(colorSlider.getMaximum()));
        tickSpacingField.setText(String.valueOf(colorSlider.getMajorTickSpacing()));
    }


    private void okButtonActionPerformed(ActionEvent e) {
        int max = Integer.parseInt(maximumField.getText());
        int min = Integer.parseInt(minimumField.getText());
        int tickSpacing = Integer.parseInt(tickSpacingField.getText());
        colorSlider.setMinimum(min);
        colorSlider.setMaximum(max);
        colorSlider.setMajorTickSpacing(tickSpacing);
        colorSlider.setMinorTickSpacing(tickSpacing);
        setVisible(false);
    }

    private void cancelButtonActionPerformed(ActionEvent e) {
        setVisible(false);
    }

    private void initComponents() {
        // JFormDesigner - Component initialization - DO NOT MODIFY  //GEN-BEGIN:initComponents
        // Generated using JFormDesigner non-commercial license
        dialogPane = new JPanel();
        panel3 = new JPanel();
        label1 = new JLabel();
        contentPanel = new JPanel();
        panel2 = new JPanel();
        label4 = new JLabel();
        minimumField = new JTextField();
        panel1 = new JPanel();
        label2 = new JLabel();
        maximumField = new JTextField();
        panel4 = new JPanel();
        label3 = new JLabel();
        tickSpacingField = new JTextField();
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

            //======== panel3 ========
            {
                panel3.setLayout(new FlowLayout(FlowLayout.LEFT, 5, 25));

                //---- label1 ----
                label1.setText("Set color slider control range");
                panel3.add(label1);
            }
            dialogPane.add(panel3, BorderLayout.NORTH);

            //======== contentPanel ========
            {
                contentPanel.setLayout(new BoxLayout(contentPanel, BoxLayout.Y_AXIS));

                //======== panel2 ========
                {
                    panel2.setLayout(new FlowLayout(FlowLayout.LEFT));

                    //---- label4 ----
                    label4.setText("Minimum:");
                    panel2.add(label4);

                    //---- minimumField ----
                    minimumField.setText("text");
                    minimumField.setMaximumSize(new Dimension(100, 16));
                    minimumField.setMinimumSize(new Dimension(100, 16));
                    minimumField.setPreferredSize(new Dimension(100, 16));
                    panel2.add(minimumField);
                }
                contentPanel.add(panel2);

                //======== panel1 ========
                {
                    panel1.setLayout(new FlowLayout(FlowLayout.LEADING));

                    //---- label2 ----
                    label2.setText("Maximum");
                    panel1.add(label2);

                    //---- maximumField ----
                    maximumField.setText("text");
                    maximumField.setMaximumSize(new Dimension(100, 16));
                    maximumField.setMinimumSize(new Dimension(100, 16));
                    maximumField.setPreferredSize(new Dimension(100, 16));
                    panel1.add(maximumField);
                }
                contentPanel.add(panel1);

                //======== panel4 ========
                {
                    panel4.setLayout(new FlowLayout(FlowLayout.LEADING));

                    //---- label3 ----
                    label3.setText("Tick spacing:");
                    panel4.add(label3);

                    //---- tickSpacingField ----
                    tickSpacingField.setText("text");
                    tickSpacingField.setMaximumSize(new Dimension(100, 16));
                    tickSpacingField.setMinimumSize(new Dimension(100, 16));
                    tickSpacingField.setPreferredSize(new Dimension(100, 16));
                    panel4.add(tickSpacingField);
                }
                contentPanel.add(panel4);
            }
            dialogPane.add(contentPanel, BorderLayout.CENTER);

            //======== buttonBar ========
            {
                buttonBar.setBorder(new EmptyBorder(12, 0, 0, 0));
                buttonBar.setLayout(new GridBagLayout());
                ((GridBagLayout)buttonBar.getLayout()).columnWidths = new int[] {0, 85, 80};
                ((GridBagLayout)buttonBar.getLayout()).columnWeights = new double[] {1.0, 0.0, 0.0};

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
        contentPane.add(dialogPane, BorderLayout.CENTER);
        pack();
        setLocationRelativeTo(getOwner());
        // JFormDesigner - End of component initialization  //GEN-END:initComponents
    }

    // JFormDesigner - Variables declaration - DO NOT MODIFY  //GEN-BEGIN:variables
    // Generated using JFormDesigner non-commercial license
    private JPanel dialogPane;
    private JPanel panel3;
    private JLabel label1;
    private JPanel contentPanel;
    private JPanel panel2;
    private JLabel label4;
    private JTextField minimumField;
    private JPanel panel1;
    private JLabel label2;
    private JTextField maximumField;
    private JPanel panel4;
    private JLabel label3;
    private JTextField tickSpacingField;
    private JPanel buttonBar;
    private JButton okButton;
    private JButton cancelButton;
    // JFormDesigner - End of variables declaration  //GEN-END:variables
}
