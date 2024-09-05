package org.broad.igv.variant;

import com.intellij.uiDesigner.core.GridConstraints;
import com.intellij.uiDesigner.core.GridLayoutManager;
import com.jidesoft.swing.AutoCompletion;
import com.jidesoft.swing.ListSearchable;
import htsjdk.variant.vcf.VCFCompoundHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broad.igv.renderer.ColorScale;
import org.broad.igv.renderer.ContinuousColorScale;
import org.broad.igv.ui.IGVDialog;
import org.broad.igv.ui.color.PaletteColorTable;
import org.broad.igv.ui.legend.ContinuousLegendPanel;
import org.broad.igv.ui.legend.DiscreteLegendPanel;
import org.broad.igv.ui.legend.LegendPanel;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.List;
import java.util.Optional;

public class SelectVcfFieldDialog extends IGVDialog {
    private JPanel contentPane;
    private JButton buttonOK;
    private JButton buttonCancel;
    private JTextField inputField;
    private JList<VCFCompoundHeaderLine> possibleValueList;
    private JTextArea description;
    private JPanel buttonPanel;
    private JPanel inputPanel;
    private JPanel displayPanel;
    private JTextField idLabel;
    private JTextField typeLabel;
    private JTextField countLabel;
    private JPanel colorLegendPanel;
    private JButton saveScaleButton;
    private LegendPanel numericalPanel;

    private String value;
    private ColorScale colorScale;
    private final AttributeColorManager.Type fieldType;

    public SelectVcfFieldDialog(Frame parent, String title, String defaultValue, List<? extends VCFCompoundHeaderLine> possibleValues) {
        super(parent, title, true);
        this.fieldType = switch(possibleValues.getFirst()){
            case VCFInfoHeaderLine info -> AttributeColorManager.Type.INFO;
            case VCFFormatHeaderLine format -> AttributeColorManager.Type.FORMAT;
            default -> AttributeColorManager.Type.INFO;
        };
        $$$setupUI$$$();
        setContentPane(contentPane);
        getRootPane().setDefaultButton(buttonOK);

        this.setTitle(title);

        //set cell renderer to change the value that gets displayed, there is probably a more elegant way to do this
        possibleValueList.setCellRenderer(new DefaultListCellRenderer() {
            @Override
            public Component getListCellRendererComponent(final JList<?> list, final Object value, final int index, final boolean isSelected, final boolean cellHasFocus) {
                String id = ((VCFCompoundHeaderLine) value).getID();
                return super.getListCellRendererComponent(list, id, index, isSelected, cellHasFocus);
            }
        });

        possibleValueList.setListData(possibleValues.toArray(new VCFCompoundHeaderLine[0]));
        final AutoCompletion autoCompletion = new AutoCompletion(inputField, new ListSearchable(possibleValueList) {
            @Override
            protected String convertElementToString(Object element) {
                return ((VCFCompoundHeaderLine) element).getID();
            }
        });
        autoCompletion.setStrict(false);

        buttonOK.addActionListener(e -> onOK());
        buttonCancel.addActionListener(e -> onCancel());

        // call onCancel() when cross is clicked
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        addWindowListener(new WindowAdapter() {
            @Override
            public void windowClosing(WindowEvent e) {
                onCancel();
            }
        });

        // call onCancel() on ESCAPE
        contentPane.registerKeyboardAction(e -> onCancel(), KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);

        possibleValueList.addListSelectionListener(e -> onSelectionChange());
        //set this after the listeners are set up so that the list is correctly populated
        inputField.setText(defaultValue);

        saveScaleButton.setEnabled(false);
        saveScaleButton.addActionListener(action -> onPersist());

        colorLegendPanel.revalidate();
        colorLegendPanel.repaint();
    }

    private static AttributeColorManager.ColorScaleType headerLineToDisplayType(VCFHeaderLineType type) {
        return switch (type) {
            case Integer, Float -> AttributeColorManager.ColorScaleType.CONTINUOUS;
            case String, Character -> AttributeColorManager.ColorScaleType.DISCRETE;
            case Flag -> AttributeColorManager.ColorScaleType.FLAG;
        };
    }

    private void onSelectionChange() {
        final VCFCompoundHeaderLine selectedValue = possibleValueList.getSelectedValue();
        if (selectedValue == null) {
            setToNoSelection();
        } else {
            setSelectedValue(selectedValue);
        }
        colorLegendPanel.add(saveScaleButton, BorderLayout.EAST);
        updateValues();
        colorLegendPanel.revalidate();
        colorLegendPanel.repaint();
    }

    private void setToNoSelection() {
        idLabel.setText("");
        countLabel.setText("");
        typeLabel.setText("");
        description.setText("");
        colorLegendPanel.removeAll();
        saveScaleButton.setEnabled(false);
    }

    private void setSelectedValue(VCFCompoundHeaderLine selectedValue) {
        String id = selectedValue.getID();
        idLabel.setText(id);
        countLabel.setText(getCountAsString(selectedValue));
        typeLabel.setText(selectedValue.getType().toString());
        description.setText(selectedValue.getDescription());
        colorLegendPanel.removeAll();
        numericalPanel = getLegendPanel(id, selectedValue.getType());
        colorLegendPanel.add(numericalPanel, BorderLayout.CENTER);
        saveScaleButton.setEnabled(true);
    }

    private LegendPanel getLegendPanel(String id, VCFHeaderLineType type) {
        ColorScale scale = AttributeColorManager.getColorTable(fieldType, id, headerLineToDisplayType(type));
        LegendPanel panel = switch (scale) {
            case ContinuousColorScale continuous -> new ContinuousLegendPanel(id, continuous);
            case PaletteColorTable discrete -> new DiscreteLegendPanel(discrete);
            default -> throw new IllegalArgumentException("Unexpected ColorScale type " + scale.asString());
        };

        return panel;
    }


    /**
     * @param headerLine
     * @return
     */
    private static String getCountAsString(VCFCompoundHeaderLine headerLine) {
        return headerLine.isFixedCount() ? Integer.toString(headerLine.getCount()) : headerLine.getCountType().toString();
    }

    private void onPersist() {
        if (value != null && !value.isBlank() && colorScale != null) {
            //   PreferencesManager.updateAll( Map.of(AttributeColorManager.getPreferencesKey(AttributeColorManager.Type.INFO, value), colorScale.asString()););
        }
    }

    private void onOK() {
        updateValues();
        dispose();
    }

    private void updateValues() {
        value = inputField.getText();
        colorScale = numericalPanel.getColorScale();
    }

    private void onCancel() {
        value = null;
        colorScale = null;
        dispose();
    }

    public static Optional<ColorResult> showValueChooseDialog(Frame parent, String title, String defaultValue, List<? extends VCFCompoundHeaderLine> choices) {
        SelectVcfFieldDialog dialog = new SelectVcfFieldDialog(parent, title, defaultValue, choices);
        dialog.pack();
        dialog.setLocationRelativeTo(parent);
        dialog.setVisible(true);
        return dialog.getResult();
    }

    private Optional<ColorResult> getResult() {
        return value == null || colorScale == null
                ? Optional.empty()
                : Optional.of(new ColorResult(value, colorScale));
    }

    public static void main(String[] args) {
        SelectVcfFieldDialog dialog = new SelectVcfFieldDialog(new JFrame(), "Select Info Field to Display", "", List.of());

        dialog.pack();
        dialog.setVisible(true);
        System.exit(0);
    }

    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * >>> IMPORTANT!! <<<
     * DO NOT edit this method OR call it in your code!
     *
     * @noinspection ALL
     */
    private void $$$setupUI$$$() {
        contentPane = new JPanel();
        contentPane.setLayout(new GridLayoutManager(2, 1, new Insets(10, 10, 10, 10), -1, -1));
        final JSplitPane splitPane1 = new JSplitPane();
        splitPane1.setOrientation(1);
        contentPane.add(splitPane1, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        final JPanel panel1 = new JPanel();
        panel1.setLayout(new GridLayoutManager(1, 1, new Insets(0, 0, 0, 0), -1, -1));
        splitPane1.setLeftComponent(panel1);
        inputPanel = new JPanel();
        inputPanel.setLayout(new GridLayoutManager(2, 4, new Insets(0, 0, 0, 0), -1, -1));
        panel1.add(inputPanel, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        final JScrollPane scrollPane1 = new JScrollPane();
        inputPanel.add(scrollPane1, new GridConstraints(1, 0, 1, 4, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null, 0, false));
        possibleValueList = new JList();
        final DefaultListModel defaultListModel1 = new DefaultListModel();
        possibleValueList.setModel(defaultListModel1);
        possibleValueList.setSelectionMode(0);
        scrollPane1.setViewportView(possibleValueList);
        inputField = new JTextField();
        inputPanel.add(inputField, new GridConstraints(0, 0, 1, 4, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, new Dimension(150, -1), null, 0, false));
        final JSplitPane splitPane2 = new JSplitPane();
        splitPane2.setContinuousLayout(false);
        splitPane2.setDividerSize(1);
        splitPane2.setEnabled(true);
        splitPane2.setOrientation(0);
        splitPane1.setRightComponent(splitPane2);
        displayPanel = new JPanel();
        displayPanel.setLayout(new GridLayoutManager(2, 3, new Insets(0, 0, 0, 0), -1, -1));
        splitPane2.setLeftComponent(displayPanel);
        final JPanel panel2 = new JPanel();
        panel2.setLayout(new GridLayoutManager(3, 3, new Insets(0, 0, 0, 0), -1, -1));
        displayPanel.add(panel2, new GridConstraints(0, 0, 1, 3, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        final JLabel label1 = new JLabel();
        label1.setText("ID");
        panel2.add(label1, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        final JLabel label2 = new JLabel();
        label2.setText("Type");
        panel2.add(label2, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        final JLabel label3 = new JLabel();
        label3.setText("Count");
        panel2.add(label3, new GridConstraints(0, 2, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_NONE, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        idLabel = new JTextField();
        idLabel.setEditable(false);
        idLabel.setText(" ");
        panel2.add(idLabel, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        typeLabel = new JTextField();
        typeLabel.setEditable(false);
        typeLabel.setEnabled(true);
        typeLabel.setText(" ");
        panel2.add(typeLabel, new GridConstraints(1, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        countLabel = new JTextField();
        countLabel.setEditable(false);
        countLabel.setText(" ");
        panel2.add(countLabel, new GridConstraints(1, 2, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        description = new JTextArea();
        description.setEditable(false);
        description.setLineWrap(true);
        description.setText("");
        description.setWrapStyleWord(true);
        displayPanel.add(description, new GridConstraints(1, 0, 1, 3, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, new Dimension(500, 50), null, 0, false));
        colorLegendPanel = new JPanel();
        colorLegendPanel.setLayout(new BorderLayout(0, 0));
        colorLegendPanel.setToolTipText("heatmap");
        splitPane2.setRightComponent(colorLegendPanel);
        saveScaleButton = new JButton();
        saveScaleButton.setEnabled(false);
        saveScaleButton.setText("Persist");
        saveScaleButton.setToolTipText("Save the mapping between field and color for future sessions.");
        colorLegendPanel.add(saveScaleButton, BorderLayout.EAST);
        buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridLayoutManager(1, 2, new Insets(0, 0, 0, 0), -1, -1, true, false));
        contentPane.add(buttonPanel, new GridConstraints(1, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null, 0, false));
        buttonOK = new JButton();
        buttonOK.setText("OK");
        buttonPanel.add(buttonOK, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
        buttonCancel = new JButton();
        buttonCancel.setText("Cancel");
        buttonPanel.add(buttonCancel, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null, 0, false));
    }

    /**
     * @noinspection ALL
     */
    public JComponent $$$getRootComponent$$$() {
        return contentPane;
    }

    public record ColorResult(String value, ColorScale colors) {
    }
}
