package org.igv.ucsc.hub;

import org.igv.Globals;
import org.igv.ui.util.IconFactory;

import javax.swing.*;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableRowSorter;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.net.URI;
import java.util.HashMap;
import java.util.Map;

class HubTable extends JTable {

    public static Color ALT_BACKGROUND = new Color(245, 245, 245);
    private final Map<Integer, Integer> rowHeights = new HashMap<>();

    /**
     * In single-select mode the checkbox column is hidden and a single row is chosen by row selection.  Rows whose
     * descriptor is flagged "selected" are shown disabled (greyed out) and cannot be selected.  In the default
     * (checkbox) mode the "selected" flag is the checkbox state.
     */
    private final boolean singleSelect;

    public HubTable(HubTableModel model) {
        this(model, false);
    }

    public HubTable(HubTableModel model, boolean singleSelect) {

        super(model);
        this.singleSelect = singleSelect;

        getTableHeader().setReorderingAllowed(false);

        setRowSorter(new TableRowSorter<>(model));

        getColumnModel().getColumn(0).setMaxWidth(30);
        getColumnModel().getColumn(0).setCellRenderer(new DefaultTableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                JCheckBox checkBox = new JCheckBox();
                checkBox.setSelected(value != null && (Boolean) value);
                checkBox.setHorizontalAlignment(SwingConstants.CENTER);
                checkBox.setVerticalAlignment(SwingConstants.TOP);
                checkBox.setOpaque(true);
                checkBox.setBackground(table.getBackground());
                return checkBox;
            }
        });

        // Set a cell editor for column 0 to make the checkbox interactive
        getColumnModel().getColumn(0).setCellEditor(new DefaultCellEditor(new JCheckBox()));

        getColumnModel().getColumn(1).setCellRenderer(getMultiLineCellRenderer(Integer.MAX_VALUE));
        getColumnModel().getColumn(2).setCellRenderer(getMultiLineCellRenderer(Integer.MAX_VALUE));

        getColumnModel().getColumn(3).setPreferredWidth(150);
        getColumnModel().getColumn(3).setMaxWidth(200);
        getColumnModel().getColumn(3).setCellRenderer(getMultiLineCellRenderer(100));

        getColumnModel().getColumn(4).setMaxWidth(30);
        getColumnModel().getColumn(4).setCellRenderer(new DefaultTableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {

                if (value != null) {
                    ImageIcon icon = IconFactory.getInstance().getIcon(IconFactory.IconID.INFO);
                    JLabel label = new JLabel(icon);
                    label.setToolTipText(value.toString());
                    label.setOpaque(true);
                    label.setBackground(table.getBackground());
                    return label;
                } else {
                    return new JLabel();
                }
            }
        });

        // Add a MouseMotionListener to handle mouse movement over the info column (model column 4)
        addMouseMotionListener(new MouseMotionAdapter() {
            @Override
            public void mouseMoved(MouseEvent e) {
                int column = columnAtPoint(e.getPoint());
                if (column >= 0 && convertColumnIndexToModel(column) == 4) {
                    setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
                } else {
                    setCursor(Cursor.getDefaultCursor());
                }
            }
        });

        // Add a MouseListener to handle clicks on the info column (model column 4)
        addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                int row = rowAtPoint(e.getPoint());
                int column = columnAtPoint(e.getPoint());
                if (row >= 0 && column >= 0 && convertColumnIndexToModel(column) == 4) {
                    Object value = getModel().getValueAt(convertRowIndexToModel(row), 4);
                    if (value != null) {
                        try {
                            Desktop.getDesktop().browse(new URI(value.toString()));
                        } catch (Exception ex) {
                            // Handle exception
                        }
                    }
                }
            }
        });

        if (Globals.isDarkMode()) {
            setBackground(Color.DARK_GRAY);
            setForeground(Color.WHITE);
            setSelectionBackground(Color.GRAY);
            setSelectionForeground(Color.WHITE);
            ALT_BACKGROUND = new Color(60, 60, 60); // Darker alternate background
            // Update renderers
            for (int i = 0; i < getColumnModel().getColumnCount(); i++) {
                TableCellRenderer renderer = getColumnModel().getColumn(i).getCellRenderer();
                if (renderer instanceof DefaultTableCellRenderer) {
                    ((DefaultTableCellRenderer) renderer).setBackground(getBackground());
                    ((DefaultTableCellRenderer) renderer).setForeground(getForeground());
                }
            }
        }

        if (singleSelect) {
            // Hide the checkbox column.  Model column indices are unchanged.
            removeColumn(getColumnModel().getColumn(0));
            setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
            setRowSelectionAllowed(true);
            setColumnSelectionAllowed(false);
        }
    }

    private static DefaultTableCellRenderer getMultiLineCellRenderer(int maxLength) {
        return new DefaultTableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                String v = value != null ? value.toString() : "";
                if (v.length() > maxLength) {
                    v = v.substring(0, maxLength) + "...";
                }
                JTextArea textArea = new JTextArea(v);
                textArea.setToolTipText(value.toString());
                textArea.setLineWrap(true);
                textArea.setWrapStyleWord(true);
                textArea.setOpaque(true);
                textArea.setFont(table.getFont()); // Ensure font matches the table's font
                textArea.setBackground(table.getBackground());
                textArea.setForeground(table.getForeground());

                // Set the width of the JTextArea to match the column width
                int columnWidth = table.getColumnModel().getColumn(column).getWidth();
                textArea.setSize(columnWidth, Short.MAX_VALUE);

                // Adjust row height to fit the content if neccessary
                int preferredHeight = textArea.getPreferredSize().height + 5;
                if (preferredHeight > table.getRowHeight(row)) {
                    table.setRowHeight(row, preferredHeight);
                }

                return textArea;
            }
        };
    }

    /**
     * A row is disabled, in single-select mode, when its descriptor is flagged selected (already added).
     */
    private boolean isDisabledRow(int viewRow) {
        return singleSelect && ((HubTableModel) getModel()).getRecords().get(convertRowIndexToModel(viewRow)).isSelected();
    }

    @Override
    public void changeSelection(int rowIndex, int columnIndex, boolean toggle, boolean extend) {
        if (isDisabledRow(rowIndex)) {
            return;
        }
        super.changeSelection(rowIndex, columnIndex, toggle, extend);
    }

    @Override
    public Class<?> getColumnClass(int column) {
        int modelColumn = convertColumnIndexToModel(column);
        return modelColumn == 0 ? Boolean.class : (modelColumn == 4 ? JLabel.class : String.class);
    }

    @Override
    public boolean isCellSelected(int row, int column) {
        return singleSelect ? super.isCellSelected(row, column) : false;
    }

    @Override
    public boolean isColumnSelected(int column) {
        return singleSelect ? super.isColumnSelected(column) : false;
    }

    @Override
    public boolean isRowSelected(int row) {
        return singleSelect ? super.isRowSelected(row) : false;
    }

    @Override
    public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
        Component component = super.prepareRenderer(renderer, row, column);
        if (isRowSelected(row)) {
            component.setBackground(getSelectionBackground());
            component.setForeground(getSelectionForeground());
        } else {
            component.setBackground(row % 2 == 0 ? ALT_BACKGROUND : Globals.isDarkMode() ? Color.black : Color.WHITE);
            if (isDisabledRow(row)) {
                component.setForeground(UIManager.getColor("Label.disabledForeground"));
            }
        }
        return component;
    }

    //
    @Override
    public boolean isCellEditable(int row, int column) {
        return !singleSelect && column == 0;
    }


    private DefaultTableCellRenderer getDefaultTableCellRenderer() {
        return new DefaultTableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                JTextArea textArea = new JTextArea(value != null ? value.toString() : "");
                textArea.setLineWrap(true);
                textArea.setWrapStyleWord(true);
                textArea.setOpaque(true);
                textArea.setFont(table.getFont()); // Ensure font matches the table's font
                textArea.setBackground(table.getBackground());
                textArea.setForeground(table.getForeground());

                // Set the width of the JTextArea to match the column width
                int columnWidth = table.getColumnModel().getColumn(column).getWidth();
                textArea.setSize(columnWidth, Short.MAX_VALUE);

                // Adjust row height to fit the content
                int preferredHeight = textArea.getPreferredSize().height;
                if (table.getRowHeight(row) != preferredHeight) {
                    table.setRowHeight(row, preferredHeight);
                }

                return textArea;
            }
        };
    }

}
