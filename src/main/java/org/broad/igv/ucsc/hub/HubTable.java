package org.broad.igv.ucsc.hub;

import org.broad.igv.Globals;
import org.broad.igv.ui.util.IconFactory;

import javax.swing.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.net.URI;
import java.util.HashMap;
import java.util.Map;

class HubTable extends JTable {

    public static  Color ALT_BACKGROUND = new Color(245, 245, 245);
    private final Map<Integer, Integer> rowHeights = new HashMap<>();

    public HubTable(HubTableModel model) {

        super(model);

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
        getColumnModel().getColumn(0).setCellEditor(new DefaultCellEditor(new JCheckBox()) {
            @Override
            public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
                JCheckBox checkBox = (JCheckBox) super.getTableCellEditorComponent(table, value, isSelected, row, column);
                //checkBox.setSelected(value != null && (Boolean) value);
                checkBox.addActionListener(e -> {
                    setValueAt(checkBox.isSelected(), row, column);
                });
                return checkBox;
            }
        });

        getColumnModel().getColumn(1).setCellRenderer(getMultiLineCellRenderer(Integer.MAX_VALUE));
        getColumnModel().getColumn(2).setCellRenderer(getMultiLineCellRenderer(Integer.MAX_VALUE));

        getColumnModel().getColumn(3).setPreferredWidth(150);
        getColumnModel().getColumn(3).setMaxWidth(200);
        getColumnModel().getColumn(3).setCellRenderer(getMultiLineCellRenderer(100));

        getColumnModel().getColumn(4).setMaxWidth(30);
        getColumnModel().getColumn(4).setCellRenderer(new DefaultTableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                ImageIcon icon = IconFactory.getInstance().getIcon(IconFactory.IconID.INFO);
                JLabel label = new JLabel(icon);
                label.setToolTipText(value.toString());
                label.setOpaque(true);
                label.setBackground(table.getBackground());
                return label;
            }
        });

        // Add a MouseMotionListener to handle mouse movement over column 3
        addMouseMotionListener(new MouseMotionAdapter() {
            @Override
            public void mouseMoved(MouseEvent e) {
                int column = columnAtPoint(e.getPoint());
                if (column == 4) {
                    setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
                } else {
                    setCursor(Cursor.getDefaultCursor());
                }
            }
        });

        // Add a MouseListener to handle clicks on column 3
        addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                int row = rowAtPoint(e.getPoint());
                int column = columnAtPoint(e.getPoint());
                if (column == 4) {
                    Object value = getValueAt(row, column);
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

        if(Globals.isDarkMode()) {
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
    }

    private static DefaultTableCellRenderer getMultiLineCellRenderer(int maxLength) {
        return new DefaultTableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                String v = value != null ? value.toString() : "";
                if(v.length() > maxLength) {
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
                int preferredHeight = textArea.getPreferredSize().height;
                if (preferredHeight > table.getRowHeight(row)) {
                    table.setRowHeight(row, preferredHeight);
                }

                return textArea;
            }
        };
    }

    @Override
    public Class<?> getColumnClass(int column) {
        return column == 0 ? Boolean.class : (column == 4 ? JLabel.class : String.class);
    }

    @Override
    public boolean isCellSelected(int row, int column) {
        return false;
    }

    @Override
    public boolean isColumnSelected(int column) {
        return false;
    }

    @Override
    public boolean isRowSelected(int row) {
        return false;
    }

    @Override
    public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
        Component component = super.prepareRenderer(renderer, row, column);
        component.setBackground(row % 2 == 0 ? ALT_BACKGROUND : Globals.isDarkMode() ? Color.black : Color.WHITE);
        return component;
    }

    //
    @Override
    public boolean isCellEditable(int row, int column) {
        return column == 0;
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
