package org.broad.igv.ucsc.hub;

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

    public static final Color ALT_BACKGROUND = new Color(245, 245, 245);
    private final Map<String, JTextArea> cellRendererCache = new HashMap<>();

    public HubTable(Object[][] data, String[] columnNames) {

        super(data, columnNames);

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

        getColumnModel().getColumn(1).setCellRenderer(getDefaultTableCellRenderer());
        getColumnModel().getColumn(2).setCellRenderer(getDefaultTableCellRenderer());
        getColumnModel().getColumn(3).setCellRenderer(getDefaultTableCellRenderer());
        getColumnModel().getColumn(3).setMaxWidth(150);
        getColumnModel().getColumn(3).setPreferredWidth(150);

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
                if (column == 3) {
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
                if (column == 3) {
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
    }


    @Override
    public Class<?> getColumnClass(int column) {
        return column == 0 ? Boolean.class : (column == 4 ? JLabel.class : String.class);
    }

    @Override
    public void changeSelection(int rowIndex, int columnIndex, boolean toggle, boolean extend) {
        // Suppress both row and column selection
        super.changeSelection(-1, -1, toggle, extend);
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


    private DefaultTableCellRenderer getDefaultTableCellRenderer() {
        return new DefaultTableCellRenderer() {
            @Override
            public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                String cacheKey = row + ":" + column;
                JTextArea textArea = cellRendererCache.computeIfAbsent(cacheKey, key -> {
                    JTextArea newTextArea = new JTextArea(value != null ? value.toString() : "");
                    newTextArea.setLineWrap(true);
                    newTextArea.setWrapStyleWord(true);
                    newTextArea.setOpaque(true);
                    newTextArea.setFont(table.getFont()); // Ensure font matches the table's font
                    return newTextArea;
                });

                textArea.setBackground(isSelected ? table.getSelectionBackground() : table.getBackground());
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

    @Override
    public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
        Component component = super.prepareRenderer(renderer, row, column);
        component.setBackground(row % 2 == 0 ? ALT_BACKGROUND : Color.WHITE);
        return component;
    }


    @Override
    public boolean isCellEditable(int row, int column) {
        return column == 0;
    }

}
