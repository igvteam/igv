package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.Genome;

import javax.swing.*;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

public class HubSelectionDialog extends JDialog {

    private static final long serialVersionUID = 1L;

    private static Map<String, HubSelectionDialog> instances = new java.util.HashMap<>();

    String ucscId;
    List<HubDescriptor> hubs;
    private boolean canceled;
    private JTable table;

    public HubSelectionDialog(Frame owner, Genome genome) {

        super(owner, "Select Track Hubs", true);
        this.ucscId = genome.getUCSCId();

        setSize(owner.getWidth() - 20, 600);

        if(owner != null) {
            setLocationRelativeTo(owner);
           // setLocation(getX() - getWidth() / 2, getY() - getHeight() / 2);
        }


        final String ucscId = genome.getUCSCId();
        List<HubDescriptor> hubUrls = HubRegistry.getPublicHubs (ucscId);
        if (hubUrls != null && !hubUrls.isEmpty()) {
            this.hubs = hubUrls;
        } else {
            new ArrayList<>();
        }

        Collection<Hub> loadedHubs = genome.getTrackHubs();

        initComponents(loadedHubs);
    }


    public List<Hub> getSelectedHubs() {
        List<Hub> selectedHubs = new ArrayList<>();
        for (int i = 0; i < table.getRowCount(); i++) {
            Boolean isSelected = (Boolean) table.getValueAt(i, 0);
            if (isSelected) {
                try {
                    selectedHubs.add(HubParser.loadHub(hubs.get(i).getUrl(), ucscId ));
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
        }
        return selectedHubs;
    }

    public Set<String> getUnselectedURLs() {
        Set<String> unselectedURLs = new HashSet<>();
        for (int i = 0; i < table.getRowCount(); i++) {
            Boolean isSelected = (Boolean) table.getValueAt(i, 0);
            if (!isSelected) {
                unselectedURLs.add(hubs.get(i).getUrl());
            }
        }
        return unselectedURLs;
    }


    private void initComponents(Collection<Hub> loadedHubs) {

        this.setLayout(new BorderLayout());

        Set<String> loadedURLs = loadedHubs.stream().map(Hub::getUrl).collect(Collectors.toSet());

        String[] columnNames = {"", "Name", "Description"};
        Object[][] data = new Object[this.hubs.size()][3];
        for (int i = 0; i < this.hubs.size(); i++) {
            HubDescriptor hub = this.hubs.get(i);
            data[i][0] = loadedURLs.contains(hub.getUrl()); // Default unchecked checkbox
            data[i][1] = hub.getShortLabel();
            data[i][2] = hub.getLongLabel();
        }

        table = new HubTable(data, columnNames);

        JScrollPane scrollPane = new JScrollPane(table);
        add(scrollPane, BorderLayout.CENTER);

        JPanel buttonPanel = new JPanel();
        ((FlowLayout) buttonPanel.getLayout()).setAlignment(FlowLayout.RIGHT);

        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(e -> {
            canceled = true;
            setVisible(false);
            dispose();
        });

        JButton okButton = new JButton("OK");
        okButton.addActionListener(e -> {
            canceled = false;
            setVisible(false);
            dispose();
        });

//        if (Globals.IS_MAC) {
//            buttonPanel.add(cancelButton);
//            buttonPanel.add(okButton);
//        } else {
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);
        //}
        add(buttonPanel, BorderLayout.SOUTH);
    }
}

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
                checkBox.setBackground(isSelected ? table.getSelectionBackground() : table.getBackground());
                return checkBox;
            }
        });
    }

    @Override
    public Class<?> getColumnClass(int column) {
        return column == 0 ? Boolean.class : String.class;
    }

    @Override
    public TableCellRenderer getCellRenderer(int row, int column) {
        if (column == 1 || column == 2) { // Columns for "Name" and "Description"
            return new DefaultTableCellRenderer() {
                @Override
                public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
                    String cacheKey = row + ":" + column;
                    JTextArea textArea = cellRendererCache.computeIfAbsent(cacheKey, key -> {
                        JTextArea newTextArea = new JTextArea(value != null ? value.toString() : "");
                        newTextArea.setLineWrap(true);
                        newTextArea.setWrapStyleWord(true);
                        newTextArea.setOpaque(true);
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
        return super.getCellRenderer(row, column);
    }

    @Override
    public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
        Component component = super.prepareRenderer(renderer, row, column);
        if (!isRowSelected(row)) {
            component.setBackground(row % 2 == 0 ? ALT_BACKGROUND : Color.WHITE);
        }
        return component;
    }

    @Override
    public boolean isCellSelected(int row, int column) {
        return false;
    }

    @Override
    public boolean isCellEditable(int row, int column) {
        return column == 0;
    }



}