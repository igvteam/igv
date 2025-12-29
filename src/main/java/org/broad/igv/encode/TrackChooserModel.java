package org.broad.igv.encode;

import org.broad.igv.ui.IGV;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import javax.swing.table.TableStringConverter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * @author jrobinso
 * Date: 10/31/13
 * Time: 10:09 PM
 */
public class TrackChooserModel extends AbstractTableModel {

    private String[] columnHeadings;
    private List<FileRecord> records;
    private final TableRowSorter<TrackChooserModel> sorter;

    public TrackChooserModel(List<String> headings, List<FileRecord> records) {

        this.records = records;

        List<String> tmp = new ArrayList<String>();
        tmp.add("");  // Checkbox heading
        for (String h : headings) {
            String heading = h.trim();
            if (heading.length() > 0 && !"path".equals(heading)) {
                tmp.add(heading);
            }
        }

        columnHeadings = tmp.toArray(new String[tmp.size()]);

        sorter = new TableRowSorter<>(this);

        sorter.setStringConverter(new TableStringConverter() {
            @Override
            public String toString(TableModel model, int row, int column) {
                final Object value = model.getValueAt(row, column);
                return value == null ? "" : value.toString();
            }
        });
    }

    public TableRowSorter<TrackChooserModel> getSorter() {
        return sorter;
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
        return columnIndex == 0 ? Boolean.class : String.class;
    }

    @Override
    public String getColumnName(int column) {
        return columnHeadings[column];
    }

    @Override
    public int getRowCount() {
        return records.size();
    }

    @Override
    public int getColumnCount() {
        return columnHeadings.length;
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {

        if (rowIndex >= records.size() || columnIndex >= columnHeadings.length) {
            return null;
        }

        FileRecord record = records.get(rowIndex);
        if (columnIndex == 0) {
            return record.isSelected();
        } else {
            String att = columnHeadings[columnIndex];
            return record.getAttributeValue(att);
        }

    }

    @Override
    public boolean isCellEditable(int rowIndex, int columnIndex) {
        return columnIndex == 0;
    }

    @Override
    public void setValueAt(Object value, int row, int col) {
        if (col == 0) {
            records.get(row).setSelected((Boolean) value);
        }
        fireTableCellUpdated(row, col);
    }

    public void updateSelections() {

        Set<String> loadedPaths = IGV.hasInstance() ? IGV.getInstance().getLoadedPaths() : Collections.emptySet();

        for (int row = 0; row < records.size(); row++) {
            FileRecord record = records.get(row);
            if (loadedPaths.contains(record.getPath())) {
                record.setSelected(true);
            }
        }
    }

    public List<FileRecord> getRecords() {
        return records;
    }
}