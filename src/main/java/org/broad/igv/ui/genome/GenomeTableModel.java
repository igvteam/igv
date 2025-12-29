package org.broad.igv.ui.genome;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import javax.swing.table.TableStringConverter;
import java.util.List;

/**
 * @author jrobinso
 */
public class GenomeTableModel extends AbstractTableModel {

    private String[] columnHeadings;
    private List<GenomeListItem> records;
    private final TableRowSorter<GenomeTableModel> sorter;

    public GenomeTableModel(String [] headers, List<GenomeListItem> records) {

        this.records = records;

        //tmp.add("path");
        columnHeadings = headers;

        sorter = new TableRowSorter<>(this);

        sorter.setStringConverter(new TableStringConverter() {
            @Override
            public String toString(TableModel model, int row, int column) {
                final Object value = model.getValueAt(row, column);
                return value == null ? "" : value.toString();
            }
        });
    }

    public TableRowSorter<GenomeTableModel> getSorter() {
        return sorter;
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

        GenomeListItem record = records.get(rowIndex);
            String att = columnHeadings[columnIndex];
            return record.getAttributeValue(att);

    }

    public List<GenomeListItem> getRecords() {
        return records;
    }
}