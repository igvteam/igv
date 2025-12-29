package org.igv.ucsc.hub;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import javax.swing.table.TableStringConverter;
import java.util.List;

/**
 * @author jrobinso
 *
 */
public class HubTableModel extends AbstractTableModel {

    private String[] columnHeadings;
    private List<HubDescriptor> records;

    public HubTableModel(List<HubDescriptor> records) {

        this.records = records;

        //tmp.add("path");
        columnHeadings = new String[]{"", "Name", "Description", "DB List", "Info"};

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
            return ""; // Return an empty string for invalid indices to avoid null issues
        }

        HubDescriptor record = records.get(rowIndex);
        switch (columnIndex) {
            case 0:
                return record.isSelected(); // Boolean for sorting checkboxes
            case 1:
                return record.getShortLabel(); // String
            case 2:
                return record.getLongLabel(); // String
            case 3:
                return record.getDbList(); // String
            case 4:
                return record.getDescriptionUrl(); // String
            default:
                return ""; // Default to an empty string for unsupported columns
        }
    }

    @Override
    public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
        if(columnIndex == 0) {
            records.get(rowIndex).setSelected((Boolean) aValue);
        }
    }

    public List<HubDescriptor> getRecords() {
        return records;
    }


}