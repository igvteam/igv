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

package org.broad.igv.ucsc.hub;

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