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

package org.broad.igv.encode;

import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;
import javax.swing.table.TableStringConverter;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jrobinso
 *         Date: 10/31/13
 *         Time: 10:09 PM
 */
public class EncodeTableModel extends AbstractTableModel {

    private String[] columnHeadings;
    private List<EncodeFileRecord> records;
    private final TableRowSorter<EncodeTableModel> sorter;

    public EncodeTableModel(List<String> headings, List<EncodeFileRecord> records) {

        this.records = records;

        List<String> tmp = new ArrayList<String>();
        tmp.add("");  // Checkbox heading
        for(String h : headings) {
            String heading = h.trim();
            if(heading.length() > 0 && !"path".equals(heading)) {
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

    public TableRowSorter<EncodeTableModel> getSorter() {
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

        EncodeFileRecord record = records.get(rowIndex);
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
        if(col == 0) {
            records.get(row).setSelected((Boolean) value);
        }
        fireTableCellUpdated(row, col);
    }



    public List<EncodeFileRecord> getRecords() {
        return records;
    }
}