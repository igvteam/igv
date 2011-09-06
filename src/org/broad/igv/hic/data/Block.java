package org.broad.igv.hic.data;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Aug 10, 2010
 */
public class Block {

    int number;

    // Temporary map used during parsing only TODO -- remove from this class
    Map<Point, ContactRecord> contactRecordMap;

    ContactRecord[] records;

    public Block(int number) {
        this.number = number;
        contactRecordMap = new HashMap();
    }

    public Block(int number, ContactRecord[] records) {
        this.number = number;
        this.records = records;
    }

    public void incrementCount(int col, int row) {
        Point p = new Point(col, row);
        ContactRecord rec = contactRecordMap.get(p);
        if (rec == null) {
            rec = new ContactRecord(number, col, row, (short) 1);
            contactRecordMap.put(p, rec);
        } else {
            rec.incrementCount();
        }
    }

    public void parsingComplete() {
        if (contactRecordMap.size() > 0) {
            records = new ContactRecord[contactRecordMap.size()];
            int i = 0;
            for (ContactRecord rec : contactRecordMap.values()) {
                records[i] = rec;
                i++;
            }
        }
    }


    public ContactRecord[] getContactRecords() {
        return records;
    }
}
