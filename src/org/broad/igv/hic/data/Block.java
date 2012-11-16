package org.broad.igv.hic.data;

import java.awt.*;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * @author jrobinso
 * @date Aug 10, 2010
 */
public class Block {

    private int number;

   ContactRecord[] records;


    public Block(int number, ContactRecord[] records) {
        this.number = number;
        this.records = records;
    }

    public int getNumber() {
        return number;
    }


    public ContactRecord[] getContactRecords() {
        return records;
    }

}
