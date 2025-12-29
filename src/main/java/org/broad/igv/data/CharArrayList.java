/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.data;

/**
 * @author jrobinso
 */
public class CharArrayList {

    static final int maxGrowIncrement = Integer.MAX_VALUE / 10;

    int size = 0;
    char[] values;

    public CharArrayList(int maxSize) {
        values = new char[maxSize];
    }

    public void add(char v) {
        if (size >= values.length) {
            grow();
        }
        values[size] = v;
        size++;
    }

    public char[] toArray() {
        trim();
        return values;
    }

    public int getSize() {
        return size;
    }

    private void grow() {
        if (values.length >= Integer.MAX_VALUE) {
            throw new RuntimeException("Maximum array size exceeded");
        }
        int increment = (int) (Math.max(1000, 0.2 * values.length));
        int newSize = Math.min(Integer.MAX_VALUE, values.length + increment);
        resize(newSize);

    }

    private void resize(int newSize) {
        char[] tmp = new char[newSize];
        System.arraycopy(values, 0, tmp, 0, Math.min(tmp.length, values.length));
        values = tmp;
    }

    private void trim() {
        resize(size);
    }

}
