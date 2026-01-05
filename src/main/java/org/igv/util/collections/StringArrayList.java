
package org.igv.util.collections;

/**
 * @author jrobinso
 */
public class StringArrayList {

    static final int maxGrowIncrement = Integer.MAX_VALUE / 10;

    int size = 0;
    String[] values;

    public StringArrayList(int maxSize) {
        values = new String[maxSize];
    }

    public void add(String v) {
        if (size >= values.length) {
            grow();
        }
        values[size] = v;
        size++;
    }

       public String get(int idx) {
        if (idx >= size) {
            throw new IndexOutOfBoundsException("Index out of bounds.  Size = " + size + " Index = " + idx);
        } else {
            return values[idx];
        }
    }

    public String[] toArray() {
        trim();
        return values;
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
        String[] tmp = new String[newSize];
        System.arraycopy(values, 0, tmp, 0, Math.min(tmp.length, values.length));
        values = tmp;
    }

    private void trim() {
        resize(size);
    }

}
