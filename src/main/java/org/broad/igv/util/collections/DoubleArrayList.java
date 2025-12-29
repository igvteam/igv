package org.broad.igv.util.collections;

import java.util.Arrays;

/**
 * Author: jrobinso
 * Date: Jul 22, 2010
 * <p/>
 * ArrayList type collection for int types.  Purpose is to avoid the need to create an object for each entry
 * in the standard java collections.
 */
public class DoubleArrayList {


    private transient double[] elements;

    private int size;


    public DoubleArrayList() {
        this(100);
    }

    public DoubleArrayList(int initialCapacity) {
        if (initialCapacity < 0)
            throw new IllegalArgumentException("Illegal Capacity: " + initialCapacity);
        this.elements = new double[initialCapacity];
    }

    public DoubleArrayList(double[] elements) {
        this.elements = elements;
        size = elements.length;
    }

    public void add(double e) {
        if (size + 1 >= elements.length) {
            grow();
        }
        elements[size++] = e;
    }

    public void addAll(double[] args) {
        double[] newElements = new double[size + args.length];
        System.arraycopy(elements, 0, newElements, 0, size);
        System.arraycopy(args, 0, newElements, size, args.length);
        elements = newElements;
        size += args.length;
    }

    public void addAll(DoubleArrayList aList) {
        addAll(aList.toArray());
    }


    public double get(int idx) {
        return elements[idx];
    }


    public void set(int idx, double d) {

        if(idx >= size) {
            throw new IndexOutOfBoundsException("Index " + idx + " out of bounds.  Size = " + size);
        }
        elements[idx] = d;

    }

    public int size() {
        return size;
    }

    public boolean isEmpty() {
        return size == 0;
    }

    /**
     * Empty all elements.  This logically clears the collection but does not free up any space.
     */
    public void clear() {
        size = 0;
    }

    private void grow() {
        int oldCapacity = elements.length;
        int newCapacity;
        if (oldCapacity < 10000000) {
            newCapacity = oldCapacity * 2;
        } else {
            newCapacity = (oldCapacity * 3) / 2 + 1;
        }
        double[] tmp = new double[newCapacity];
        System.arraycopy(elements, 0, tmp, 0, elements.length);
        elements = tmp;
    }


    public double[] toArray() {
        trimToSize();
        return elements;
    }

    public double [] toArray(int start, int end) {

        if(end > size) throw new IndexOutOfBoundsException();

        return Arrays.copyOfRange(elements, start, end);


    }

    private void trimToSize() {
        int oldCapacity = elements.length;
        if (size < oldCapacity) {
            double[] tmp = new double[size];
            System.arraycopy(elements, 0, tmp, 0, size);
            elements = tmp;
        }
    }

}
