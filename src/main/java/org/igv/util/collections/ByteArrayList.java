package org.igv.util.collections;


/**
 * Author: jrobinso
 * Date: Jul 22, 2010
 * <p/>
 * ArrayList type collection for int types.  Purpose is to avoid the need to create an object for each entry
 * in the standard java collections.
 */
public class ByteArrayList {


    private transient byte[] elements;

    private int size;


    public ByteArrayList() {
        this(100);
    }

    public ByteArrayList(int initialCapacity) {
        if (initialCapacity < 0)
            throw new IllegalArgumentException("Illegal Capacity: " + initialCapacity);
        this.elements = new byte[initialCapacity];
    }

    public ByteArrayList(byte[] elements) {
        this.elements = elements;
        size = elements.length;
    }

    public void add(byte e) {
        if (size + 1 >= elements.length) {
            grow();
        }
        elements[size++] = e;
    }

    public void addAll(byte[] args) {
        byte[] newElements = new byte[size + args.length];
        System.arraycopy(elements, 0, newElements, 0, size);
        System.arraycopy(args, 0, newElements, size, args.length);
        elements = newElements;
        size += args.length;
    }

    public void addAll(ByteArrayList aList) {
        addAll(aList.toArray());
    }

    public byte get(int idx) {
        return elements[idx];
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
        byte[] tmp = new byte[newCapacity];
        System.arraycopy(elements, 0, tmp, 0, elements.length);
        elements = tmp;
    }


    public byte[] toArray() {
        trimToSize();
        return elements;
    }


    private void trimToSize() {
        int oldCapacity = elements.length;
        if (size < oldCapacity) {
            byte[] tmp = new byte[size];
            System.arraycopy(elements, 0, tmp, 0, size);
            elements = tmp;
        }
    }

    public void set(int idx, byte i) {
        while(idx >= elements.length) {
            grow();
        }
        elements[idx] = i;
        idx++;
        if(idx > size) size = idx;  // Tried Math.max here, it showed up in cpu profiles!
    }
}
