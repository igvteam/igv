package org.igv.ui;

import java.util.*;

/**
 * This is a very specific collection which is designed to support "recent item" lists.
 *
 * It behaves like a stack, the most recently added item is always the first, second most recent is the second, and so forth
 * It disallows duplicate items. Adding a duplicate item will move it to the top of the list.
 * There is also size limit, adding additional items beyond the maximum size will remove the oldest items from the collection
 * to make room.
 *
 * @implNote This is implemented in an extremely inefficient way, do not use this for large collections.
 */
public class StackSet<T> extends AbstractCollection<T> implements Set<T>, SequencedCollection<T> {
    final private LinkedList<T> values;
    final int maxSize;

    /**
     * Create an empty StackSet with maxSize
     */
    public StackSet(int maxSize) {
        this(new LinkedList<>(), maxSize);
    }

    /**
     * Constructs a StackSet containing the elements of the specified collection, in the order
     * they are returned by the collection's iterator.  (The first element returned by the collection's
     * iterator becomes the first element, or <i>top</i> of the stack.)
     *
     * Duplicate values appear at their earliest position.
     *
     * This order consistent with the way {@link ArrayDeque} is implemented;
     *
     * @param initialValues
     * @param maxSize
     */
    public StackSet(Collection<T> initialValues, int maxSize) {
        this(new LinkedList<>(), maxSize);
        Collection<T> unique = new LinkedHashSet<>(initialValues);
        List<T> limited = unique.stream().limit(maxSize).toList();
        values.addAll(limited);
    }


    /**
     * private constructor allowing direct setting of the values in order to support reverse views
     */
    private StackSet(LinkedList<T> values, int maxSize) {
        if (maxSize <= 0) {
            throw new IllegalArgumentException("maxSize must be > 0");
        }
        this.maxSize = maxSize;
        this.values = values;
    }

    /**
     * Add an element to the top of the stack.  If it is already present it will be
     * moved to the top.
     * @param t element to add
     * @return always true
     */
    @Override
    public boolean add(T t) {
        values.remove(t);
        if (values.size() >= maxSize) {
            values.removeLast();
        }
        values.addFirst(t);
        return true;
    }

    /**
     * Add the elements of this collection to the top of the stack in the order of the collections
     * iterator.
     *
     * This is different from the implementation of {@link Deque} which adds to the end of the queue
     * @param c collection containing elements to be added to this collection
     * @return always true
     */
    @Override
    public boolean addAll(Collection<? extends T> c) {
        ArrayList<? extends T> list = new ArrayList<>(c);
        list.reversed().forEach(this::add);
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean remove(Object o) {
        return values.remove(o);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void clear() {
        values.clear();
    }

    /**
     * @return an iterator over the elements of the Stack from most recently to oldest
     */
    @Override
    public Iterator<T> iterator() {
        return values.iterator();
    }

    /**
    * {@inheritDoc}
    */
    @Override
    public int size() {
        return values.size();
    }

    /**
     * @return the maximum allowed size for this collection
     */
    public int getMaxSize(){
        return maxSize;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public StackSet<T> reversed() {
        return new StackSet<>(values.reversed(), maxSize);
    }

}


