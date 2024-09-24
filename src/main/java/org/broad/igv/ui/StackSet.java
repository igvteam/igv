package org.broad.igv.ui;

import java.util.*;

/**
 *
 */
public class StackSet<T> extends AbstractCollection<T> implements Set<T>, SequencedCollection<T> {
    final private LinkedList<T> values;
    final int maxSize;

    public StackSet(int maxSize) {
        this(new LinkedList<>(), maxSize);
    }

    public StackSet(Collection<T> c, int maxSize) {
        this(new LinkedList<>(), maxSize);
        Collection<T> unique = new LinkedHashSet<>(c);
        List<T> limited = unique.stream().limit(maxSize).toList();
        values.addAll(limited);
    }

    private StackSet(LinkedList<T> values, int maxSize) {
        if (maxSize <= 0) {
            throw new IllegalArgumentException("maxSize must be > 0");
        }
        this.maxSize = maxSize;
        this.values = values;
    }

    @Override
    public boolean add(T t) {
        values.remove(t);
        if (values.size() >= maxSize) {
            values.removeLast();
        }
        values.addFirst(t);
        return true;
    }

    @Override
    public boolean addAll(Collection<? extends T> c) {
        ArrayList<? extends T> list = new ArrayList<>(c);
        list.reversed().forEach(this::add);
        return true;
    }

    @Override
    public boolean remove(Object o) {
        return values.remove(o);
    }

    @Override
    public void clear() {
        values.clear();
    }

    @Override
    public Iterator<T> iterator() {
        return values.iterator();
    }

    @Override
    public int size() {
        return values.size();
    }

    @Override
    public StackSet<T> reversed() {
        return new StackSet<>(values.reversed(), maxSize);
    }

}


