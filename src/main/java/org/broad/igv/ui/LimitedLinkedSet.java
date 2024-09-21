package org.broad.igv.ui;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 *
 */
public class LimitedLinkedSet<T> extends AbstractCollection<T>{
    final protected LinkedList<T> values = new LinkedList<>();
    final int maxSize;

    public LimitedLinkedSet(int maxSize) {
        if (maxSize <= 0) {
            throw new IllegalArgumentException("maxSize must be > 0");
        }
        this.maxSize = maxSize;
    }

    public LimitedLinkedSet(Collection<T> c, int maxSize) {
        this(maxSize);
        Collection<T> unique = new LinkedHashSet<>(c);
        List<T> limited = unique.stream().limit(maxSize).toList();
        values.addAll(limited);
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

}


