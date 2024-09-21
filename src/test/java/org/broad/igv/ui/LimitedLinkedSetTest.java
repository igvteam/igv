package org.broad.igv.ui;

import org.junit.Assert;
import org.junit.Test;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

public class LimitedLinkedSetTest {


    @Test
    public void testAdds(){
        LimitedLinkedSet<Integer> set = new LimitedLinkedSet<>(5);
        set.add(1);
        assertEquals(set, List.of(1));
        set.add(2);
        assertEquals(set, List.of(2, 1));
        set.add(3);
        assertEquals(set, List.of(3, 2, 1));
        set.add(4);
        assertEquals(set, List.of(4, 3, 2, 1));
        set.add(5);
        assertEquals(set, List.of(5, 4, 3, 2, 1));
        set.add(6);
        assertEquals(set, List.of(6, 5, 4, 3, 2));
        set.add(7);
        assertEquals(set, List.of(7, 6, 5, 4, 3));
        set.add(7);
        assertEquals(set, List.of(7, 6, 5, 4, 3));
        set.add(1);
        assertEquals(set, List.of(1, 7, 6, 5, 4));
        set.add(7);
        assertEquals(set, List.of(7, 1, 6, 5, 4));
        set.remove(1);
        assertEquals(set, List.of(7, 6, 5, 4));
        set.add(7);
        assertEquals(set, List.of(7, 6, 5, 4));
    }

    @Test
    public void testNewCollection(){
        LimitedLinkedSet<Integer> set = new LimitedLinkedSet<>(List.of(1, 1, 2, 3, 4, 5, 6, 7, 3), 5);
        assertEquals(set, List.of(1,2,3,4,5));
        set = new LimitedLinkedSet<>(List.of(1,1,1,1,1,1,1,1,1,1), 10);
        assertEquals(set, List.of(1));
        set = new LimitedLinkedSet<>(Collections.emptySet(), 1);
        assertEquals(set, List.of());
        set.add(1);
        assertEquals(set, List.of(1));
        set.add(2);
        assertEquals(set, List.of(2));
    }

    public static <T> void assertEquals(Collection<T> actual, List<T> expected){
        Assert.assertEquals(expected.size(), actual.size());
        Assert.assertArrayEquals(expected.toArray(), actual.toArray());
    }
}