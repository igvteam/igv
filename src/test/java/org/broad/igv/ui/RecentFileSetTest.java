package org.broad.igv.ui;

import org.junit.Assert;
import org.junit.Test;

import java.util.List;

public class RecentFileSetTest {

    @Test
    public void testAdds() {
        RecentFileSet set1 = new RecentFileSet(3);
        set1.add("a");
        set1.add("b");
        set1.add("c");
        set1.add("d");
        set1.add("e");
        LimitedLinkedSetTest.assertEquals(set1, List.of("e","d","c"));
    }

    @Test
    public void testAsListRoundTrip(){
        RecentFileSet set = new RecentFileSet(List.of("a", "b", "c", "d", "b"), 3);
        LimitedLinkedSetTest.assertEquals(set, List.of("a", "b", "c"));
        String string = set.asString();
        Assert.assertEquals("a;b;c", string);
        RecentFileSet set2 = RecentFileSet.fromString(string, 5);
        LimitedLinkedSetTest.assertEquals(set2, List.of("a", "b", "c"));
    }

    @Test
    public void testNullString(){
        RecentFileSet set = RecentFileSet.fromString(null, 5);
        LimitedLinkedSetTest.assertEquals(set, List.of());
    }

    @Test
    public void testEmptyString(){
        RecentFileSet set = RecentFileSet.fromString("", 5);
        LimitedLinkedSetTest.assertEquals(set, List.of());
    }

    @Test
    public void testWhiteSpaceString(){
        RecentFileSet set = RecentFileSet.fromString(" a; b ; c; ", 5);
        LimitedLinkedSetTest.assertEquals(set, List.of("a", "b", "c"));
    }

    @Test
    public void testInternalWhiteSpaceString(){
        RecentFileSet set = RecentFileSet.fromString("this file has spaces;thisonedoesnt", 5);
        LimitedLinkedSetTest.assertEquals(set, List.of("this file has spaces", "thisonedoesnt"));
    }


}