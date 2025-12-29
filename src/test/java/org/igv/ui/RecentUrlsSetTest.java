package org.igv.ui;

import org.igv.util.ResourceLocator;
import org.junit.Assert;
import org.junit.Test;

import java.util.Iterator;
import java.util.List;


public class RecentUrlsSetTest {

    @Test
    public void testRoundTrip() {
        ResourceLocator withIndex1 = new ResourceLocator("path");
        withIndex1.setIndexPath("index");

        ResourceLocator withIndex2 = new ResourceLocator("i have an index");
        withIndex2.setIndexPath("yes i do");

        ResourceLocator noIndex1 = new ResourceLocator("path/no/index");
        ResourceLocator noIndex2 = new ResourceLocator("path 2");

        RecentUrlsSet resourceLocators = new RecentUrlsSet(List.of(withIndex1, withIndex2, noIndex1, noIndex2, withIndex1), 5);
        String serialized = resourceLocators.asString();
        Assert.assertEquals(serialized, "path index:index|i have an index index:yes i do|path/no/index|path 2");

        RecentUrlsSet roundTrip = RecentUrlsSet.fromString(serialized, 5);
        Assert.assertEquals(roundTrip.size(), resourceLocators.size());
        Iterator<ResourceLocator> actual = roundTrip.iterator();
        Iterator<ResourceLocator> expected = resourceLocators.iterator();
        while (actual.hasNext() && expected.hasNext()) {
            assertEqualLocator(expected.next(), actual.next());
        }
    }

    @Test
    public void testEmpty(){
        RecentUrlsSet empty = new RecentUrlsSet(10);
        Assert.assertEquals("", empty.asString());
        RecentUrlsSet fromEmptyString = RecentUrlsSet.fromString("", 5);
        Assert.assertEquals(0, fromEmptyString.size());
    }

    @Test
    public void testNull(){
        RecentUrlsSet empty = RecentUrlsSet.fromString(null, 5);
        Assert.assertEquals(empty.size(),0);
    }

    public static void assertEqualLocator(ResourceLocator expected, ResourceLocator actual){
        Assert.assertEquals(expected.getPath(), actual.getPath());
        Assert.assertEquals(expected.getIndexPath(), actual.getIndexPath());
    }


}