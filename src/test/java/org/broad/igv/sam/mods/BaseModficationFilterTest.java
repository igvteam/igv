package org.broad.igv.sam.mods;

import org.junit.Test;

import static org.junit.Assert.*;

public class BaseModficationFilterTest {

    @Test
    public void pass() {
        BaseModficationFilter filter = new BaseModficationFilter(null, 'C');
        assertTrue(filter.pass(null, 'C'));

        String nomod = "NONE_C";
        assertTrue(filter.pass(nomod));
    }

    @Test
    public void fromString() {

        String str = "hm,C";
        BaseModficationFilter filter = BaseModficationFilter.fromString(str);

        assertEquals("hm", filter.modification);
        assertEquals('C', filter.base);

        assertTrue(filter.pass("m"));
        assertTrue(filter.pass("h"));
        assertTrue(filter.pass("NONE_C"));
        assertFalse(filter.pass("a"));
    }
}