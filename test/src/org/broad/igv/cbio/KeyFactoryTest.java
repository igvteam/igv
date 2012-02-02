/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *  
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *  
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.cbio;

import org.junit.Test;

import static org.junit.Assert.*;

/**
 * User: jacob
 * Date: 2012/02/02
 */
public class KeyFactoryTest {

    public KeyFactoryTest() {
    }

    @Test
    public void testNoDuplicateKeys() throws Exception {
        //We store each key that is generated
        //to create the schema. Test that this is 
        //done properly
        String[] keys = "aoeu,wj,dulllife,yeahyeahyeah,mapswait,skippydoo,".split(",");
        KeyFactory factory = new KeyFactory("test");
        int trials = 10;
        for (int ii = 0; ii < trials; ii++) {
            for (String key : keys) {
                factory.getDataKey(key, key.hashCode());
            }

            assertEquals(keys.length, factory.getKeySet().size());
        }
    }

    @Test
    public void testValidTypes() throws Exception {
        KeyFactory factory = new KeyFactory("test");
        factory.getDataKey("int", 1);
        factory.getDataKey("long", 12345678923748L);
        factory.getDataKey("double", 1.01236548d);
        factory.getDataKey("float", 2.4e6f);
        factory.getDataKey("boolean", true);
        factory.getDataKey("string", "teststring");
        assertEquals(6, factory.getKeySet().size());
    }

    @Test
    public void testInvalidTypes() throws Exception {
        boolean ex = false;
        KeyFactory factory = new KeyFactory("test");
        try {
            factory.getDataKey("char", 'a');
        } catch (IllegalArgumentException e) {
            ex = true;
        }
        assertTrue("Factory did not throw exception", ex);

        ex = false;
        try {
            factory.getDataKey("short", (short) 5613);
        } catch (IllegalArgumentException e) {
            ex = true;
        }
        assertTrue("Factory did not throw exception", ex);

        ex = false;
        try {
            factory.getDataKey("byte", (byte) 2);
        } catch (IllegalArgumentException e) {
            ex = true;
        }
        assertTrue("Factory did not throw exception", ex);

        assertEquals(0, factory.getKeySet().size());
    }

    @Test
    public void testDataKeyEquals() {
        DataKey k1 = new DataKey("test", "double");
        DataKey k2 = new DataKey("test", "string");
        assertEquals(k1, k2);

        DataKey k3 = new DataKey("test1", "double");
        assertNotSame(k1, k3);
        assertFalse(k1.equals(k3));

        assertEquals(k1.equals(k2), k2.equals(k1));
        assertEquals(k1.hashCode(), k2.hashCode());
    }
}
