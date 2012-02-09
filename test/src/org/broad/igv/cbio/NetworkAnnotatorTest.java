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

import org.apache.commons.collections.Predicate;
import org.broad.igv.util.TestUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/02/09
 */
public class NetworkAnnotatorTest {

    private String localpath = TestUtils.DATA_DIR + "/tp53network.xml";
    private NetworkAnnotator annotator;

    @Before
    public void setUp() {
        annotator = new NetworkAnnotator();
    }

    @After
    public void tearDown() {
        annotator = null;
    }

    @Test
    public void testLoadLocal() throws Exception {
        assertTrue("Failed to load network", annotator.loadNetwork(localpath));
    }

    @Test
    public void testFilter() throws Exception {
        Predicate tPred = new Predicate() {

            public boolean evaluate(Object object) {
                Node node = (Node) object;
                NamedNodeMap map = node.getAttributes();
                if (map == null) {
                    return false;
                }
                int id = Integer.parseInt(map.getNamedItem("id").getTextContent());
                return id % 2 == 0;
            }
        };

        annotator.loadNetwork(localpath);
        int removed = annotator.filterNodes(tPred);
        assertTrue(removed > 0);
    }


}
