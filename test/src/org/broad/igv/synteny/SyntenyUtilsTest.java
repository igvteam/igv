/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.synteny;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import org.junit.Test;

import java.util.List;
import java.util.Map;

/**
 * User: jrobinso
 * Date: Feb 19, 2010
 */
public class SyntenyUtilsTest {

    String testMappings = "test/data/hg18_to_mm8.regions";

    @Test
    public void testLoadMappings() {

        //region R:chr2:chr2:D10 chr2 139839284 176464687 + chr2 39414197 74293149 +

        // Test forward mapping
        String humanChr = "chr2";
        int humanPosition = 139839284;

        String mouseChr = "chr2";
        int mousePosition = 39414197;

        Map<String, List<SyntenyMapping>> mappings = SyntenyUtils.loadMappings(testMappings, false);
        SyntenyMapping mapping = SyntenyUtils.getMappingContaining(mappings.get(humanChr), humanPosition);
        assertNotNull(mapping);

        double toPosition = mapping.mapPosition(humanPosition);
        assertEquals("Forward mapping", mousePosition, toPosition, 1.0);

        mappings = SyntenyUtils.loadMappings(testMappings, true);
        mapping = SyntenyUtils.getMappingContaining(mappings.get(mouseChr), mousePosition);
        assertNotNull(mapping);

        toPosition = mapping.mapPosition(mousePosition);
        assertEquals("Forward mapping", humanPosition, toPosition, 1.0);

    }
}
