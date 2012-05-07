/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.synteny;

import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

/**
 * User: jrobinso
 * Date: Feb 19, 2010
 */
public class SyntenyUtilsTest {

    String testMappings = TestUtils.DATA_DIR + "igv/hg18_to_mm8.regions";

    @Test
    public void testLoadMappings() throws IOException {

        //region R:chr2:chr2:D10 chr2 139839284 176464687 + chr2 39414197 74293149 +

        // Test forward mapping
        String humanChr = "chr2";
        int humanPosition = 139839284;

        String mouseChr = "chr2";
        int mousePosition = 39414197;


        Map<String, List<Mapping>> mappings = SyntenyUtils.loadMappings(testMappings, false);
        Mapping mapping = SyntenyUtils.getMappingContaining(mappings.get(humanChr), humanPosition);
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
