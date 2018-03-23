/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
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
