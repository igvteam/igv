/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.broad.igv.feature;

import htsjdk.tribble.NamedFeature;
import org.broad.igv.AbstractHeadlessTest;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

import static junit.framework.Assert.assertNotNull;

/**
 * @author jrobinso
 */
public class TestLoadGeneManager extends AbstractHeadlessTest {

    @Test
    public void main() throws IOException {
        NamedFeature feature = genome.getFeatureDB().getFeature("EGFR");
        assertNotNull(feature);
    }

}
