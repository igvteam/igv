package org.broad.igv.feature.genome.load;

import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.Track;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

public class JsonGenomeLoaderTest {

    @Test
    public void loadGenome() throws IOException {
        String path = TestUtils.DATA_DIR + "genomes/json/hg38.json";
        JsonGenomeLoader loader = new JsonGenomeLoader(path);
        Genome genome = loader.loadGenome();

        assertEquals("hg38", genome.getId());
        assertEquals("Human (GRCh38/hg38)", genome.getDisplayName());

        List<ResourceLocator> annotationResources = genome.getAnnotationResources();
        assertEquals(1, annotationResources.size());
        ResourceLocator loc = annotationResources.get(0);
        assertEquals("Refseq Genes", loc.getName());
        assertEquals("https://s3.amazonaws.com/igv.org.genomes/hg38/ncbiRefSeq.sorted.txt.gz", loc.getPath());
        assertTrue( -1 ==  loc.getVisibilityWindow());
    }

    @Test
    public void loadDescriptor() {
    }
}
