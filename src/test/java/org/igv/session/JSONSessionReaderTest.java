package org.igv.session;

import org.igv.data.TestDataSource;
import org.igv.renderer.DataRange;
import org.igv.track.DataSourceTrack;
import org.igv.track.WindowFunction;
import org.igv.util.TestUtils;
import org.json.JSONObject;
import org.junit.Test;

import java.nio.file.Files;
import java.nio.file.Paths;

import static org.junit.Assert.*;

public class JSONSessionReaderTest {



    @Test

    /**
     * Test unmarshalling the following wig track
     *
     * {
     *   "name": "Homo sapiens GM12878 ATF2",
     *   "url": "https://www.encodeproject.org/files/ENCFF126YRP/@@download/ENCFF126YRP.bigWig",
     *   "metadata": {
     *     "Biosample": "Homo sapiens GM12878",
     *     "AssayType": "ChIP-seq",
     *     "Target": "ATF2",
     *     "BioRep": "2",
     *     "TechRep": "2_1",
     *     "OutputType": "signal p-value",
     *     "Format": "bigWig",
     *     "Lab": "Michael Snyder, Stanford",
     *     "Accession": "ENCFF126YRP",
     *     "Experiment": "ENCSR961PPA"
     *   },
     *   "order": 6,
     *   "format": "bigwig",
     *   "type": "wig",
     *   "color": "rgb(69, 132, 1)",
     *   "logScale": true,
     *   "min": 0,
     *   "max": 100,
     *   "windowFunction": "max"
     * }
     */
    public void testUnmarshalWigTrack() throws Exception {
        String jsonFilePath = TestUtils.DATA_DIR + "sessions/json/wigTrack.json";

        byte[] bytes = Files.readAllBytes(Paths.get(jsonFilePath));
        String jsonString = new String(bytes);
        JSONObject jsonObject = new JSONObject(jsonString);
        assertNotNull(jsonObject);

        DataSourceTrack track = new DataSourceTrack();
        track.setDatasource(new TestDataSource());
        track.unmarshalJSON(jsonObject);

        assertEquals("Homo sapiens GM12878 ATF2", track.getName());
        assertEquals(WindowFunction.max, track.getWindowFunction());

        DataRange dataRange = track.getDataRange();
        assertEquals(100, dataRange.getMaximum(), 0.00000001f);
        assertTrue(dataRange.isLog());
    }
}