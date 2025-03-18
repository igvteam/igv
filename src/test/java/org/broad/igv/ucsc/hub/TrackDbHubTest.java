package org.broad.igv.ucsc.hub;

import org.broad.igv.feature.genome.load.TrackConfig;
import org.broad.igv.util.TestUtils;
import org.junit.Test;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

public class TrackDbHubTest {

    //metadata differentiation="10 hour" treatment=X donor=A lab="List Meta Lab" data_set_id=ucscTest1 access=group assay=long-RNA-seq enriched_in=exon life_stage=postpartum species="Homo sapiens" ucsc_db=hg38

    @Test
    public void parseMetadata() {

        String metadata = "differentiation=\"10 hour\" treatment=X donor=A lab=\"List Meta Lab\" data_set_id=ucscTest1 access=group assay=long-RNA-seq enriched_in=exon life_stage=postpartum species=\"Homo sapiens\" ucsc_db=hg38";
        Map<String, String> metadataMap = TrackDbHub.parseMetadata(metadata);
        assertEquals(11, metadataMap.size());
        assertEquals("X", metadataMap.get("Treatment"));
        assertEquals("10 hour", metadataMap.get("Differentiation"));
        assertEquals("A", metadataMap.get("Donor"));
        assertEquals("group", metadataMap.get("Access"));
        assertEquals("hg38", metadataMap.get("Ucsc_db"));
        assertEquals("Homo sapiens", metadataMap.get("Species"));
        assertEquals("postpartum", metadataMap.get("Life_stage"));
        assertEquals("long-RNA-seq", metadataMap.get("Assay"));
        assertNotNull(metadataMap);
    }

    //        metadata "Epigenome_Mnemonic"="GI.STMC.GAST" "Standardized_Epigenome_name"="Gastric" "EDACC_Epigenome_name"="Gastric" "Group"="<span style="color:#C58DAA">Digestive</span>" "Age"="34Y" "Lab"="UCSD" "Sex"="Male" "Anatomy"="GI_STOMACH" "EID"="E094" "Type"="PrimaryTissue" "Order"="100" "Ethnicity"="Caucasian"

    @Test
    public void parseMetadata2() {

        String metadata = "\"Epigenome_Mnemonic\"=\"GI.STMC.GAST\" \"Standardized_Epigenome_name\"=\"Gastric\" \"EDACC_Epigenome_name\"=\"Gastric\" \"Group\"=\"<span style=\"color:#C58DAA\">Digestive</span>\" \"Age\"=\"34Y\" \"Lab\"=\"UCSD\" \"Sex\"=\"Male\" \"Anatomy\"=\"GI_STOMACH\" \"EID\"=\"E094\" \"Type\"=\"PrimaryTissue\" \"Order\"=\"100\" \"Ethnicity\"=\"Caucasian\"";
        Map<String, String> metadataMap = TrackDbHub.parseMetadata(metadata);
        assertEquals("GI.STMC.GAST", metadataMap.get("Epigenome_Mnemonic"));
        assertEquals("Digestive", metadataMap.get("Group"));
        assertEquals("Caucasian", metadataMap.get("Ethnicity"));

        assertNotNull(metadataMap);
    }

    /**
     * group roadmap
     * .. track CompRoadmapbySample
     * .... track CompRoadmap_HUES64
     * ...... track CompRoadmap_HUES64ViewCOV
     * ........ track ENCFF002FAZ_HUES64_ChIP-seq_bySample
     * ........ track ENCFF021OWW_HUES64_ChIP-seq_bySample
     *
     * @throws IOException
     */
    @Test
    public void testSuperTrack() throws IOException {

        List<Stanza> stanzas = HubParser.loadStanzas(TestUtils.DATA_DIR + "hubs/supertrack_trackDb.txt");

        TrackDbHub trackDbHub = new TrackDbHub(stanzas, null);

        List<TrackConfigContainer> containers = trackDbHub.getGroupedTrackConfigurations("test");
        assertEquals(1, containers.size());

        TrackConfigContainer superTrack = containers.get(0);
        assertEquals("CompRoadmapbySample", superTrack.name);
        assertEquals(1, superTrack.children.size());

        TrackConfigContainer compositeTrack = superTrack.children.get(0);
        assertEquals("HUES64 tracks", compositeTrack.label);
        assertEquals(1, compositeTrack.children.size());

        TrackConfigContainer view = compositeTrack.children.get(0);
        assertEquals("Coverage", view.label);
        assertEquals(0, view.children.size());
        assertEquals(2, view.tracks.size());

        TrackConfig track = view.tracks.get(0);
        Map<String, String> meta = track.getAttributes();
        assertEquals(4, meta.size());
    }

    @Test
    public void testGroupedContainers() throws IOException {

        List<Stanza> stanzas = HubParser.loadStanzas(TestUtils.DATA_DIR + "hubs/gtexCoverage.txt");
        List<Stanza> groupStanzas = HubParser.loadStanzas(TestUtils.DATA_DIR + "hubs/groups.txt");
        TrackDbHub trackDbHub = new TrackDbHub(stanzas, groupStanzas);
        List<TrackConfigContainer> containers = trackDbHub.getGroupedTrackConfigurations("test");
        assertEquals(1, containers.size());

        TrackConfigContainer groupContainer = containers.get(0);
        assertEquals("expression", groupContainer.name);
        assertEquals(1, groupContainer.children.size());
        assertEquals(0, groupContainer.tracks.size());

        TrackConfigContainer compositeContainer = groupContainer.children.get(0);
        assertEquals("gtexCov", compositeContainer.name);
        assertEquals(0, compositeContainer.children.size());
        assertEquals(2, compositeContainer.tracks.size());

    }
}