package org.igv.ucsc.hub;

import org.igv.feature.genome.load.TrackConfig;
import org.igv.prefs.Constants;
import org.igv.prefs.PreferencesManager;
import org.igv.util.TestUtils;
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
     * Hub heights below the minimum track height preference are raised to it; larger heights are left alone.
     */
    @Test
    public void testMinimumHeight() throws IOException {

        TestUtils.resetPrefsFile();
        int minimumHeight = PreferencesManager.getPreferences().getAsInt(Constants.TRACK_HEIGHT_KEY);
        assertTrue(minimumHeight > 12);

        List<Stanza> stanzas = HubParser.loadStanzas(TestUtils.DATA_DIR + "hubs/min_height_trackDb.txt");
        TrackDbHub trackDbHub = new TrackDbHub(stanzas, null);
        List<TrackConfig> tracks = trackDbHub.getGroupedTrackConfigurations("test").get(0).tracks;
        assertEquals(2, tracks.size());

        TrackConfig small = tracks.get(0);
        assertEquals(64, (int) small.maxHeight);
        assertEquals(minimumHeight, (int) small.height);
        assertEquals(minimumHeight, (int) small.minHeight);

        TrackConfig large = tracks.get(1);
        assertEquals(50, (int) large.height);
        assertEquals(30, (int) large.minHeight);
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

        // No groups.txt, so "group roadmap" produces an implicit top-level container
        TrackConfigContainer group = containers.get(0);
        assertEquals("roadmap", group.name);
        assertEquals(1, group.children.size());

        TrackConfigContainer superTrack = group.children.get(0);
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
        Map<String, String> meta = track.attributes;
        assertEquals(4, meta.size());
    }

    @Test
    public void testUnsupportedTrackCount() throws IOException {

        // groups_trackDb.txt's only data tracks are two halSnake tracks, a type IGV does not support, so both are
        // counted and nothing is left to show
        List<Stanza> stanzas = HubParser.loadStanzas(TestUtils.DATA_DIR + "hubs/groups_trackDb.txt");
        List<Stanza> groupStanzas = HubParser.loadStanzas(TestUtils.DATA_DIR + "hubs/groups.txt");
        TrackDbHub trackDbHub = new TrackDbHub(stanzas, groupStanzas);
        List<TrackConfigContainer> containers = trackDbHub.getGroupedTrackConfigurations("test");
        assertEquals(2, trackDbHub.getUnsupportedTrackCount());
        assertTrue(containers.isEmpty());
    }

    /**
     * A data track that also declares "view" (the Weng lab ENCODE hubs do this on every data track) is a track, not
     * a view container.  Treating it as a container left the whole hub empty.
     */
    @Test
    public void testViewOnDataTrack() throws IOException {

        List<Stanza> stanzas = HubParser.loadStanzas(TestUtils.DATA_DIR + "hubs/view_on_datatrack_trackDb.txt");
        TrackDbHub trackDbHub = new TrackDbHub(stanzas, null);
        List<TrackConfigContainer> containers = trackDbHub.getGroupedTrackConfigurations("test");
        assertEquals(0, trackDbHub.getUnsupportedTrackCount());
        assertEquals(1, containers.size());

        TrackConfigContainer composite = containers.get(0).children.get(0);
        assertEquals("largeIntestine", composite.name);
        // The bigWig with "view bigWig" sits directly in the composite
        assertEquals(1, composite.tracks.size());
        assertEquals("ENCFF514ZYW", composite.tracks.get(0).id);
        // The bigBed sits in its real view container
        assertEquals(1, composite.children.size());
        assertEquals("ENCSR091QXPview", composite.children.get(0).name);
        assertEquals(1, composite.children.get(0).tracks.size());
        assertEquals("ENCFF948PMS", composite.children.get(0).tracks.get(0).id);
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