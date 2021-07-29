package org.broad.igv.feature.genome.load;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import org.apache.log4j.Logger;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.fasta.FastaBlockCompressedSequence;
import org.broad.igv.feature.genome.fasta.FastaIndexedSequence;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class JsonGenomeLoader extends GenomeLoader {

    private static Logger log = Logger.getLogger(JsonGenomeLoader.class);

    private String genomePath;

    public JsonGenomeLoader(String genomePath) {
        this.genomePath = genomePath;
    }

    @Override
    public Genome loadGenome() throws IOException {

        BufferedReader reader = null;

        try {
            reader = ParsingUtils.openBufferedReader(genomePath);
            JsonParser parser = new JsonParser();
            JsonObject json = parser.parse(reader).getAsJsonObject();

            String id = json.get("id").getAsString();
            String name = json.get("name").getAsString();

            String fastaPath;
            String indexPath = null;
            if (json.has("compressedFastaURL")) {
                JsonElement fastaElement = json.has("compressedFastaURL") ?
                        json.get("compressedFastaURL") :
                        json.get("fastaURL");
                fastaPath = fastaElement.getAsString();
                // index path ignored for bgzipped fasta
            } else {
                fastaPath = json.get("fastaURL").getAsString();
                JsonElement indexPathObject = json.get("indexURL");
                indexPath = indexPathObject == null ? null : indexPathObject.getAsString();
            }

            fastaPath = FileUtils.getAbsolutePath(fastaPath, genomePath);
            if (indexPath != null) {
                indexPath = FileUtils.getAbsolutePath(indexPath, genomePath);
            }

            FastaIndexedSequence sequence = fastaPath.endsWith(".gz") ?
                    new FastaBlockCompressedSequence(fastaPath, indexPath) :
                    new FastaIndexedSequence(fastaPath, indexPath);

            JsonElement orderedElement = json.get("ordered");
            boolean ordered = orderedElement != null && orderedElement.getAsBoolean();

            ArrayList<ResourceLocator> tracks = new ArrayList<>();
            ArrayList<ResourceLocator> hiddenTracks = new ArrayList<>();
            JsonArray annotations = json.getAsJsonArray("tracks");
            if (annotations == null) {
                annotations = json.getAsJsonArray("annotations");
            }
            if (annotations != null) {
                annotations.forEach((JsonElement jsonElement) -> {
                    JsonObject obj = jsonElement.getAsJsonObject();
                    String trackPath = obj.get("url").getAsString();
                    JsonElement trackName = obj.get("name");
                    JsonElement trackIndex = obj.get("indexURL");
                    JsonElement indexedElement = obj.get("indexed");
                    JsonElement hiddenElement = obj.get("hidden");
                    boolean hidden = hiddenElement != null && hiddenElement.getAsBoolean();
                    boolean indexed = indexedElement != null && indexedElement.getAsBoolean();

                    String trackIndexPath = null;
                    if (trackPath != null) {
                        trackPath = FileUtils.getAbsolutePath(trackPath, genomePath);
                    }
                    if (trackIndex != null) {
                        trackIndexPath = FileUtils.getAbsolutePath(trackIndex.getAsString(), genomePath);
                    }

                    ResourceLocator res = new ResourceLocator(trackPath);
                    if (trackName != null) res.setName(trackName.getAsString());
                    if (trackIndexPath != null) res.setIndexPath(trackIndexPath);
                    if (indexedElement != null) res.setIndexed(indexed);

                    if (hidden) {
                        if (indexed || trackIndex != null) {
                            log.info("Hidden tracks cannot be indexed.  Ignoring " + trackPath);
                        } else {
                            hiddenTracks.add(res);
                        }
                    } else {
                        tracks.add(res);
                    }
                });
            }

            Genome newGenome = new Genome(id, name, sequence, ordered);
            newGenome.setAnnotationResources(tracks);

            JsonElement ucscIDElement = json.get("ucscID");
            if (ucscIDElement != null) {
                newGenome.setUcscID( ucscIDElement.getAsString());
            }
            JsonElement blatDB = json.get("blatDB");
            if (blatDB != null) {
                newGenome.setUcscID( blatDB.getAsString());
            }
            JsonElement aliasURL = json.get("aliasURL");
            if (aliasURL != null) {
                newGenome.addChrAliases(GenomeLoader.loadChrAliases(aliasURL.getAsString()));
            }
            if (hiddenTracks.size() > 0) {
                addToFeatureDB(hiddenTracks, newGenome);
            }
            return newGenome;
        } finally {
            reader.close();
        }
    }

    public GenomeDescriptor loadDescriptor() throws IOException {
        BufferedReader reader = null;
        try {
            reader = ParsingUtils.openBufferedReader(genomePath);
            JsonParser parser = new JsonParser();
            JsonObject json = parser.parse(reader).getAsJsonObject();
            String id = json.get("id").getAsString();
            String name = json.get("name").getAsString();
            String fastaPath = json.get("fastaURL").getAsString();
            return new GenomeDescriptor(id, name, fastaPath);
        } finally {
            reader.close();
        }
    }


    private void addToFeatureDB(List<ResourceLocator> locators, Genome genome) {
        for (ResourceLocator locator : locators) {
            try {
                FeatureReader featureReader = TribbleFeatureSource.getBasicReader(locator, genome);
                CloseableTribbleIterator<Feature> iter = featureReader.iterator();
                while (iter.hasNext()) {
                    Feature f = iter.next();
                    if (f instanceof NamedFeature) {
                        FeatureDB.addFeature((NamedFeature) f, genome);
                    }
                }
            } catch (IOException e) {
                log.error("Error loading " + locator.getPath());
            }
        }
    }


    public static class GenomeDescriptor {
        String id;
        String name;
        String fastaURL;

        public GenomeDescriptor(String id, String name, String fastaURL) {
            this.id = id;
            this.name = name;
            this.fastaURL = fastaURL;
        }

        public String getId() {
            return id;
        }

        public String getName() {
            return name;
        }

        public String getFastaURL() {
            return fastaURL;
        }
    }

}
