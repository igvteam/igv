package org.broad.igv.feature.genome.load;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureReader;
import org.broad.igv.logging.*;
import org.broad.igv.Globals;
import org.broad.igv.feature.CytoBandFileParser;
import org.broad.igv.feature.FeatureDB;
import org.broad.igv.feature.NamedFeature;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.feature.genome.fasta.FastaBlockCompressedSequence;
import org.broad.igv.feature.genome.fasta.FastaIndexedSequence;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackProperties;
import org.broad.igv.track.TribbleFeatureSource;
import org.broad.igv.ui.color.ColorUtilities;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.liftover.Liftover;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;

public class JsonGenomeLoader extends GenomeLoader {

    private static Logger log = LogManager.getLogger(JsonGenomeLoader.class);

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
            fastaPath = json.get("fastaURL").getAsString();

            JsonElement indexPathObject = json.get("indexURL");
            String indexPath = indexPathObject == null ? null : indexPathObject.getAsString();

            JsonElement gziObject = json.get("gziIndexURL");
            if (gziObject == null) {
                gziObject = json.get("compressedIndexURL");
            }
            String gziIndexPath = gziObject == null ? null : gziObject.getAsString();

            fastaPath = FileUtils.getAbsolutePath(fastaPath, genomePath);
            if (indexPath != null) {
                indexPath = FileUtils.getAbsolutePath(indexPath, genomePath);
            }
            if (gziIndexPath != null) {
                gziIndexPath = FileUtils.getAbsolutePath(gziIndexPath, genomePath);
            }

            FastaIndexedSequence sequence = fastaPath.endsWith(".gz") ?
                    new FastaBlockCompressedSequence(fastaPath, gziIndexPath, indexPath) :
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
                    if (trackPath != null) {
                        trackPath = FileUtils.getAbsolutePath(trackPath, genomePath);
                    }

                    ResourceLocator res = new ResourceLocator(trackPath);

                    JsonElement trackName = obj.get("name");
                    if (trackName != null) {
                        res.setName(trackName.getAsString());
                    }

                    JsonElement trackIndex = obj.get("indexURL");
                    if (trackIndex != null) {
                        res.setIndexPath(FileUtils.getAbsolutePath(trackIndex.getAsString(), genomePath));
                    }

                    JsonElement format = obj.get("format");
                    if (format != null) {
                        res.setFormat(format.getAsString().toLowerCase());
                    }

//                    JsonElement colorBy = obj.get("colorBy");
//                    if (altColor != null) {
//                        try {
//                            properties.colorBy(ColorUtilities.stringToColor(colorBy.toString()));
//                        } catch (Exception e) {
//                            log.error("Error parsing color string: " + altColor.toString(), e);
//                        }
//                    }

                    JsonElement infoURL = obj.get("infoURL");
                    if (infoURL != null) {
                        res.setFeatureInfoURL(infoURL.getAsString());
                    }


                    // Track properties
                    TrackProperties properties = new TrackProperties();
                    JsonElement color = obj.get("color");
                    if (color != null) {
                        try {
                            properties.setColor(ColorUtilities.stringToColor(color.toString()));
                        } catch (Exception e) {
                            log.error("Error parsing color string: " + color, e);
                        }
                    }
                    JsonElement altColor = obj.get("altColor");
                    if (altColor != null) {
                        try {
                            properties.setAltColor(ColorUtilities.stringToColor(altColor.toString()));
                        } catch (Exception e) {
                            log.error("Error parsing color string: " + altColor, e);
                        }
                    }
                    JsonElement displayMode = obj.get("displayMode");
                    if (displayMode != null) {
                        try {
                            Track.DisplayMode dp = Track.DisplayMode.valueOf(stripQuotes(displayMode.toString()));
                            properties.setDisplayMode(dp);
                        } catch (Exception e) {
                            log.error("Error parsing displayMode " + displayMode, e);
                        }
                    }
                    JsonElement vizwindow = obj.get("visibilityWindow");
                    if (vizwindow != null) {
                        properties.setFeatureVisibilityWindow(obj.get("visibilityWindow").getAsInt());
                    } else {
                        // If not explicitly set, assume whole chromosome viz window for annotations
                        properties.setFeatureVisibilityWindow(-1);
                    }
                    res.setTrackProperties(properties);

                    JsonElement indexedElement = obj.get("indexed");
                    JsonElement hiddenElement = obj.get("hidden");
                    boolean hidden = hiddenElement != null && hiddenElement.getAsBoolean();
                    boolean indexed = indexedElement != null && indexedElement.getAsBoolean();
                    if (indexedElement != null) {
                        res.setIndexed(indexed);
                    }
                    if (hidden) {
                        if (indexed || trackIndex != null) {
                            log.warn("Hidden tracks cannot be indexed.  Ignoring " + trackPath);
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

            JsonElement cyobandElement = json.get("cytobandURL");
            if (cyobandElement != null) {
                String cytobandPath = FileUtils.getAbsolutePath(cyobandElement.getAsString(), genomePath);
                BufferedReader br = null;
                try {
                    br = ParsingUtils.openBufferedReader(cytobandPath);
                    newGenome.setCytobands(CytoBandFileParser.loadData(br));
                } finally {
                    try {
                        br.close();
                    } catch (IOException e) {
                        // ignore
                    }
                }
            }

            JsonElement blatDB = json.get("blatDB");
            if (blatDB != null) {
                newGenome.setUcscID(blatDB.getAsString());
            }
            JsonElement ucscIDElement = json.get("ucscID");
            if (ucscIDElement != null) {
                newGenome.setUcscID(ucscIDElement.getAsString());
                if(blatDB == null) {
                    newGenome.setBlatDB(ucscIDElement.getAsString());
                }
            }

            JsonElement aliasURL = json.get("aliasURL");
            if (aliasURL != null) {
                String aliasPath = FileUtils.getAbsolutePath(aliasURL.getAsString(), genomePath);
                newGenome.addChrAliases(GenomeLoader.loadChrAliases(aliasPath));
            }
            if (hiddenTracks.size() > 0) {
                addToFeatureDB(hiddenTracks, newGenome);
            }
            JsonElement wholeGenomeView = json.get("wholeGenomeView");
            if (wholeGenomeView != null) {
                newGenome.setShowWholeGenomeView(wholeGenomeView.getAsBoolean());
            }
            JsonElement chromosomeOrder = json.get("chromosomeOrder");
            if (chromosomeOrder != null) {
                List<String> chrs;
                if (chromosomeOrder.isJsonArray()) {
                    JsonArray a = chromosomeOrder.getAsJsonArray();
                    chrs = new ArrayList<>();
                    for (JsonElement e : a) {
                        chrs.add(e.getAsString());
                    }
                } else {
                    // Assume comma delimited stream
                    String[] c = Globals.commaPattern.split(chromosomeOrder.getAsString());
                    chrs = new ArrayList<>(c.length);
                    for (String t : c) chrs.add(t.trim());
                }
                newGenome.setLongChromosomeNames(chrs);
            }

            // Load liftover "chain" files.  This enables navigating by coordinates of another genome.
            // Not a common option.

            JsonElement chains = json.get("chains");
            if (chains != null) {
                Map<String, Liftover> liftoverMap = new HashMap<>();
                JsonObject chainsObj = chains.getAsJsonObject();
                for (Map.Entry<String, JsonElement> entry : chainsObj.entrySet()) {
                    String chainsPath = FileUtils.getAbsolutePath(entry.getValue().getAsString(), genomePath);
                    liftoverMap.put(entry.getKey(), Liftover.load(chainsPath));

                }
                newGenome.setLiftoverMap(liftoverMap);
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

    private String stripQuotes(String str) {
        if(str.startsWith("\"")) {
            return str.substring(1, str.length() - 1);  // Assume also ends with
        }
        else {
            return str;
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
