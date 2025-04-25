package org.broad.igv.ucsc.hub;

import org.broad.igv.Globals;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class HubParser {

    private static Logger log = LogManager.getLogger(HubParser.class);

    private static Set<String> urlProperties = new HashSet<>(Arrays.asList("descriptionUrl", "desriptionUrl",
            "twoBitPath", "blat", "chromAliasBb", "twoBitBptURL", "twoBitBptUrl", "htmlPath", "bigDataUrl",
            "genomesFile", "trackDb", "groups", "include", "html", "searchTrix", "groups",
            "chromSizes"));

    private static final Map<String, Hub> hubCache = new HashMap<>();

    public static Hub loadAssemblyHub(String url) throws IOException {
        return loadHub(url, null);
    }

    public static Hub loadHub(String url, String genomeId) throws IOException {

        // Check if the Hub is already cached
        String key = url + (genomeId == null ? "" : genomeId);
        if (hubCache.containsKey(key)) {

            return hubCache.get(key);
        }

        log.info("Loading Hub: " + url);

        // Load stanzas from the hub.txt file, which might be all stanzas if 'useOneFile' is on
        List<Stanza> stanzas = HubParser.loadStanzas(url);

        // Validation checks
        if (stanzas.size() < 1) {
            throw new RuntimeException("Expected at least 1 stanza");
        }
        if (!"hub".equals(stanzas.get(0).type)) {
            throw new RuntimeException("Unexpected hub.txt file -- does the first line start with 'hub'?");
        }
        Stanza hubStanza = stanzas.get(0);
        List<Stanza> trackStanzas = null;
        List<Stanza> groupStanzas = null;
        List<Stanza> genomeStanzas = null;

        if ("on".equals(stanzas.get(0).getProperty("useOneFile"))) {

            if (!"genome".equals(stanzas.get(1).type)) {
                throw new RuntimeException("Unexpected hub file -- expected 'genome' stanza but found " + stanzas.get(1).type);
            }
            Stanza genomeStanza = stanzas.get(1);
            genomeStanzas = Arrays.asList(genomeStanza);

            final boolean isAssemblyHub = genomeStanza.hasProperty("twoBitPath");
            if (!isAssemblyHub && genomeId != null) {
                // Not an assembly hub, validate genome
                if (!genomeStanza.getProperty("genome").equals(genomeId)) {
                    throw new RuntimeException("Hub file does not contain tracks for genome " + genomeId);
                }
            }

            trackStanzas = new ArrayList<>(stanzas.subList(2, stanzas.size()));

        } else {
            if (!hubStanza.hasProperty("genomesFile")) {
                throw new RuntimeException("hub.txt must specify 'genomesFile'");
            }
            genomeStanzas = HubParser.loadStanzas(hubStanza.getProperty("genomesFile"));
        }

        // Load group files for all genomes, if any.
        groupStanzas = new ArrayList<>();
        Set<String> uniqGroupURLs = genomeStanzas.stream().map(s -> s.getProperty("groups")).filter(Objects::nonNull).collect(Collectors.toSet());
        for (String groupTxtURL : uniqGroupURLs) {
            groupStanzas.addAll(HubParser.loadStanzas(groupTxtURL));
        }

        Hub hub = new Hub(url, hubStanza, genomeStanzas, trackStanzas, groupStanzas);

        // Cache the loaded Hub
        hubCache.put(key, hub);

        return hub;
    }

    /**
     * Load track hubs in parallel, but set a timeout to prevent a non-responsive hub server from preventing genome load
     * @param ucscId
     * @param hubUrls
     * @return
     */
    public static List<Hub> loadHubs(String ucscId, List<String> hubUrls) {

        List<Hub> trackHubs = new ArrayList<>();

        List<CompletableFuture<Object>> futures = IntStream.range(0, hubUrls.size())
                .parallel()
                .mapToObj(i -> CompletableFuture.supplyAsync(() -> {
                    try {
                        int order = i + 1;
                        final Hub hub = HubParser.loadHub(hubUrls.get(i), ucscId);
                        hub.setOrder(order);
                        trackHubs.add(hub);
                    } catch (Exception e) {
                        log.error("Error loading hub " + hubUrls.get(i), e);
                    }
                    return null;
                }))
                .collect(Collectors.toList());
        CompletableFuture<Void> combinedFuture = CompletableFuture.allOf(futures.toArray(new CompletableFuture[0]));
        try {
            combinedFuture.get(10, TimeUnit.SECONDS);
        } catch (Exception e) {
            log.error("Error loading hubs", e);
        }
        Collections.sort(trackHubs, (h1, h2) -> h1.getOrder() - h2.getOrder());
        return trackHubs;
    }

    static List<Stanza> loadStanzas(String url) throws IOException {

        int idx = url.lastIndexOf("/");
        String baseURL = url.substring(0, idx + 1);
        String host = getHost(url);

        List<Stanza> nodes = new ArrayList<>();
        Stanza currentNode = null;
        boolean startNewNode = true;
        int order = 0;
        try (BufferedReader br = ParsingUtils.openBufferedReader(url)) {
            String line;
            while ((line = br.readLine()) != null) {

                if (line.startsWith("#")) {
                    continue;
                }

                while (line.endsWith("\\")) {
                    String continuation = br.readLine();
                    if (continuation == null) {
                        break;
                    } else {
                        line = line.substring(0, line.length() - 1) + " " + continuation.trim();
                    }
                }

                if (line.startsWith("include")) {
                    String relativeURL = line.substring(8).trim();
                    String includeURL = getDataURL(relativeURL, host, baseURL);
                    List<Stanza> includeStanzas = HubParser.loadStanzas(includeURL);
                    nodes.addAll(includeStanzas);
                }

                int indent = indentLevel(line);
                int i = line.indexOf(" ", indent);
                if (i < 0) {
                    // Break - start a new node
                    startNewNode = true;
                } else {
                    String key = line.substring(indent, i);
                    if (key.startsWith("#")) continue;

                    String value = line.substring(i + 1).trim();

                    if ("type".equals(key)) {
                        // The "type" property contains format and sometimes other information.  For example, data range
                        // on a bigwig "type bigWig 0 .5"
                        String[] tokens = Globals.whitespacePattern.split(value);
                        value = tokens[0];
                        if ("bigWig".equals(value) && tokens.length == 3) {
                            // This is a bigWig with a range
                            String min = tokens[1];
                            String max = tokens[2];
                            if (currentNode != null) {
                                currentNode.properties.put("min", min);
                                currentNode.properties.put("max", max);
                            }
                        }

                    } else if (!("shortLabel".equals(key) || "longLabel".equals(key) || "metadata".equals(key) || "label".equals(key))) {
                        String[] tokens = Globals.whitespacePattern.split(value);
                        value = tokens[0];
                    }

                    if (urlProperties.contains(key) || value.endsWith("URL") || value.endsWith("Url")) {
                        value = getDataURL(value, host, baseURL);
                    }

                    if (startNewNode) {
                        // Start a new node -- indent is currently ignored as igv.js does not support sub-tracks,
                        // so track stanzas are flattened
                        Stanza newNode = new Stanza(key, value);
                        nodes.add(newNode);
                        currentNode = newNode;
                        startNewNode = false;
                    }
                    currentNode.properties.put(key, value);
                }
            }
        }

        return resolveParents(nodes);
    }

    private static int indentLevel(String str) {
        int level;
        for (level = 0; level < str.length(); level++) {
            char c = str.charAt(level);
            if (c != ' ' && c != '\t') break;
        }
        return level;
    }

    private static List<Stanza> resolveParents(List<Stanza> nodes) {
        Map<String, Stanza> nodeMap = new HashMap<>();
        for (Stanza n : nodes) {
            nodeMap.put(n.name, n);
        }
        for (Stanza n : nodes) {
            if (n.properties.containsKey("parent")) {
                String parentName = firstWord(n.properties.get("parent"));
                n.parent = nodeMap.get(parentName);
                n.properties.put("parent", parentName);
            }
        }
        return nodes;
    }

    static String firstWord(String str) {
        return Globals.whitespacePattern.split(str)[0];
    }

    private static String getDataURL(String url, String host, String baseURL) {
        return url.startsWith("http://") || url.startsWith("https://") ? url :
                url.startsWith("/") ? host + url : baseURL + url;
    }

    private static String getHost(String url) {
        String host;
        if (url.startsWith("https://") || url.startsWith("http://")) {
            try {
                URL tmp = new URL(url);
                host = tmp.getProtocol() + "://" + tmp.getHost();
            } catch (MalformedURLException e) {
                // This should never happen
                log.error("Error parsing base URL host", e);
                throw new RuntimeException(e);
            }
        } else {
            // Local file, no host
            host = "";
        }
        return host;
    }


}
