package org.igv.ucsc.hub;

import org.igv.Globals;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.util.FileUtils;
import org.igv.util.ParsingUtils;

import java.io.BufferedReader;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class HubParser {

    private static Logger log = LogManager.getLogger(HubParser.class);

    private static Set<String> urlProperties = new HashSet<>(Arrays.asList("descriptionUrl", "desriptionUrl",
            "twoBitPath", "blat", "chromAliasBb", "twoBitBptURL", "twoBitBptUrl", "htmlPath", "bigDataUrl",
            "genomesFile", "trackDb", "groups", "include", "html", "searchTrix", "linkDataUrl", "groups",
            "chromSizes"));

    /**
     * Total time to wait for the hubs of a genome to load.  Hubs that have not loaded within this budget are skipped
     * so that a slow or non-responsive hub server cannot block genome loading.  This is a budget for the genome's
     * hubs as a whole, not per hub, so it can be generous without the wait scaling with the number of hubs.
     */
    private static final long HUB_LOAD_BUDGET_SECONDS = 20;

    // Accessed from the hub loading threads started in loadHubs(), as well as the calling thread.
    private static final Map<String, Hub> hubCache = new ConcurrentHashMap<>();

    public static Hub loadHub(String url) throws IOException {

        // Check if the Hub is already cached
        String key = url;
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

        if ("on".equals(hubStanza.getProperty("useOneFile"))) {

            // This is a "onefile" hub, all stanzas are in the same file
            if (!"genome".equals(stanzas.get(1).type)) {
                throw new RuntimeException("Unexpected hub file -- expected 'genome' stanza but found " + stanzas.get(1).type);
            }
            Stanza genomeStanza = stanzas.get(1);
            genomeStanzas = Arrays.asList(genomeStanza);
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
     * The result of loading a list of hubs -- those that loaded, and the urls of those that did not.  The failures
     * are reported so that the caller can make them visible;  a hub that quietly disappears from the menu looks
     * like a bug in IGV.
     */
    public record HubLoadResult(List<Hub> hubs, List<String> failedUrls) {
    }

    /**
     * Load track hubs in parallel, subject to an overall time budget (HUB_LOAD_BUDGET_SECONDS, shared by all of the
     * genome's hubs) so that a slow or non-responsive hub server cannot block genome loading.  Hubs that fail to
     * load, or that do not load within the budget, are reported as failures -- the url and reason are logged in
     * both cases.
     * <p>
     * A hub that exceeds the budget is not abandoned, it is left to finish on its daemon thread and enters the hub
     * cache, so it will be available the next time the genome is loaded.
     *
     * @param hubUrls
     * @return the hubs that loaded, in the order given by hubUrls, and the urls of those that did not
     */
    public static HubLoadResult loadHubs(List<String> hubUrls) {

        if (hubUrls == null || hubUrls.isEmpty()) {
            return new HubLoadResult(Collections.emptyList(), Collections.emptyList());
        }

        final long deadline = System.nanoTime() + TimeUnit.SECONDS.toNanos(HUB_LOAD_BUDGET_SECONDS);

        ExecutorService executor = Executors.newFixedThreadPool(hubUrls.size(), runnable -> {
            Thread thread = new Thread(runnable, "hub-loader");
            thread.setDaemon(true);
            return thread;
        });

        try {
            List<Future<Hub>> futures = hubUrls.stream()
                    .map(url -> executor.submit(() -> {
                        Hub hub = loadHub(url);
                        if (System.nanoTime() > deadline) {
                            log.warn("Track hub loaded after the " + HUB_LOAD_BUDGET_SECONDS + " second budget: " + url +
                                    ".  It will be available the next time this genome is loaded.");
                        }
                        return hub;
                    }))
                    .collect(Collectors.toList());

            List<Hub> hubs = new ArrayList<>(futures.size());
            List<String> failedUrls = new ArrayList<>();
            for (int i = 0; i < futures.size(); i++) {
                String url = hubUrls.get(i);
                try {
                    // The timeout is a budget shared by all hubs, so the total wait is bounded regardless of their number
                    Hub hub = futures.get(i).get(Math.max(0, deadline - System.nanoTime()), TimeUnit.NANOSECONDS);
                    hub.setOrder(i + 1);
                    hubs.add(hub);
                } catch (TimeoutException e) {
                    log.error("Track hub not loaded within the " + HUB_LOAD_BUDGET_SECONDS +
                            " second budget for this genome's hubs: " + url);
                    failedUrls.add(url);
                } catch (ExecutionException e) {
                    log.error("Error loading track hub: " + url, e.getCause());
                    failedUrls.add(url);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    log.error("Interrupted loading track hubs: " + hubUrls.subList(i, hubUrls.size()));
                    failedUrls.addAll(hubUrls.subList(i, hubUrls.size()));
                    break;
                }
            }
            return new HubLoadResult(hubs, failedUrls);

        } finally {
            // Does not interrupt hubs still loading;  they complete on their daemon threads and populate the cache.
            executor.shutdown();
        }
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
                    List<Stanza> includeStanzas = loadStanzas(includeURL);
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
                        currentNode = new Stanza(key, value);
                        nodes.add(currentNode);
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
                final Stanza parent = nodeMap.get(parentName);
                n.setParent(parent);
                n.properties.put("parent", parentName);
            }
        }
        return nodes;
    }

    static String firstWord(String str) {
        return Globals.whitespacePattern.split(str)[0];
    }

    private static String getDataURL(String url, String host, String baseURL) {
        return FileUtils.isRemote(url) ? url :
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
