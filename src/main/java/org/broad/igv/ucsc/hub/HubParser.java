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

public class HubParser {

    private static Logger log = LogManager.getLogger(HubParser.class);

    private static Set<String> urlProperties = new HashSet<>(Arrays.asList("descriptionUrl", "desriptionUrl",
            "twoBitPath", "blat", "chromAliasBb", "twoBitBptURL", "twoBitBptUrl", "htmlPath", "bigDataUrl",
            "genomesFile", "trackDb", "groups", "include", "html", "searchTrix"));

    public static Hub loadAssemblyHub(String url) throws IOException {
        return loadHub(url, null);
    }

    public static Hub loadHub(String url, String genomeId) throws IOException {

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
        Stanza genomeStanza = null;
        List<Stanza> trackStanzas = null;
        List<Stanza> groupStanzas = null;
        String trackDbURL = null;

        if ("on".equals(stanzas.get(0).getProperty("useOneFile"))) {

            if (!"genome".equals(stanzas.get(1).type)) {
                throw new RuntimeException("Unexpected hub file -- expected 'genome' stanza but found " + stanzas.get(1).type);
            }
            genomeStanza = stanzas.get(1);

            if(!genomeStanza.hasProperty("twoBitPath")) {
                // Not an assembly hub, validate genome
                if(!genomeStanza.getProperty("genome").equals(genomeId)) {
                    throw new RuntimeException("Hub file does not contain tracks for genome " + genomeId);
                }
            }

            trackStanzas = new ArrayList<>(stanzas.subList(2, stanzas.size()));

        } else {

            if (!hubStanza.hasProperty("genomesFile")) {
                throw new RuntimeException("hub.txt must specify 'genomesFile'");
            }

            // Load genome stanza(s)
            List<Stanza> genomeStanzaList = HubParser.loadStanzas(hubStanza.getProperty("genomesFile"));

            // Find stanza for requested genome, or if no requested genome the first.

            for (Stanza s : genomeStanzaList) {
                if (genomeId == null || genomeId.equals(s.getProperty("genome"))) {
                    stanzas.add(s);
                    genomeStanza = s;
                    trackDbURL = genomeStanza.getProperty("trackDb");
                    break;
                }
            }
            if(trackDbURL == null) {
                throw new RuntimeException("Hub file does not contain tracks for genome " + genomeId);
            }
        }

        // Assembly hub validation
        if (genomeId == null ) {
            if (!genomeStanza.hasProperty("twoBitPath")) {
                throw new RuntimeException("Assembly hubs must specify 'twoBitPath'");
            }
        }

        if (genomeStanza.hasProperty("groups")) {
            if (genomeStanza.hasProperty("groups")) {
                String groupsTxtURL = genomeStanza.getProperty("groups");
                groupStanzas = HubParser.loadStanzas(groupsTxtURL);
            }
        }

        return new Hub(url, trackDbURL, hubStanza, genomeStanza, trackStanzas, groupStanzas);
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

                if(line.startsWith("include")) {
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

                    if (!("shortLabel".equals(key) || "longLabel".equals(key) || "metadata".equals(key))) {
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

}
