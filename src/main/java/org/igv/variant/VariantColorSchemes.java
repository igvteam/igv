package org.igv.variant;

import org.igv.DirectoryManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.AbstractColorScale;
import org.igv.renderer.ContinuousColorScale;
import org.igv.util.HttpUtils;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringReader;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * The color schemes available for coloring variants by a VCF INFO attribute.
 * <p>
 * IGV keeps one scheme per attribute, in a file named for it: "CLNSIG.txt" in the IGV directory
 * ("variantcolors").  The directory is the store, scanned at startup -- the attribute is the file name, so there
 * is never more than one scheme for an attribute, and a user's scheme for CLNSIG replaces IGV's entirely.
 * Schemes are files rather than individual preferences because the set of INFO keys, and of values for each, is
 * unbounded, so there is nothing to enumerate in a preferences grid.  Schemes shipped with IGV are read from the
 * classpath and are never copied, so upgrades can revise them.
 * <p>
 * An imported file may cover any number of attributes.  It is read once and written into the store as one file
 * per attribute, replacing whatever the store held for each; the original is not kept or referred to again.
 */
public class VariantColorSchemes {

    private static Logger log = LogManager.getLogger(VariantColorSchemes.class);

    static final String SCHEME_DIRECTORY = "variantcolors";

    /**
     * Schemes shipped with IGV, one per INFO attribute.
     */
    private static final String[] BUILTIN_RESOURCES = {
            "resources/clnsig.txt",
            "resources/svtype.txt",
            "resources/vt.txt"
    };

    /**
     * INFO key (upper case) -> scheme, sorted by key.
     */
    private static Map<String, VariantColorScheme> userSchemes;
    private static Map<String, VariantColorScheme> builtinSchemes;

    /**
     * The scheme in effect for each attribute.  Colors are resolved through this for every variant drawn, so it
     * is built once and discarded only when the schemes change.
     */
    private static Map<String, VariantColorScheme> allSchemes;

    /**
     * The IGV directory the cached user schemes were read from.  The IGV directory can move (Preferences >
     * Advanced), which leaves cached schemes holding paths into a directory that no longer exists -- saving them
     * would fail and removals would come back on restart.  Compared by identity with DirectoryManager's own
     * File, so checking it allocates nothing.
     */
    private static File cacheIgvDirectory;

    private VariantColorSchemes() {
    }

    /**
     * @return the scheme in effect for each attribute, sorted by attribute: the user's where there is one,
     * otherwise the one shipped with IGV.
     */
    public static synchronized List<VariantColorScheme> getSchemes() {
        return new ArrayList<>(getSchemeMap().values());
    }

    private static synchronized Map<String, VariantColorScheme> getSchemeMap() {
        loadUserSchemesIfNeeded();
        if (allSchemes == null) {
            Map<String, VariantColorScheme> schemes = new TreeMap<>(getBuiltinSchemeMap());
            schemes.putAll(userSchemes);
            allSchemes = Collections.unmodifiableMap(schemes);
        }
        return allSchemes;
    }

    public static synchronized List<VariantColorScheme> getUserSchemes() {
        loadUserSchemesIfNeeded();
        return new ArrayList<>(userSchemes.values());
    }

    /**
     * Load the user schemes if not yet loaded, or if the IGV directory has moved.  Called on every color lookup,
     * so it allocates nothing on the usual path.
     */
    private static synchronized void loadUserSchemesIfNeeded() {
        File igvDirectory = DirectoryManager.getIgvDirectory();
        if (userSchemes == null || igvDirectory != cacheIgvDirectory) {
            cacheIgvDirectory = igvDirectory;
            userSchemes = loadUserSchemes();
            allSchemes = null;
        }
    }

    public static synchronized List<VariantColorScheme> getBuiltinSchemes() {
        return new ArrayList<>(getBuiltinSchemeMap().values());
    }

    private static synchronized Map<String, VariantColorScheme> getBuiltinSchemeMap() {
        if (builtinSchemes == null) {
            builtinSchemes = loadBuiltinSchemes();
        }
        return builtinSchemes;
    }

    /**
     * @return the scheme in effect for an attribute, or null if none covers it.
     */
    public static VariantColorScheme getScheme(String infoKey) {
        return infoKey == null ? null : getSchemeMap().get(infoKey.toUpperCase());
    }

    /**
     * Return the color assigned to an INFO attribute value, or null if no scheme covers it, in which case the
     * caller assigns a color from a palette.
     */
    public static Color getColor(String infoKey, String value) {
        VariantColorScheme scheme = getScheme(infoKey);
        return scheme == null ? null : scheme.getColor(infoKey, value);
    }

    /**
     * Return the scale for a numeric INFO attribute, or null if there is none.  The scheme covering the attribute
     * decides.  Otherwise a well-known allele frequency field (GMAF, gnomAD_AF, ...) is colored by rarity -- see
     * {@link AlleleFrequencyColors} -- and other numeric attributes have no scale, as there is no sensible default.
     */
    public static AbstractColorScale getScale(String infoKey) {
        VariantColorScheme scheme = getScheme(infoKey);
        if (scheme != null) {
            return scheme.getScale(infoKey);
        }
        return AlleleFrequencyColors.isFrequencyField(infoKey) ? AlleleFrequencyColors.getScale() : null;
    }

    /**
     * @return the scheme providing the color scale for an INFO attribute, or null if none does.
     */
    public static VariantColorScheme getSchemeForScale(String infoKey) {
        VariantColorScheme scheme = getScheme(infoKey);
        return scheme != null && scheme.getScale(infoKey) != null ? scheme : null;
    }

    /**
     * @return the colors the scheme in effect assigns to values of an INFO attribute.  Used to keep colors
     * assigned from the palette distinguishable from them.
     */
    public static Collection<Color> getColors(String infoKey) {
        VariantColorScheme scheme = getScheme(infoKey);
        return scheme == null ? Collections.emptyList() : scheme.getColors(infoKey).values();
    }

    /**
     * @return every INFO key covered by some scheme.  Used to offer numeric attributes that would otherwise not
     * be colorable.
     */
    public static Set<String> getKeys() {
        return getSchemeMap().keySet();
    }

    /**
     * Import color schemes from a file.  The file is read once and written into the store as one file per
     * attribute it covers, so the schemes survive the original being moved or deleted and are there on the next
     * startup.  The store's colors for each of those attributes are replaced.
     *
     * @return the imported schemes, one per attribute
     */
    public static synchronized List<VariantColorScheme> importFile(File file) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8))) {
            return importSchemes(reader, file.getAbsolutePath());
        }
    }

    /**
     * Import color schemes from a URL, as {@link #importFile} does from a file.
     *
     * @return the imported schemes, one per attribute
     */
    public static List<VariantColorScheme> importUrl(String url) throws IOException {
        // Not synchronized: the fetch can block for the full connect and read timeouts, and the event thread
        // takes this class's lock on every repaint through the scheme accessors.  Fetch first, and hold the lock
        // only to parse and write.
        String contents = HttpUtils.getInstance().getContentsAsString(new URL(url));
        return importContents(contents, url);
    }

    private static synchronized List<VariantColorScheme> importContents(String contents, String source)
            throws IOException {
        return importSchemes(new BufferedReader(new StringReader(contents)), source);
    }

    private static List<VariantColorScheme> importSchemes(BufferedReader reader, String source) throws IOException {

        VariantColorScheme imported = VariantColorScheme.parse(reader, null);
        if (imported.getKeys().isEmpty()) {
            throw new IOException("No INFO attribute colors found in " + source);
        }

        List<VariantColorScheme> schemes = new ArrayList<>();
        for (String key : imported.getKeys()) {
            schemes.add(save(imported.forKey(key)));
        }
        return schemes;
    }

    /**
     * Write the given value -> color assignments as the scheme for an attribute, so they apply to every VCF with
     * this attribute rather than only the track they were chosen on.  The store's colors for the attribute are
     * replaced.
     *
     * @param infoKey  the INFO attribute the colors are for
     * @param colors   value -> color
     * @return the saved scheme
     */
    public static synchronized VariantColorScheme saveScheme(String infoKey, Map<String, Color> colors)
            throws IOException {
        return saveScheme(infoKey, colors, false);
    }

    /**
     * @param categorical true to record that this numeric attribute holds categories rather than quantities.  The
     *                    declaration is what persists the answer -- the colors only cover the values seen so far.
     */
    public static synchronized VariantColorScheme saveScheme(String infoKey, Map<String, Color> colors,
                                                             boolean categorical) throws IOException {

        VariantColorScheme scheme = new VariantColorScheme(infoKey);
        if (categorical) {
            scheme.setCategorical(infoKey);
        }
        for (Map.Entry<String, Color> entry : colors.entrySet()) {
            scheme.setColor(infoKey, entry.getKey(), entry.getValue());
        }
        return save(scheme);
    }

    /**
     * Write a scheme into the store, replacing the store's scheme for its attribute, and make it current.  A
     * scheme shipped with IGV cannot be written, so saving one writes a copy the user owns, which then replaces
     * it.
     *
     * @return the saved scheme -- for a scheme with no file yet this is a copy, so use it rather than the argument
     */
    public static synchronized VariantColorScheme save(VariantColorScheme scheme) throws IOException {

        Set<String> keys = scheme.getKeys();
        if (keys.size() != 1) {
            throw new IOException("A color scheme file covers one INFO attribute, not " + keys);
        }
        String key = keys.iterator().next();

        if (scheme.getFile() == null) {
            // New, or shipped with IGV.  Work on a copy: giving a cached built-in a file would make it user owned,
            // and removable, while it is still listed among the built-ins -- so it would appear twice.
            scheme = scheme.copy();
        }
        scheme.setFile(new File(createSchemeDirectory(), key + ".txt"));

        try (PrintWriter writer = new PrintWriter(scheme.getFile(), StandardCharsets.UTF_8)) {
            scheme.write(writer);
        }

        loadUserSchemesIfNeeded();
        userSchemes.put(key, scheme);
        allSchemes = null;
        return scheme;
    }

    /**
     * Write a color scale for a numeric INFO attribute as its scheme.  The scale is serialized in the form IGV
     * uses for heat map scales elsewhere, so there is one representation of a color scale, not two.
     *
     * @return the saved scheme
     */
    public static synchronized VariantColorScheme saveScale(String infoKey, ContinuousColorScale scale)
            throws IOException {

        if (scale == null || !VariantColorScheme.isUsableRange(scale.getMinimum(), scale.getMaximum())) {
            throw new IOException("Color scale range must be finite and increasing");
        }
        VariantColorScheme scheme = new VariantColorScheme(infoKey);
        scheme.setScale(infoKey, scale);
        return save(scheme);
    }

    /**
     * Delete a user's scheme.  The scheme shipped with IGV for the attribute, if there is one, applies again.
     * Schemes shipped with IGV cannot be removed.
     *
     * @return true if the scheme was removed
     */
    public static synchronized boolean remove(VariantColorScheme scheme) {

        if (scheme == null || scheme.isBuiltIn()) {
            return false;
        }
        if (scheme.getFile().exists() && !scheme.getFile().delete()) {
            log.error("Could not delete variant color scheme: " + scheme.getFile().getAbsolutePath());
            return false;
        }
        loadUserSchemesIfNeeded();
        userSchemes.values().remove(scheme);
        allSchemes = null;
        return true;
    }

    /**
     * The directory holding the user's schemes.  It is the store -- every file in it is loaded at startup.  It
     * is not created here, only when a scheme is actually written.
     */
    public static File getSchemeDirectory() {
        return new File(DirectoryManager.getIgvDirectory(), SCHEME_DIRECTORY);
    }

    private static File createSchemeDirectory() throws IOException {
        File directory = getSchemeDirectory();
        if (!directory.exists() && !directory.mkdirs()) {
            throw new IOException("Could not create directory " + directory.getAbsolutePath());
        }
        return directory;
    }

    /**
     * Load the store.  The file name is the attribute, so a file is read for that attribute only; rows for any
     * other are dropped with a warning, as is a file that says nothing about its attribute.  Two files can name
     * the same attribute on a case sensitive file system ("clnsig.txt" beside the "CLNSIG.txt" IGV writes); the
     * first by name wins, the same one after every restart, and the other is logged.
     */
    private static Map<String, VariantColorScheme> loadUserSchemes() {

        Map<String, VariantColorScheme> schemes = new TreeMap<>();
        File directory = new File(DirectoryManager.getIgvDirectory(), SCHEME_DIRECTORY);
        File[] files = directory.isDirectory() ? directory.listFiles() : null;
        if (files == null) {
            return schemes;
        }
        Arrays.sort(files, Comparator.comparing(File::getName, String.CASE_INSENSITIVE_ORDER));

        for (File file : files) {
            if (file.isDirectory() || file.isHidden()) {
                continue;
            }
            String key = stripExtension(file.getName()).toUpperCase();
            if (schemes.containsKey(key)) {
                log.warn("Variant color scheme file " + file.getAbsolutePath() + " ignored, " + key
                        + " is already read from " + schemes.get(key).getFile().getName());
                continue;
            }
            try (BufferedReader reader = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8))) {
                VariantColorScheme scheme = VariantColorScheme.parse(reader, key);
                if (!scheme.covers(key)) {
                    log.warn("Variant color scheme file " + file.getAbsolutePath() + " is named for " + key
                            + " but colors " + scheme.getKeys() + ", skipping it.  Use Import to add it.");
                    continue;
                }
                if (scheme.getKeys().size() > 1) {
                    log.warn("Variant color scheme file " + file.getAbsolutePath() + " colors " + scheme.getKeys()
                            + ", reading only " + key);
                    scheme = scheme.forKey(key);
                }
                scheme.setFile(file);
                schemes.put(key, scheme);
            } catch (IOException e) {
                log.error("Error loading variant color scheme: " + file.getAbsolutePath(), e);
            }
        }
        return schemes;
    }

    private static Map<String, VariantColorScheme> loadBuiltinSchemes() {

        Map<String, VariantColorScheme> schemes = new TreeMap<>();
        for (String resource : BUILTIN_RESOURCES) {
            try (InputStream is = VariantColorSchemes.class.getResourceAsStream(resource)) {
                if (is == null) {
                    log.error("Built in variant color scheme not found: " + resource);
                    continue;
                }
                BufferedReader reader = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8));
                VariantColorScheme scheme = VariantColorScheme.parse(reader, null);
                schemes.put(scheme.getKey(), scheme);
            } catch (IOException e) {
                log.error("Error loading built in variant color scheme " + resource, e);
            }
        }
        return schemes;
    }

    private static String stripExtension(String fileName) {
        int idx = fileName.lastIndexOf('.');
        return idx > 0 ? fileName.substring(0, idx) : fileName;
    }

    /**
     * Discard cached schemes, so the next access rereads them.  For tests and for the preferences editor.
     */
    public static synchronized void reset() {
        userSchemes = null;
        builtinSchemes = null;
        allSchemes = null;
        cacheIgvDirectory = null;
    }
}
