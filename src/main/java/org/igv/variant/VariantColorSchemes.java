package org.igv.variant;

import org.igv.DirectoryManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.AbstractColorScale;
import org.igv.renderer.ContinuousColorScale;
import org.igv.util.FileUtils;

import java.awt.Color;
import org.igv.ui.color.ColorUtilities;
import org.igv.util.HttpUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringReader;
import java.net.URL;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * The color schemes available for coloring variants by a VCF INFO attribute.
 * <p>
 * Schemes are files, not individual preferences -- the set of INFO keys, and of values for each key, is unbounded,
 * so there is nothing to enumerate in a preferences grid.  Importing a scheme copies it into the IGV directory
 * ("variantcolors"), which IGV then owns: the directory is the registry, scanned at startup, in the same way gene
 * lists are managed (see {@link org.igv.lists.GeneListManager}).  Schemes shipped with IGV are read from the
 * classpath and are never copied, so upgrades can revise them.
 * <p>
 * User schemes are searched before the built-in ones, so a user scheme covering CLNSIG shadows IGV's.
 */
public class VariantColorSchemes {

    private static Logger log = LogManager.getLogger(VariantColorSchemes.class);

    static final String SCHEME_DIRECTORY = "variantcolors";

    private static final String BUILTIN_RESOURCE = "resources/variant_colors.txt";

    private static List<VariantColorScheme> userSchemes;
    private static List<VariantColorScheme> builtinSchemes;

    private VariantColorSchemes() {
    }

    /**
     * @return the user imported schemes, followed by those shipped with IGV.
     */
    public static synchronized List<VariantColorScheme> getSchemes() {
        List<VariantColorScheme> schemes = new ArrayList<>(getUserSchemes());
        schemes.addAll(getBuiltinSchemes());
        return schemes;
    }

    public static synchronized List<VariantColorScheme> getUserSchemes() {
        if (userSchemes == null) {
            userSchemes = loadUserSchemes();
        }
        return Collections.unmodifiableList(userSchemes);
    }

    public static synchronized List<VariantColorScheme> getBuiltinSchemes() {
        if (builtinSchemes == null) {
            builtinSchemes = loadBuiltinSchemes();
        }
        return Collections.unmodifiableList(builtinSchemes);
    }

    /**
     * Return the color assigned to an INFO attribute value, or null if no scheme covers it, in which case the
     * caller assigns a color from a palette.
     */
    public static Color getColor(String infoKey, String value) {
        for (VariantColorScheme scheme : getSchemes()) {
            Color color = scheme.getColor(infoKey, value);
            if (color != null) {
                return color;
            }
        }
        return null;
    }

    /**
     * Return the continuous scale for a numeric INFO attribute, or null if no scheme defines one.  Numeric
     * attributes are only colorable if a scheme gives them a range -- there is no sensible default.
     */
    public static AbstractColorScale getScale(String infoKey) {
        for (VariantColorScheme scheme : getSchemes()) {
            AbstractColorScale scale = scheme.getScale(infoKey);
            if (scale != null) {
                return scale;
            }
        }
        return null;
    }

    /**
     * @return true if some scheme declares this numeric attribute to hold categories rather than quantities.
     */
    public static boolean isCategorical(String infoKey) {
        for (VariantColorScheme scheme : getSchemes()) {
            if (scheme.isCategorical(infoKey)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return the scheme providing the color scale for an INFO attribute, or null if none does.
     */
    public static VariantColorScheme getSchemeForScale(String infoKey) {
        for (VariantColorScheme scheme : getSchemes()) {
            if (scheme.getScale(infoKey) != null) {
                return scheme;
            }
        }
        return null;
    }

    /**
     * @return the colors any scheme assigns to values of an INFO attribute.  Used to keep colors assigned from
     * the palette distinguishable from them.
     */
    public static Collection<Color> getColors(String infoKey) {
        List<Color> colors = new ArrayList<>();
        for (VariantColorScheme scheme : getSchemes()) {
            colors.addAll(scheme.getColors(infoKey).values());
        }
        return colors;
    }

    /**
     * @return every INFO key covered by some scheme.  Used to offer numeric attributes that would otherwise not
     * be colorable.
     */
    public static Set<String> getKeys() {
        Set<String> keys = new LinkedHashSet<>();
        for (VariantColorScheme scheme : getSchemes()) {
            keys.addAll(scheme.getKeys());
        }
        return keys;
    }

    /**
     * Import a color scheme.  The file is copied into the IGV directory unless it is already there, so the
     * scheme survives the original being moved or deleted, and is available on the next startup.
     *
     * @return the imported scheme
     */
    public static synchronized VariantColorScheme importFile(File file) throws IOException {

        VariantColorScheme scheme;
        try (BufferedReader reader = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8))) {
            scheme = VariantColorScheme.parse(reader, stripExtension(file.getName()));
        }

        File directory = createSchemeDirectory();
        if (directory.equals(file.getParentFile())) {
            scheme.setFile(file);
        } else {
            if (scheme.getSource() == null) {
                scheme.setSource(file.getAbsolutePath());
            }
            File copy = new File(directory, file.getName());
            FileUtils.copyFile(file, copy);
            scheme.setFile(copy);
        }

        register(scheme);
        return scheme;
    }

    /**
     * Import a color scheme from a URL.  The contents are fetched once and saved in the IGV directory, so the
     * scheme is available offline and on the next startup.  The URL is recorded so it can be fetched again.
     *
     * @return the imported scheme
     */
    public static synchronized VariantColorScheme importUrl(String url) throws IOException {

        String contents = HttpUtils.getInstance().getContentsAsString(new URL(url));

        VariantColorScheme scheme = VariantColorScheme.parse(
                new BufferedReader(new StringReader(contents)), stripExtension(fileNameFromUrl(url)));
        if (scheme.getSource() == null) {
            scheme.setSource(url);
        }

        File file = new File(createSchemeDirectory(), fileNameFromUrl(url));
        try (PrintWriter writer = new PrintWriter(file, StandardCharsets.UTF_8)) {
            if (!contents.toLowerCase().contains("#source=")) {
                writer.println("#source=" + url);
            }
            writer.print(contents);
        }
        scheme.setFile(file);

        register(scheme);
        return scheme;
    }

    /**
     * Write the given value -> color assignments as a new scheme in the IGV directory, so they apply to every VCF
     * with this attribute rather than only the track they were chosen on.  An existing scheme of the same name is
     * replaced.
     *
     * @param name     the scheme name, also the basis for the file name
     * @param infoKey  the INFO attribute the colors are for
     * @param colors   value -> color
     * @return the saved scheme
     */
    public static synchronized VariantColorScheme saveScheme(String name, String infoKey, Map<String, Color> colors)
            throws IOException {
        return saveScheme(name, infoKey, colors, false);
    }

    /**
     * @param categorical true to record that this numeric attribute holds categories rather than quantities.  The
     *                    declaration is what persists the answer -- the colors only cover the values seen so far.
     */
    public static synchronized VariantColorScheme saveScheme(String name, String infoKey, Map<String, Color> colors,
                                                             boolean categorical) throws IOException {

        VariantColorScheme scheme = new VariantColorScheme(name);
        if (categorical) {
            scheme.setCategorical(infoKey);
        }
        for (Map.Entry<String, Color> entry : colors.entrySet()) {
            scheme.setColor(infoKey, entry.getKey(), entry.getValue());
        }
        return save(scheme);
    }

    /**
     * Write a scheme to the IGV directory and make it current.  A scheme shipped with IGV cannot be written, so
     * saving one writes a copy the user owns, which then shadows it.
     *
     * @return the saved scheme
     */
    public static synchronized VariantColorScheme save(VariantColorScheme scheme) throws IOException {

        File file = scheme.getFile();
        if (file == null) {
            file = new File(createSchemeDirectory(), getLegalFileName(scheme.getName()) + ".txt");
        } else {
            createSchemeDirectory();
        }

        try (PrintWriter writer = new PrintWriter(file, StandardCharsets.UTF_8)) {
            scheme.write(writer);
        }
        scheme.setFile(file);

        register(scheme);
        return scheme;
    }

    /**
     * Write a color scale for a numeric INFO attribute as a scheme, in the same human editable form as the rest
     * of the file: a "min:max" or "min:mid:max" range followed by a color per stop.
     *
     * @return the saved scheme
     */
    public static synchronized VariantColorScheme saveScale(String name, String infoKey, ContinuousColorScale scale)
            throws IOException {

        VariantColorScheme scheme = new VariantColorScheme(name);
        scheme.setScale(infoKey, scale);
        return save(scheme);
    }

    /**
     * Delete an imported scheme.  Schemes shipped with IGV cannot be removed.
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
        getUserSchemes();
        userSchemes.remove(scheme);
        return true;
    }

    /**
     * Add a scheme to the head of the user list, replacing any scheme already loaded from the same file.
     */
    private static void register(VariantColorScheme scheme) {
        getUserSchemes();   // ensure loaded before adding
        userSchemes.removeIf(s -> scheme.getFile().equals(s.getFile()));
        userSchemes.add(0, scheme);
    }

    /**
     * The last path segment of a URL, defaulting to something usable if there isn't one.
     */
    private static String fileNameFromUrl(String url) {
        String path = url;
        int query = path.indexOf('?');
        if (query > 0) {
            path = path.substring(0, query);
        }
        int slash = path.lastIndexOf('/');
        String name = slash < 0 ? path : path.substring(slash + 1);
        name = name.replaceAll("[^A-Za-z0-9._-]", "_");
        return name.isEmpty() || name.startsWith(".") ? "colors.txt" : name;
    }

    /**
     * The directory holding imported schemes.  It is the registry -- every file in it is loaded at startup.  It is
     * not created here, only when a scheme is actually imported.
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

    private static List<VariantColorScheme> loadUserSchemes() {

        List<VariantColorScheme> schemes = new ArrayList<>();
        File directory = new File(DirectoryManager.getIgvDirectory(), SCHEME_DIRECTORY);
        File[] files = directory.isDirectory() ? directory.listFiles() : null;
        if (files == null) {
            return schemes;
        }

        for (File file : files) {
            if (file.isDirectory() || file.isHidden()) {
                continue;
            }
            try (BufferedReader reader = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8))) {
                VariantColorScheme scheme = VariantColorScheme.parse(reader, stripExtension(file.getName()));
                scheme.setFile(file);
                schemes.add(scheme);
            } catch (IOException e) {
                log.error("Error loading variant color scheme: " + file.getAbsolutePath(), e);
            }
        }
        return schemes;
    }

    private static List<VariantColorScheme> loadBuiltinSchemes() {

        List<VariantColorScheme> schemes = new ArrayList<>();
        try (InputStream is = VariantColorSchemes.class.getResourceAsStream(BUILTIN_RESOURCE)) {
            if (is == null) {
                log.error("Built in variant color scheme not found: " + BUILTIN_RESOURCE);
            } else {
                BufferedReader reader = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8));
                schemes.add(VariantColorScheme.parse(reader, "IGV defaults"));
            }
        } catch (IOException e) {
            log.error("Error loading built in variant color scheme", e);
        }
        return schemes;
    }

    /**
     * Encode a scheme name so it can be used as a file name, as gene lists do.
     */
    private static String getLegalFileName(String name) {
        return URLEncoder.encode(name, StandardCharsets.UTF_8);
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
    }
}
