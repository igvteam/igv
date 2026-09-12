package org.igv.variant;

import org.igv.DirectoryManager;
import org.igv.logging.LogManager;
import org.igv.logging.Logger;
import org.igv.renderer.AbstractColorScale;
import org.igv.util.FileUtils;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
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
        try (BufferedReader reader = new BufferedReader(new java.io.FileReader(file, StandardCharsets.UTF_8))) {
            scheme = VariantColorScheme.parse(reader, stripExtension(file.getName()));
        }

        File directory = getSchemeDirectory();
        if (!directory.equals(file.getParentFile())) {
            if (scheme.getSource() == null) {
                scheme.setSource(file.getAbsolutePath());
            }
            FileUtils.copyFile(file, new File(directory, file.getName()));
        }

        getUserSchemes();   // ensure loaded before adding
        userSchemes.add(0, scheme);
        return scheme;
    }

    /**
     * The directory holding imported schemes.  It is the registry -- every file in it is loaded at startup.
     */
    public static File getSchemeDirectory() {
        File directory = new File(DirectoryManager.getIgvDirectory(), SCHEME_DIRECTORY);
        if (!directory.exists()) {
            directory.mkdir();
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
            try (BufferedReader reader = new BufferedReader(new java.io.FileReader(file, StandardCharsets.UTF_8))) {
                schemes.add(VariantColorScheme.parse(reader, stripExtension(file.getName())));
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
