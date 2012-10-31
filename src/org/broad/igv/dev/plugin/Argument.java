/*
 * Copyright (c) 2007-2012 The Broad Institute, Inc.
 * SOFTWARE COPYRIGHT NOTICE
 * This software and its documentation are the copyright of the Broad Institute, Inc. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. The Broad Institute is not responsible for its use, misuse, or functionality.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 */

package org.broad.igv.dev.plugin;

import org.apache.log4j.Logger;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.util.FileUtils;
import org.broad.igv.util.Utilities;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import java.net.URL;
import java.util.List;

public class Argument {

    private static final Logger log = Logger.getLogger(Argument.class);

    private String name;
    private InputType type;

    /**
     * Text which goes before Argument value on command line
     */
    private String cmdArg;
    private String defaultValue;

    /**
     * Whether the argument gets written to the command line
     * This is true by default (usually it will be). Some only exist
     * to take user input, and feed the result to another argument
     */
    private boolean output;

    /**
     * id used by spec by which this argument can be referred.
     * Does not need to be human readable, must be unique
     * within a command
     */
    private String id;

    /**
     * Full class name of encoding codec to be used
     * In addition to default classpath, will
     * also search {@link #libURLs}
     */
    private String encodingCodec;

    /**
     * URLs to search for encoding codec class, in addition
     * to current classpath
     */
    private URL[] libURLs;

    public String getEncodingCodec() {
        return encodingCodec;
    }

    public enum InputType {
        TEXT,
        LONGTEXT,
        FEATURE_TRACK,
        DATA_TRACK,
        MULTI_FEATURE_TRACK,
        ALIGNMENT_TRACK
    }

    Argument(String name, InputType type, String cmdArg, String defaultValue, String encodingCodec,
             String libs, String specPath, boolean isOutput, String id) {
        this.name = name;
        this.type = type;
        this.cmdArg = cmdArg != null ? cmdArg : "";
        this.defaultValue = defaultValue;
        this.encodingCodec = encodingCodec;
        this.libURLs = FileUtils.getURLsFromString(libs, specPath);
        this.output = isOutput;
        this.id = id;

        if (!output && id == null) {
            log.info(String.format("Argument %s is not output but it also has no id. This argument" +
                    " will have no effect", (name)));
        }
    }

    static Argument parseFromNode(Node node, String specPath) {
        NamedNodeMap attrs = node.getAttributes();
        String name = Utilities.getNullSafe(attrs, "name");
        InputType type = InputType.valueOf(Utilities.getNullSafe(attrs, "type").toUpperCase());
        String cmdArg = Utilities.getNullSafe(attrs, "cmd_arg");
        String defVal = Utilities.getNullSafe(attrs, "default");
        String encCodec = Utilities.getNullSafe(attrs, "encoding_codec");
        String libString = Utilities.getNullSafe(attrs, "libs");
        libString = libString != null ? libString : "";

        String soutput = Utilities.getNullSafe(attrs, "output");
        boolean output = soutput == null || Boolean.parseBoolean(soutput);

        String id = Utilities.getNullSafe(attrs, "id");
        return new Argument(name, type, cmdArg, defVal, encCodec, libString, specPath, output, id);
    }

    boolean isValidValue(Object value) {
        switch (this.type) {
            case TEXT:
                return value instanceof String || value == null;
            case FEATURE_TRACK:
                return value instanceof FeatureTrack;
            case MULTI_FEATURE_TRACK:
                if (!(value instanceof List)) return false;
                try {
                    List<FeatureTrack> lVal = (List<FeatureTrack>) value;
                    FeatureTrack fVal = lVal.get(0);
                    return true;
                } catch (Exception e) {
                    return false;
                }
            default:
                return false;
        }
    }

    public String getName() {
        return name;
    }

    public String getCmdArg() {
        return cmdArg;
    }

    public InputType getType() {
        return type;
    }

    public String getDefaultValue() {
        return defaultValue;
    }

    public URL[] getLibURLs() {
        return libURLs;
    }

    public boolean isOutput() {
        return output;
    }

    public String getId() {
        return id;
    }

}