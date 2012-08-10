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
    private String cmdArg;
    private String defaultValue;

    /**
     * Full class name of encoding codec to be used
     */
    private String encodingCodec;

    /**
     * URLs to search for encoding codec, in addition
     * to default
     */
    private URL[] libURLs;

    public String getEncodingCodec() {
        return encodingCodec;
    }

    public enum InputType {
        TEXT,
        FEATURE_TRACK,
        MULTI_FEATURE_TRACK
    }

    Argument(String name, InputType type, String cmdArg, String defaultValue, String encodingCodec,
             String libs, String specPath) {
        this.name = name;
        this.type = type;
        this.cmdArg = cmdArg != null ? cmdArg : "";
        this.defaultValue = defaultValue;
        this.encodingCodec = encodingCodec;
        this.libURLs = FileUtils.getURLsFromString(libs, specPath);
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
        return new Argument(name, type, cmdArg, defVal, encCodec, libString, specPath);
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

}