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

package org.broad.igv.cli_plugin;

import org.apache.log4j.Logger;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.FeatureTrack;

import javax.xml.bind.annotation.*;
import java.util.List;

@XmlAccessorType(XmlAccessType.NONE)
public class Argument{

    public static final String CMD_ARG = "cmd_arg";
    public static final String LIBS = "libs";

    private static final Logger log = Logger.getLogger(Argument.class);

    @XmlAttribute
    private String name;

    @XmlAttribute
    private InputType type;

    /**
     * Text which goes before Argument value on command line
     */
    @XmlAttribute(name = CMD_ARG)
    private String cmdArg;

    @XmlAttribute(name="default")
    private String defaultValue;

    /**
     * Whether the argument gets written to the command line
     * This is true by default (usually it will be). Some only exist
     * to take user input, and feed the result to another argument
     */
    @XmlAttribute
    private boolean output = false;

    /**
     * id used by spec by which this argument can be referred.
     * Does not need to be human readable, must be unique
     * within a command
     */
    @XmlAttribute
    private String id;

    /**
     * Full class name of encoding codec to be used
     * In addition to default classpath, will
     * also search {@link #libPaths}
     */
    @XmlAttribute
    private String encodingCodec;

    /**
     * URLs to search for encoding codec class, in addition
     * to current classpath
     */
    @XmlElement(name = LIBS)
    private String[] libPaths;

    @XmlEnum
    public enum InputType {
        TEXT,
        LONGTEXT,
        FEATURE_TRACK,
        DATA_TRACK,
        MULTI_FEATURE_TRACK,
        ALIGNMENT_TRACK
    }

    /**
     * Used reading/writing to session files
     */
    @XmlElement
    List<String> value;

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

    public String[] getLibPaths(){
        return this.libPaths;
    }

    public boolean isOutput() {
        return output;
    }

    public String getId() {
        return id;
    }

    public String getEncodingCodec() {
        return encodingCodec;
    }

    @SubtlyImportant
    private Argument(){}

    //This constructor is here largely for testing, should consider getting rid of it
    Argument(String name, InputType type, String cmdArg, String defaultValue, String encodingCodec,
             String[] libPaths, boolean isOutput, String id) {
        this.name = name;
        this.type = type;
        this.cmdArg = cmdArg != null ? cmdArg : "";
        this.defaultValue = defaultValue;
        this.encodingCodec = encodingCodec;
        this.libPaths = libPaths;
        this.output = isOutput;
        this.id = id;

        if (!output && id == null) {
            log.info(String.format("Argument %s is not output but it also has no id. This argument" +
                    " will have no effect", (name)));
        }
    }

}