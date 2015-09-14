/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2007-2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package org.broad.igv.cli_plugin;

import org.apache.log4j.Logger;
import org.broad.igv.sam.AlignmentTrack;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.track.FeatureTrack;
import org.broad.igv.variant.VariantTrack;

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

    @XmlAttribute
    private String defaultValue;

    /**
     * Whether the argument gets written to the command line
     * This is true by default (usually it will be). Some only exist
     * to take user input, and feed the result to another argument
     */
    @XmlAttribute
    private boolean output = false;

    /**
     * Whether this argument is displayed to the user.
     * If we just want to have a value in the XML spec
     * that the user can't change, set this to false
     */
    @XmlAttribute
    private boolean visible = true;

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

    /**
     * Key used in paths to indicate same directory as tool
     */
    public static final CharSequence TOOL_DIR_KEY = "$toolDir";

    @XmlAttribute
    private boolean remembered;

    public boolean isRemembered() {
        return remembered;
    }

    @XmlEnum
    public enum InputType {
        BOOL,
        TEXT,
        LONGTEXT,
        FEATURE_TRACK,
        DATA_TRACK,
        MULTI_FEATURE_TRACK,
        ALIGNMENT_TRACK,
        VARIANT_TRACK,
        LOCUS
    }

    /**
     * Used reading/writing to session files
     */
    @XmlElement
    List<String> value;

    boolean isValidValue(Object value) {
        switch (this.type) {
            case BOOL:
                return value instanceof Boolean;
            case TEXT:
            case LONGTEXT:
                return value instanceof String || value == null;
            case VARIANT_TRACK:
                return value instanceof VariantTrack;
            case ALIGNMENT_TRACK:
                return value instanceof AlignmentTrack;
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

    public void setDefaultValue(String defaultValue){
        this.defaultValue = defaultValue;
    }

    public String[] getLibPaths(){
        return this.libPaths;
    }

    public boolean isOutput() {
        return output;
    }

    public void setOutput(boolean output){
        this.output = output;
    }

    public String getId() {
        return id;
    }

    public String getEncodingCodec() {
        return encodingCodec;
    }

    public boolean isVisible() {
        return visible;
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