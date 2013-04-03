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

package org.broad.igv.gs.dm;

import com.google.gson.*;
import org.broad.igv.session.SubtlyImportant;

import java.lang.reflect.Type;

/**
 * Represents a file or directory in GS storage.
 *
 * @author Jim Robinson
 * @date Aug 2, 2011
 */
public class GSFileMetadata {
    private boolean isDirectory;
    private String name;
    private String path;
    private String url;
    private String format;
    private String size;

    public GSFileMetadata(String name, String path, String url, String format, String size, boolean isDirectory) {
        this.isDirectory = isDirectory;
        this.name = name;
        this.path = path;
        this.url = url;
        this.format = format;
        this.size = size;
    }

    @SubtlyImportant
    private GSFileMetadata(){}

//    public GSFileMetadata(Map<String, Object> o){
//        name = (String) o.get("name");
//        path = (String) o.get("path");
//        url = (String) o.get("url");
//        isDirectory = (Boolean) o.get("isDirectory");
//        if (o.containsKey("dataFormat")) {
//            Object dataFormat = o.get("dataFormat");
//            format = dataFormat == null ? "" : (String) ((Map) dataFormat).get("name");
//            size =  o.get("size").toString();
//        }
//    }

    public String toString() {
        return getName();
    }

    public boolean isDirectory() {
        return isDirectory;
    }

    public String getName() {
        return name;
    }

    public String getPath() {
        return path;
    }


    public String getUrl() {
        return url;
    }

    public String getFormat() {
        return format;
    }

    public String getSize() {
        return size;
    }

    public static GSFileMetadata deserializeElement(JsonElement json){
        JsonObject jobj = json.getAsJsonObject();
        String name = jobj.get("name").getAsString();
        String path = jobj.get("path").getAsString();
        String url = jobj.get("url").getAsString();
        boolean isDirectory = jobj.get("isDirectory").getAsBoolean();
        String format = "";
        String size = null;

        if (jobj.has("dataFormat")) {
            JsonObject dataFormat = jobj.get("dataFormat").getAsJsonObject();
            format = dataFormat == null ? "" : dataFormat.get("name").getAsString();
            size = jobj.get("size").getAsString();
        }
        return new GSFileMetadata(name, path, url, format, size, isDirectory);
    }

    private static class Deserializer implements JsonDeserializer<GSFileMetadata> {
        public GSFileMetadata deserialize(JsonElement json, Type typeOfT, JsonDeserializationContext context)
                throws JsonParseException {
            return GSFileMetadata.deserializeElement(json);
        }

    }


}
