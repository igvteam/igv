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
