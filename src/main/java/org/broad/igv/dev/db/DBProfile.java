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

package org.broad.igv.dev.db;

import org.apache.log4j.Logger;
import org.broad.igv.session.SubtlyImportant;
import org.broad.igv.util.ParsingUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.w3c.dom.Document;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.*;
import javax.xml.bind.annotation.adapters.XmlAdapter;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.Method;
import java.util.*;

/**
 * Object representation of the location of a database,
 * specifying address, tables, etc.
 * User: jacob
 * Date: 2013-Jan-14
 */
@XmlRootElement(name = "database")
@XmlAccessorType(XmlAccessType.NONE)
@XmlSeeAlso(DBProfile.DBTable.class)
public class DBProfile {

    private static Logger log = Logger.getLogger(DBProfile.class);

    @XmlAttribute private String name;
    @XmlAttribute private String description;
    @XmlAttribute private String version;
    @XmlAttribute private String subprotocol;
    @XmlAttribute private String host;
    @XmlAttribute private String path;
    @XmlAttribute private String port;
    @XmlAttribute private String username;
    @XmlAttribute private String password;

    @XmlElement(name = "table")
    private List<DBTable> tableList;

    public List<DBTable> getTableList() {
        return tableList;
    }

    @XmlTransient
    private ResourceLocator dbLocator;

    public ResourceLocator getDBLocator() {
        if(dbLocator == null){
            dbLocator = new ResourceLocator(DBManager.createConnectionURL(subprotocol, host, path, port));
            dbLocator.setUsername(username);
            dbLocator.setPassword(password);
        }
        return dbLocator;
    }

    private static JAXBContext jc = null;
    public static JAXBContext getJAXBContext() throws JAXBException {
        if(jc == null){
            jc = JAXBContext.newInstance(DBProfile.class);
        }
        return jc;
    }

    public static DBProfile parseProfile(String profilePath){

        InputStream profileStream = null;
        try {
            profileStream = ParsingUtils.openInputStream(profilePath);
        } catch (IOException e) {
            try {
                if (profileStream != null) profileStream.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            throw new RuntimeException("Unable to open DB profile", e);
        }

        Document doc = null;
        try{
            doc = Utilities.createDOMDocumentFromXmlStream(profileStream);
        }catch (Exception e){
            throw new RuntimeException("Error parsing DB Profile", e);
        }

        try{
        JAXBContext jc = getJAXBContext();
        Unmarshaller u = jc.createUnmarshaller();

//        u.setListener(new Unmarshaller.Listener() {
//            @Override
//            public void afterUnmarshal(Object target, Object parent) {
//                if(target instanceof DBTable && parent instanceof DBProfile){
//                    ((DBTable) target).setDbLocator(((DBProfile) parent).getDBLocator());
//                }
//            }
//        });

        return u.unmarshal(doc, DBProfile.class).getValue();
        }catch(JAXBException e){
            throw new RuntimeException("Error unmarshalling DB Profile, it may be misformed", e);
        }
    }

    @SubtlyImportant
    public void afterUnmarshal(Unmarshaller u, Object parent) {
        dbLocator = getDBLocator();
        for (DBTable table : getTableList()) {
            table.setDbLocator(dbLocator);
        }
    }

    public String getHost() {
        return host;
    }

    public String getPassword() {
        return password;
    }

    public String getPort() {
        return port;
    }

    public String getSubprotocol() {
        return subprotocol;
    }

    public String getUsername() {
        return username;
    }

    public String getPath() {
        return path;
    }

    public String getName() {
        return name;
    }

    public void setHost(String host) {
        this.host = host;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setPassword(String password) {
        this.password = password;
    }

    public void setPath(String path) {
        this.path = path;
    }

    public void setPort(String port) {
        this.port = port;
    }

    public void setSubprotocol(String subprotocol) {
        this.subprotocol = subprotocol;
    }

    public void setUsername(String username) {
        this.username = username;
    }

    void addTable(DBTable newTable) {
        if(tableList == null){
            tableList = new ArrayList<DBTable>();
        }
        tableList.add(newTable);
    }


    /**
     * Checks for {@code name}, {@code host}, {@code path},
     * and {@code username}
     * See {@link #checkMissingValues(Object, String[]}
     * @return
     */
    public List<String> checkMissingValues() {
        String[] requiredGetters = new String[]{"getName", "getHost", "getPath", "getUsername"};
        return checkMissingValues(this, requiredGetters);
    }

    /**
     * Checks {@code object} for non-null values (or empty strings)
     * retrieved by getter in {@code requiredGetters}. If the getter is not found or cannot
     * be accessed, an error is logged (access restrictions are followed).
     * @param object The object to check for missing values
     * @param requiredGetters The getter method names to use. Must take no arguments
     * @return
     */
    public static List<String> checkMissingValues(Object object, String[] requiredGetters){
        List<String> missing = new ArrayList<String>(requiredGetters.length);
        for(String reqGetter: requiredGetters){
            Method method = null;
            try {
                method = object.getClass().getMethod(reqGetter);
                Object value = method.invoke(object);
                if(value instanceof String && ((String) value).length() == 0){
                    value = null;
                }
                if(value == null) missing.add(reqGetter);
            } catch (Exception e) {
                log.error(e.getMessage(), e);
                missing.add(reqGetter);
            }
        }
        return missing;
    }

    /**
     * Object representation of a single {@code table} element of
     * a database profile. Contains static method for parsing dbXML files
     * <p/>
     * User: jacob
     * Date: 2012-Oct-31
     */
    @XmlAccessorType(XmlAccessType.NONE)
    public static class DBTable {

        private static Logger log = Logger.getLogger(DBTable.class);

        @XmlAttribute private String name;
        @XmlAttribute private String format;
        @XmlAttribute private String binColName;

        @XmlAttribute private String chromoColName;
        @XmlAttribute private String posStartColName;
        @XmlAttribute private String posEndColName;
        @XmlAttribute private int startColIndex = 1;
        @XmlAttribute private int endColIndex = Integer.MAX_VALUE;
        @XmlAttribute private String baseQuery;

        @XmlElement(name = "column")
        private List<Column> columnList;

        @XmlElement(name = "header")
        private List<String> headerLines;

        @XmlTransient
        private HashMap<Integer, String> columnLabelMap;

        @XmlTransient
        private ResourceLocator dbLocator;

        /**
         * Generally just intended for testing, where all we try
         * to do is get all data from a db and don't need anything fancy
         *
         * @param dbLocator
         * @param tableName
         */
        public static DBTable build(ResourceLocator dbLocator, String tableName) {
            return new DBTable(dbLocator, tableName, null, null, null, null, null, 1, Integer.MAX_VALUE - 1, null, null, null);
        }

        private DBTable(){}

        DBTable(ResourceLocator dbLocator, String name){
            this.dbLocator = dbLocator;
            this.name = name;
        }

        public DBTable(ResourceLocator dbLocator, String name, String format, String binColName,
                       String chromoColName, String posStartColName, String posEndColName, int startColIndex, int endColIndex,
                       HashMap<Integer, String> columnLabelMap, String baseQuery, List<String> headerLines) {
            this.dbLocator = dbLocator;
            this.name = name;
            this.format = format;
            this.binColName = binColName;
            this.chromoColName = chromoColName;
            this.posStartColName = posStartColName;
            this.posEndColName = posEndColName;
            this.startColIndex = startColIndex;
            this.endColIndex = endColIndex;
            this.columnLabelMap = columnLabelMap;
            this.baseQuery = baseQuery;
            this.headerLines = headerLines;
        }

        @SubtlyImportant
        public void afterUnmarshal(Unmarshaller u, Object parent) {
            if(columnList != null){
                columnLabelMap = ColumnMapAdapter.unmarshal(columnList);
            }
        }

        public List<String> checkMissingValues(){
            String[] requiredGetters = new String[]{"getFormat", "getChromoColName", "getPosStartColName", "getPosEndColName"};
            return DBProfile.checkMissingValues(this, requiredGetters);
        }

        public ResourceLocator getDbLocator() {
            return dbLocator;
        }

        public String getBinColName() {
            return binColName;
        }

        public String getChromoColName() {
            return chromoColName;
        }

        public int getEndColIndex() {
            return endColIndex;
        }

        public String getFormat() {
            return format;
        }

        public String getPosEndColName() {
            return posEndColName;
        }

        public String getPosStartColName() {
            return posStartColName;
        }

        public int getStartColIndex() {
            return startColIndex;
        }

        public String getName() {
            return name;
        }

        public String getBaseQuery() {
            return baseQuery;
        }

        public Map<Integer, String> getColumnLabelMap() {
            return columnLabelMap;
        }

        public List<String> getHeaderLines() {
            return headerLines;
        }

        void setDbLocator(ResourceLocator dbLocator) {
            this.dbLocator = dbLocator;
        }

        public void setBinColName(String binColName) {
            this.binColName = binColName;
        }

        public void setChromoColName(String chromoColName) {
            this.chromoColName = chromoColName;
        }

        public void setEndColIndex(int endColIndex) {
            this.endColIndex = endColIndex;
        }

        public void setFormat(String format) {
            this.format = format;
        }

        public void setPosEndColName(String posEndColName) {
            this.posEndColName = posEndColName;
        }

        public void setPosStartColName(String posStartColName) {
            this.posStartColName = posStartColName;
        }

        public void setStartColIndex(int startColIndex) {
            this.startColIndex = startColIndex;
        }

        /**
         * Return an array of column labels in specified ordinal positions
         *
         * @return
         */
        public static String[] columnMapToArray(Map<Integer, String> columnLabelMap) {
            List<Integer> arrayIndexes = new ArrayList<Integer>(columnLabelMap.keySet());
            Collections.sort(arrayIndexes);
            int minArrayIndex = arrayIndexes.get(0);
            int maxArrayIndex = arrayIndexes.get(arrayIndexes.size() - 1);
            int colCount = maxArrayIndex + 1;
            String[] tokens = new String[colCount];

            for (int cc = minArrayIndex; cc <= maxArrayIndex; cc++) {
                tokens[cc] = columnLabelMap.get(cc);
            }
            return tokens;
        }

        private static class ColumnMapAdapter extends XmlAdapter<XmlMap, HashMap<Integer, String>> {
            @Override
            public HashMap<Integer, String> unmarshal(XmlMap v){
                return unmarshal(v.column);
            }

            public static HashMap<Integer, String> unmarshal(List<Column> columnList){
                HashMap<Integer, String> result = new HashMap<Integer, String>(columnList.size());
                for(Column entry: columnList){
                    result.put(entry.fileIndex, entry.colLabel);
                }
                return result;
            }

            @Override
            public XmlMap marshal(HashMap<Integer, String> v){
                XmlMap result = new XmlMap();
                for(Map.Entry<Integer, String> entry: v.entrySet()){
                    Column mapEntry = new Column();
                    mapEntry.fileIndex = entry.getKey();
                    mapEntry.colLabel = entry.getValue();
                    result.column.add(mapEntry);
                }
                return result;
            }
        }

        private static class XmlMap {
            public List<Column> column =
                    new ArrayList<Column>();
        }

        private static class Column {
            @XmlAttribute private Integer fileIndex;
            @XmlAttribute private String colLabel;
        }
    }
}
