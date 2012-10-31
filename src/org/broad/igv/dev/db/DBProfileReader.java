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

package org.broad.igv.dev.db;

import org.broad.igv.feature.genome.GenomeManager;
import org.broad.igv.feature.tribble.CodecFactory;
import org.broad.igv.ui.util.MessageUtils;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.util.Utilities;
import org.broad.tribble.AsciiFeatureCodec;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXParseException;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * User: jacob
 * Date: 2012-Oct-31
 */
public class DBProfileReader {

    /**
     * Retrieve a reader from the XML profile located at {@code profilePath}.
     * TODO If {@code tableName == null}, the user is prompted to choose a table from the list
     *
     * @param profilePath
     * @param tableName
     * @return
     */
    public static List<SQLCodecSource> getFromProfile(String profilePath, String tableName) {
        ResourceLocator dbLocator = DBManager.getStoredConnection(profilePath);
        InputStream profileStream = null;
        try {
            profileStream = new FileInputStream(profilePath);
            Document document = Utilities.createDOMDocumentFromXmlStream(profileStream);
            NodeList tableNodes = document.getElementsByTagName("table");
            List<SQLCodecSource> sources = new ArrayList<SQLCodecSource>(tableNodes.getLength());

            for (int tnum = 0; tnum < tableNodes.getLength(); tnum++) {
                Node tableNode = tableNodes.item(tnum);
                NamedNodeMap attr = tableNode.getAttributes();
                String tabName = attr.getNamedItem("name").getTextContent();
                if (tableName == null || tableName.equals(tabName)) {

                    String chromoColName = attr.getNamedItem("chromoColName").getTextContent();
                    String posStartColName = attr.getNamedItem("posStartColName").getTextContent();
                    String posEndColName = attr.getNamedItem("posEndColName").getTextContent();
                    String format = attr.getNamedItem("format").getTextContent();
                    String startColString = Utilities.getNullSafe(attr, "startColIndex");
                    String endColString = Utilities.getNullSafe(attr, "endColIndex");
                    String binColName = Utilities.getNullSafe(attr, "binColName");

                    int startColIndex = Integer.parseInt(startColString);
                    int endColIndex = Integer.parseInt(endColString);

                    //If present, retrieve list of columns
                    NodeList columns = tableNode.getChildNodes();
                    DBReader.ColumnMap columnMap = null;

                    if (columns.getLength() > 0) {
                        columnMap = new DBReader.ColumnMap();

                        for (int col = 0; col < columns.getLength(); col++) {
                            Node column = columns.item(col);
                            NamedNodeMap colAttr = column.getAttributes();
                            int fileIndex = Integer.parseInt(colAttr.getNamedItem("fileIndex").getTextContent());

                            String colLabel = Utilities.getNullSafe(colAttr, "colLabel");
                            String colIndexStr = Utilities.getNullSafe(colAttr, "colIndex");

                            if (colLabel == null && colIndexStr == null) {
                                String msg = String.format("Error parsing column %d of %s", col, tabName);
                                msg += "colName and colIndex both null, at least 1 must be specified.";
                                throw new SAXParseException(msg, null);
                            } else if (colIndexStr != null) {
                                int colIndex = Integer.parseInt(colIndexStr);
                                columnMap.put(fileIndex, colIndex);
                            }
                            columnMap.put(fileIndex, colLabel);
                        }
                    }

                    AsciiFeatureCodec codec = CodecFactory.getCodec("." + format, GenomeManager.getInstance().getCurrentGenome());
                    SQLCodecSource source = new SQLCodecSource(dbLocator, tabName, codec, columnMap,
                            chromoColName, posStartColName, posEndColName, startColIndex, endColIndex);

                    source.binColName = binColName;
                    sources.add(source);
                }
            }

            return sources;

        } catch (Exception e) {
            String msg = "Error reading profile " + profilePath + ", table " + tableName;
            MessageUtils.showErrorMessage(msg, e);
            return null;
        } finally {
            try {
                if (profileStream != null) profileStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

    }
}
