/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.  All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
 */
package org.broad.igv.h5;

//~--- non-JDK imports --------------------------------------------------------

import org.apache.http.HttpResponse;
import org.apache.log4j.Logger;
import org.broad.igv.util.IGVHttpClientUtils;
import org.broad.igv.util.IGVHttpUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.*;
import java.net.URL;
import java.net.URLConnection;
import java.util.*;

/**
 * @author jrobinso
 */
public class HDF5RemoteReader implements HDF5Reader {

    /**
     * Map of server -> list of open files
     */
    static Map<String, List<Integer>> allFileIds = new HashMap();

    /**
     * Field description
     */
    public static final int OBJECT_NOT_FOUND_CODE = -1;
    static Logger log = Logger.getLogger(HDF5RemoteReader.class);
    ResourceLocator locator;

    /**
     * Field description
     */
    public int fileId;
    private Map<String, Integer> entityCache = new HashMap(1000);

    /**
     * Method description
     */
    public static void shutdown() {

        // TODO -- send server message to close all file IDs.
        String fileIdString = new String();
        for (String server : allFileIds.keySet()) {
            for (Integer fileId : allFileIds.get(server)) {
                fileIdString += fileId.toString() + ",";
            }

            try {
                URL url = new URL(server + "?method=closeFiles&fileId=" + fileIdString);
                IGVHttpClientUtils.executeGet(url);

            }
            catch (IOException ex) {
                log.error("Error closing files on server: " + server);
            }
        }
    }

    /**
     * Record am HDF5 file handle and associate it with a server.  HDF5 files
     * remain open on the server through the life of a user session.  A message
     * is sent to the server on shutdown to close them.
     *
     * @param server
     * @param fileId
     */
    private synchronized void recordFile(String server, Integer fileId) {
        List<Integer> fileIds = allFileIds.get(server);
        if (fileIds == null) {
            fileIds = new ArrayList();
            allFileIds.put(server, fileIds);
        }
        fileIds.add(fileId);
    }

    /**
     * Constructs ...
     *
     * @param locator
     */
    public HDF5RemoteReader(ResourceLocator locator) {
        InputStream urlStream = null;
        try {
            this.locator = locator;
            URL url = new URL(locator.getServerURL() + "?method=openFile&file="
                    + locator.getPath());
            urlStream = IGVHttpClientUtils.openConnectionStream(url);
            DataInputStream is = new DataInputStream(new BufferedInputStream(urlStream));
            fileId = is.readInt();
            recordFile(locator.getServerURL(), fileId);

        }
        catch (IOException ex) {
            log.error("Error opening file", ex);
            String exceptionMessage = ex.getMessage();
            String message = ex.toString();
            if (exceptionMessage != null) {
                if (exceptionMessage.contains("response code: 403")) {
                    message = "Access forbidden";
                } else if (exceptionMessage.contains("response code: 404")) {
                    message = "Resource not found";
                } else {
                    message = exceptionMessage;
                }
            }
            log.error("Error loading remote h5 file", ex);

            throw new DataAccessException(message);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Method description
     */
    public void closeFile() {
        final String method = "closeFile";
        closeEntity(method, fileId);
    }

    /**
     * Method description
     *
     * @param name
     * @return
     */
    public int openDataset(String name) {
        final String method = "openDataset";
        return openEntity(method, fileId, name);
    }

    /**
     * Method description
     *
     * @param datasetId
     */
    public void closeDataset(int datasetId) {
        final String closeEntity = "closeDataset";
        closeEntity(closeEntity, datasetId);
    }

    /**
     * Method description
     *
     * @param groupPath
     * @return
     */
    public int openGroup(String groupPath) {
        final String method = "openGroup";
        return openEntity(method, fileId, groupPath);
    }

    /**
     * Method description
     *
     * @param groupId
     */
    public void closeGroup(int groupId) {
        final String method = "closeGroup";
        closeEntity(method, groupId);
    }

    /**
     * Method description
     *
     * @param groupId
     * @param attrName
     * @return
     * @throws ObjectNotFoundException
     */
    public int readIntegerAttribute(int groupId, String attrName) throws ObjectNotFoundException {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=readIntegerAttribute&id=" + groupId
                    + "&name=" + attrName);
            urlStream = IGVHttpClientUtils.openConnectionStream(url);
            DataInputStream is = new DataInputStream(new BufferedInputStream(urlStream));
            int value = is.readInt();
            is.close();

            //
            if (value == OBJECT_NOT_FOUND_CODE) {
                throw new ObjectNotFoundException("Attribute: " + attrName + " not found");
            }

            return value;
        }
        catch (IOException ex) {

            // log.error("Error in readIntegerAttribute", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Method description
     *
     * @param groupId
     * @param attrName
     * @return
     */
    public String readStringAttribute(int groupId, String attrName) {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=readStringAttribute&id=" + groupId + "&name=" + attrName);
            urlStream = IGVHttpClientUtils.openConnectionStream(url);

            BufferedReader reader = new BufferedReader(new InputStreamReader(urlStream));
            String value = reader.readLine().replace((char) 711, ' ').trim();
            reader.close();
            return value.length() == 0 ? null : value;
        }
        catch (IOException ex) {

            // log.error("Error in readStringAttribute", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Method description
     *
     * @param groupId
     * @param attrName
     * @return
     */
    public double readDoubleAttribute(int groupId, String attrName) {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=readDoubleAttribute&id=" + groupId
                    + "&name=" + attrName);

            urlStream = IGVHttpClientUtils.openConnectionStream(url);

            DataInputStream dis = new DataInputStream(new BufferedInputStream(urlStream));
            double value = dis.readDouble();
            dis.close();
            return value;
        }
        catch (IOException ex) {
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Method description
     *
     * @param datasetId
     * @return
     */
    public float[] readAllFloats(int datasetId) {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=readAllFloats&id=" + datasetId);

            HttpResponse response = IGVHttpClientUtils.executeGet(url);
            int byteLength = (int) IGVHttpClientUtils.getContentLength(response);

            // Compute # of floats in the stream.
            int nFloats = byteLength / 4;

            float[] data = new float[nFloats];

            urlStream = response.getEntity().getContent();

            DataInputStream is = new DataInputStream(new BufferedInputStream(urlStream));
            for (int i = 0; i < nFloats; i++) {
                data[i] = is.readFloat();
            }
            is.close();
            return data;
        }
        catch (IOException ex) {
            log.error("Error in readAllFloats", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Method description
     *
     * @param datasetId
     * @return
     */
    public int[] readAllInts(int datasetId) {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=readAllInts&id=" + datasetId);

            // Compute the # of ints in the stream.
            HttpResponse response = IGVHttpClientUtils.executeGet(url);
            int byteLength = (int) IGVHttpClientUtils.getContentLength(response);
            int nInts = byteLength / 4;

            int[] data = new int[nInts];

            urlStream = response.getEntity().getContent();
            DataInputStream is = new DataInputStream(new BufferedInputStream(urlStream));
            for (int i = 0; i < nInts; i++) {
                data[i] = is.readInt();
            }
            is.close();
            return data;
        }
        catch (IOException ex) {
            log.error("Error in readAllFloats", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Method description
     *
     * @param datasetId
     * @return
     */
    public List<String> readAllStrings(int datasetId) {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=readAllStrings&id=" + datasetId);
            ArrayList<String> strings = new ArrayList(100);

            urlStream = IGVHttpClientUtils.openConnectionStream(url);
            BufferedReader reader = new BufferedReader(new InputStreamReader(urlStream));

            String nextLine = "";
            while ((nextLine = reader.readLine()) != null) {
                String nextString = nextLine.trim();
                if (nextString.length() > 0) {
                    strings.add(nextString);
                }
            }
            reader.close();

            return strings;

        }
        catch (IOException ex) {
            log.error("Error in readAllFloats", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Method description
     *
     * @param datasetId
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public float[] readFloats(int datasetId, int fromIndex, int toIndex) {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=readFloats&id=" + datasetId
                    + "&from=" + fromIndex + "&to=" + toIndex);

            HttpResponse response = IGVHttpClientUtils.executeGet(url);
            int byteLength = (int) IGVHttpClientUtils.getContentLength(response);
            int nFloats = byteLength / 4;

            float[] data = new float[nFloats];

            urlStream = response.getEntity().getContent();
            DataInputStream is = new DataInputStream(new BufferedInputStream(urlStream));
            for (int i = 0; i < nFloats; i++) {
                data[i] = is.readFloat();
            }
            is.close();
            return data;
        }
        catch (IOException ex) {
            log.error("Error in readAllFloats", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Method description
     *
     * @param datasetId
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public int[] readInts(int datasetId, int fromIndex, int toIndex) {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=readInts&id=" + datasetId
                    + "&from=" + fromIndex + "&to=" + toIndex);

            HttpResponse response = IGVHttpClientUtils.executeGet(url);
            int byteLength = (int) IGVHttpClientUtils.getContentLength(response);
            int nInts = byteLength / 4;

            int[] data = new int[nInts];

            urlStream = response.getEntity().getContent();
            DataInputStream is = new DataInputStream(new BufferedInputStream(urlStream));
            for (int i = 0; i < nInts; i++) {
                data[i] = is.readInt();
            }
            is.close();
            return data;
        }
        catch (IOException ex) {
            log.error("Error in readAllFloats", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /*
     *
     */

    /**
     * Method description
     *
     * @param datasetId
     * @param fromIndex
     * @param toIndex
     * @return
     */
    public float[][] readDataSlice(int datasetId, int fromIndex, int toIndex) {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=readDataSlice&id=" + datasetId
                    + "&from=" + fromIndex + "&to=" + toIndex);

            urlStream = IGVHttpClientUtils.openConnectionStream(url);
            DataInputStream is = new DataInputStream(new BufferedInputStream(urlStream));
            int nRows = is.readInt();
            int nColumns = is.readInt();

            float[][] dataSlice = new float[nRows][nColumns];

            for (int i = 0; i < nRows; i++) {
                for (int j = 0; j < nColumns; j++) {
                    dataSlice[i][j] = is.readFloat();
                }
            }
            is.close();

            return dataSlice;
        }
        catch (IOException ex) {
            log.error("Error in readAllFloats", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Send the server a "close" exceptionMessage for an HDF5 entity (file, dataset, or group).
     *
     * @param method
     * @param datasetId
     */
    private void closeEntity(final String method, int datasetId) {
        InputStream urlStream = null;
        try {
            URL url = new URL(locator.getServerURL() + "?method=" + method + "&id=" + datasetId);
            urlStream = IGVHttpClientUtils.openConnectionStream(url);

            // TODO -- read result;
        }
        catch (IOException ex) {
            log.error("Error closing entity", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }
    }

    /**
     * Send the server an "open" exceptionMessage for an HDF5 entity.
     *
     * @param method
     * @param parentNodeId
     * @param name
     * @return
     */
    private int openEntity(final String method, int parentNodeId, String name) {

        String key = method + "_" + name + "_" + parentNodeId;
        Integer entityId = entityCache.get(key);
        if (entityId == null) {
            entityId = openRemoteEntity(method, parentNodeId, name);
            entityCache.put(key, entityId);
        }

        return entityId.intValue();
    }

    private Integer openRemoteEntity(final String method, int parentNodeId, String name) {
        InputStream urlStream = null;
        try {

            URL url = new URL(locator.getServerURL() + "?method=" + method + "&id=" + parentNodeId
                    + "&name=" + name);
            urlStream = IGVHttpClientUtils.openConnectionStream(url);
            DataInputStream is = new DataInputStream(new BufferedInputStream(urlStream));
            int id = is.readInt();
            if (id < 0) {
                throw new ObjectNotFoundException("Entity not found: " + name);
            }
            return new Integer(id);
        }
        catch (IOException ex) {
            log.error("Error opening entity", ex);
            throw new RuntimeException(ex);
        }
        finally {
            if (urlStream != null) {
                try {
                    urlStream.close();
                }
                catch (IOException iOException) {
                    log.error("Error closing url stream", iOException);
                }
            }
        }

    }
}
