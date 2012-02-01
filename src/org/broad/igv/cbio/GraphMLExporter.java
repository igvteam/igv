/*
 * Copyright (c) 2007-2011 by The Broad Institute of MIT and Harvard.All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL),
 * Version 2.1 which is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR
 * WARRANTIES OF ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING,
 * WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER
 * OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR RESPECTIVE
 * TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES
 * OF ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES,
 * ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER
 * THE BROAD OR MIT SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT
 * SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 */

package org.broad.igv.cbio;

import org.jgrapht.Graph;
import org.jgrapht.ext.EdgeNameProvider;
import org.jgrapht.ext.VertexNameProvider;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.Map;

/**
 * Exports an org.jgrapht.graph into graphml format.
 * The fields on each graph and node object are written (using reflection)
 * to data elements.
 * <p/>
 * User: jacob
 * Date: 2012/02/01
 */
public class GraphMLExporter<V, E extends BaseEdge> {

    private Class<V> vertexClass;
    private Class<E> edgeClass;

    /*
     * Provider which generates a unique ID (not name, ID) for each vertex.
     */
    private VertexNameProvider<V> vIDprovider = new VertexNameProvider<V>() {
        public String getVertexName(V v) {
            return "" + v.hashCode();
        }
    };

    /**
     * Provider which generates a unique ID (not name, ID) for each edge.
     */
    private EdgeNameProvider<E> eIDprovider = new EdgeNameProvider<E>() {
        public String getEdgeName(E e) {
            return "" + e.hashCode();
        }
    };


    public GraphMLExporter(Class<V> vertexClass, Class<E> edgeClass) {
        this.vertexClass = vertexClass;
        this.edgeClass = edgeClass;
    }

    /**
     * @param vertexClass
     * @param edgeClass
     * @param vIDprovider Provider used to generate IDs for each vertex. IDs must
     *                    be unique. Default implementation is hashCode of v
     * @param eIDprovider Provider used to generate IDs for each edge. IDs must
     *                    be unique. Default implementation is hashCode of v
     */
    public GraphMLExporter(Class<V> vertexClass, Class<E> edgeClass,
                           VertexNameProvider<V> vIDprovider,
                           EdgeNameProvider<E> eIDprovider) {
        this(vertexClass, edgeClass);
        this.vIDprovider = vIDprovider;
        this.eIDprovider = eIDprovider;
    }

    public void exportGraph(String outputFile, Graph<V, E> graph) {
        try {
            // Create a DOM document
            DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
            Document document = documentBuilder.newDocument();
            document.setStrictErrorChecking(false);

            // Global root element
            Element globalElement = document.createElement("graphml");

            Field[] edgeFields = edgeClass.getFields();
            addSchema(document, globalElement, edgeFields, "edge");

            Element graphEl = document.createElement("graph");
            String edgedefault = getEdgeDefault(graph);

            graphEl.setAttribute("edgedefault", edgedefault);

            //Generate nodes
            Element docNode;
            for (V v : graph.vertexSet()) {
                docNode = document.createElement("node");
                docNode.setAttribute("id", vIDprovider.getVertexName(v));
                addData(document, vertexClass.getFields(), v, docNode);
                graphEl.appendChild(docNode);
            }

            //Generate edges
            Element docEdge;
            for (E e : graph.edgeSet()) {
                docEdge = document.createElement("edge");
                docEdge.setAttribute("source", vIDprovider.getVertexName(graph.getEdgeSource(e)));
                docEdge.setAttribute("target", vIDprovider.getVertexName(graph.getEdgeTarget(e)));
                docEdge.setAttribute("directed", "" + e.isDirected());
                addData(document, e, docEdge);
                graphEl.appendChild(docEdge);
            }

            globalElement.appendChild(graphEl);
            document.appendChild(globalElement);

            writeDocument(document, outputFile);
        } catch (Exception e) {
            throw new RuntimeException("Error outputting graph");
        }
    }

    private String getEdgeDefault(Graph<V, E> graph) {
        String edgedefault = "undirected";
        Class[] interfaces = graph.getClass().getInterfaces();
        for (Class inter : interfaces) {
            //If interface name contains directed but not undirected
            //consider it directed
            if (inter.getSimpleName().toLowerCase().matches(".*(?<!un)directed.*")) {
                edgedefault = "directed";
                break;
            }
        }
        return edgedefault;
    }

    /**
     * Write document to XML at outputFile. File is deleted if there
     * is an error writing out.
     *
     * @param document
     * @param outputFile
     * @return success
     * @throws java.io.IOException
     */
    private boolean writeDocument(Document document, String outputFile) throws IOException {
        StreamResult streamResult;
        boolean success = false;
        try {

            TransformerFactory factory = TransformerFactory.newInstance();
            Transformer transformer = factory.newTransformer();

            transformer.setOutputProperty(OutputKeys.INDENT, "yes");
            transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "4");

            streamResult = new StreamResult(new StringWriter());
            DOMSource source = new DOMSource(document);
            transformer.transform(source, streamResult);
        } catch (TransformerException e) {
            return false;
        }

        String xmlString = streamResult.getWriter().toString();
        //System.out.println(xmlString);

        FileWriter fileWriter = new FileWriter(outputFile);
        try {
            fileWriter.write(xmlString);
            success = true;
        } catch (IOException e) {
            File outfile = new File(outputFile);
            if (outfile.exists()) {
                outfile.delete();
            }
        } finally {
            if (fileWriter != null) {
                fileWriter.close();
            }
        }
        return success;
    }

    /**
     * Adds schema information (name, type) from field list
     *
     * @param document
     * @param topElement
     * @param fields     Fields to add to schema. Only primitive types will be included
     * @param typeFor    Which element this schema is for. Legal values: "graph", "all", "node", "edge", null (will not
     *                   be applied.
     */
    private void addSchema(Document document, Element topElement, Field[] fields, String typeFor) {
        Element key;
        for (Field field : fields) {
            if (!includeField(field)) {
                continue;
            }
            key = document.createElement("key");
            key.setAttribute("id", field.getName());
            key.setAttribute("type", field.getType().getSimpleName().toLowerCase());

            if (typeFor != null) {
                key.setAttribute("for", typeFor);
            }
            topElement.appendChild(key);
        }
    }

    /**
     * Add data from source to DOM Element dest.
     *
     * @param document
     * @param sourceData
     * @param dest
     */
    private void addData(Document document, Map<String, GraphMLData> sourceData, Element dest) {
        Element data;
        for (String s : sourceData.keySet()) {
            data = document.createElement("data");
            data.setAttribute(s, sourceData.get(s).getValue());
            dest.appendChild(data);
        }
    }

    /**
     * Add data from source to DOM Element dest. Each field is retrieved from the source object iff
     * it's primitive.
     *
     * @param document
     * @param fields
     * @param source
     * @param dest
     * @throws IllegalAccessException
     */
    private void addData(Document document, Field[] fields, Object source, Element dest) throws IllegalAccessException {
        Element data;
        for (Field field : fields) {
            if (!includeField(field)) {
                continue;
            }
            data = document.createElement("data");
            data.setAttribute(field.getName(), "" + field.get(source));
            dest.appendChild(data);
        }
    }

    /**
     * Whether field should be included in output. By default we
     * take primitives + string, non-static, and accessible
     *
     * @param field
     * @return
     */
    protected boolean includeField(Field field) {
        return (field.getType().isPrimitive() || field.getType().equals(String.class)) &&
                !Modifier.isStatic(field.getModifiers());
    }


}
