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

package org.broad.igv.session;

import org.broad.igv.util.Utilities;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Represents a set of attributes, as well as their children.
 * Designed to model XML Node or JSONObject tags, where
 * attributes are applied at one level and elements
 * exist as children
 *
 * Contains get, put, and putAll methods, which act on the
 * attribute map only (don't affect children)
 * User: jacob
 * Date: 2012-Dec-17
 */
public class RecursiveAttributes {

    private final String name;
    private final Map<String, String> attributes;
    private List<RecursiveAttributes> children;

    public RecursiveAttributes(String name, Map<String, String> attributes){
        this(name, attributes, null);
    }

    /**
     *
     * @param name Name of the tag; needed if used as a child.
     * @param attributes
     * @param children
     */
    public RecursiveAttributes(String name, Map<String, String> attributes, List<RecursiveAttributes> children){
        this.attributes = attributes != null ? attributes : new HashMap<String, String>();
        this.children = children != null ? children : new ArrayList<RecursiveAttributes>();
        this.name = name;
    }

    /**
     * The name of this element. May be null, because it is
     * usually un-necessary at top level.
     * @return
     */
    public String getName(){
        return this.name;
    }

    /**
     * A map of attributes at this level
     * @return
     */
    public Map<String, String> getAttributes(){
        return this.attributes;
    }

    public void putAll(Map<String, String> persistentState) {
        this.attributes.putAll(persistentState);
    }

    public String put(String key, String value){
        return this.attributes.put(key, value);
    }

    public String get(String key){
        return this.attributes.get(key);
    }

    public List<RecursiveAttributes> getChildren(){
        return this.children;
    }


    /**
     * Writes the attributes from {@code recursiveAttributes.getAttributes()} into {@code parentElement}, and
     * writes children from {@code recursiveAttributes.getChildren()} as child elements.
     * @param document
     * @param parentElement
     * @param recursiveAttributes
     */
    public static void writeElement(Element parentElement, Document document, RecursiveAttributes recursiveAttributes){

        for(Map.Entry<String, String> entry: recursiveAttributes.getAttributes().entrySet()){
            parentElement.setAttribute(entry.getKey(), entry.getValue());
        }

        List<RecursiveAttributes> children = recursiveAttributes.getChildren();
        if(children != null){
            for(RecursiveAttributes childAttributes: recursiveAttributes.getChildren()){
                Element childElement = document.createElement(childAttributes.getName());
                parentElement.appendChild(childElement);
                writeElement(childElement, document, childAttributes);
            }
        }
    }


    /**
     * Writes the attributes from {@code recursiveAttributes.getAttributes()} into {@code parentElement}, and
     * writes children from {@code recursiveAttributes.getChildren()} as child elements.
     * @param parentElement
     * @return recursiveAttributes
     */
    public static RecursiveAttributes readElement(Node parentElement){

        String name = parentElement.getNodeName();
        Map<String, String> parentMap = Utilities.getAttributes(parentElement);

        List<RecursiveAttributes> childAttributeList = null;
        if (parentElement.hasChildNodes()) {
            NodeList childNodes = parentElement.getChildNodes();

            childAttributeList = new ArrayList<RecursiveAttributes>(childNodes.getLength());
            for(int ci=0; ci < childNodes.getLength(); ci++ ){
                Node childNode = childNodes.item(ci);
                RecursiveAttributes childAttributes = readElement(childNode);

                //May wish to reconsider this, but it's basically to ignore blank text nodes
                //which show up all the time
                if(childAttributes.getAttributes().size() == 0 && childAttributes.getChildren().size() == 0){
                    continue;
                }

                childAttributeList.add(childAttributes);
            }
        }

        return new RecursiveAttributes(name, parentMap, childAttributeList);
    }

}
