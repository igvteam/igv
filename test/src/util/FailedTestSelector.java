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

package util;

import org.apache.tools.ant.BuildException;
import org.apache.tools.ant.types.Parameter;
import org.junit.Ignore;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.*;
import java.util.HashSet;
import java.util.Set;

/**
 * Selector for re-running failed tests.
 * Examines the provided reports and only accepts files
 * which match tests that previously failed (or had an error)
 *
 * Reports must use standard JUnit XML.
 * User: jacob
 * Date: 2012-Jul-18
 */
@Ignore
public class FailedTestSelector extends org.apache.tools.ant.types.selectors.BaseExtendSelector{


    /**
     * The directory containing reports. Required input.
     */
    private String reportDir = null;

    /**
     * The filename of the test report. If
     * null, will look at every .xml file in the directory
     */
    private String reportFile = null;

    private Exception errorParsingReportFile = null;

    private Set<String> failingTests = new HashSet<String>();

    public FailedTestSelector(){
        super();
    }

    @Override
    public void setParameters(Parameter[] parameters){
        super.setParameters(parameters);
        for(Parameter p: getParameters()){
            if(p.getName().equalsIgnoreCase("reportDir")){
                reportDir = p.getValue();
            }

            if(p.getName().equalsIgnoreCase("reportFile")){
                reportFile = p.getValue();
            }
        }

        if(reportDir == null){
            throw new RuntimeException("reportDir is required");
        }


        if(reportFile != null){
            //Parse a single file, throw exception if fail
            File checkFile = new File(reportDir, reportFile);
            try {
                checkReportFile(checkFile.getAbsolutePath());
            } catch (IOException e) {
                e.printStackTrace();
                errorParsingReportFile = e;
            } catch (SAXException e) {
                e.printStackTrace();
                errorParsingReportFile = e;
            } catch (ParserConfigurationException e) {
                e.printStackTrace();
                errorParsingReportFile = e;
            }
        }else{
            File dir = new File(reportDir);
            File[] files = dir.listFiles(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return name.endsWith(".xml");
                }
            });
            for(File fi: files){
                try{
                    checkReportFile(fi.getAbsolutePath());
                } catch (SAXException e) {
                    //Don't care, might be other xml file
                } catch (ParserConfigurationException e) {
                    //
                } catch (IOException e) {
                    e.printStackTrace();
                    errorParsingReportFile = e;
                }
            }

        }
    }

    private void checkReportFile(String filePath) throws IOException, SAXException, ParserConfigurationException{
        InputStream inStream = new FileInputStream(filePath);
        Document doc = createDOMDocumentFromXmlStream(inStream);
        NodeList nodes = doc.getElementsByTagName("testsuite");
        for(int nn=0; nn < nodes.getLength(); nn++){
            Node node = nodes.item(nn);
            NamedNodeMap attribs = node.getAttributes();
            int errors = Integer.parseInt(attribs.getNamedItem("errors").getTextContent());
            int failures = Integer.parseInt(attribs.getNamedItem("failures").getTextContent());
            if(errors > 0 || failures > 0){
                //This may be a full package name, or just the class
                //All we really need is the class name
                String testName = attribs.getNamedItem("name").getTextContent();
                String[] pieces = testName.split("\\.");
                if(pieces.length > 1){
                    testName = pieces[pieces.length - 1];
                }
                failingTests.add(testName);
            }
        }
    }

    @Override
    public boolean isSelected(File directory, String filename, File file) throws BuildException {
        if(errorParsingReportFile != null){
            throw new BuildException("Error processing result file", errorParsingReportFile);
        }
        //Strip extension, only check file name (not full path)
        String baseFileName = file.getName().split("\\.", 2)[0];
        return failingTests.contains(baseFileName);
    }


    /**
     * Reads an xml from an input stream and creates DOM document.
     *
     * @param inStream
     * @return
     * @throws ParserConfigurationException
     * @throws IOException
     * @throws SAXException
     */
    public static Document createDOMDocumentFromXmlStream(InputStream inStream)
            throws ParserConfigurationException, IOException, SAXException {

        DocumentBuilder documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
        Document xmlDocument = documentBuilder.parse(inStream);

        return xmlDocument;
    }

}
