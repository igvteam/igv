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

package org.broad.igv.gs.atm;

import org.broad.igv.Globals;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.StringUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.lang.reflect.Field;
import java.net.Authenticator;
import java.net.PasswordAuthentication;
import java.util.*;

import static junit.framework.Assert.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: 8/8/11
 * Time: 9:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class ATMUtilsTest {

    @Before
    public void setUp() {
        Globals.setTesting(true);
        HttpUtils.getInstance().setAuthenticator(new GSTestAuthenticator());
    }

    @After
    public void tearDown() {
        HttpUtils.getInstance().resetAuthenticator();
    }

    @Test
    public void testGetWebTools() throws Exception {
        List<WebToolDescriptor> webTools = ATMUtils.getWebTools();

        // Build a map of name -> tool
        Map<String, WebToolDescriptor> toolMap = new HashMap();
        for (WebToolDescriptor wt : webTools) {
            toolMap.put(wt.getName(), wt);
        }
        String toolName = "IGV";
        WebToolDescriptor igvDesc = toolMap.get(toolName);
        assertEquals(toolName, igvDesc.getName());
        checkWebToolDescriptor(igvDesc);
    }

    @Test
    public void testGetWebTool() throws Exception {
        String toolname = "IGV";
        WebToolDescriptor igvDesc = ATMUtils.getWebTool(toolname);
        assertEquals(toolname, igvDesc.getName());
        checkWebToolDescriptor(igvDesc);
    }

    @Test
    public void testGetIGVLaunchURL() throws Exception {

        String file = "/users/igvtest/Broad.080528.subtypes.seg.gz";
        String toolname = "IGV";
        WebToolDescriptor igvDesc = ATMUtils.getWebTool(toolname);
        String url = StringUtils.decodeURL(ATMUtils.getWebtoolLaunchURL(igvDesc, file));
        assertTrue(url.startsWith(igvDesc.getBaseUrl() + "?sessionURL="));
        assertTrue(url.endsWith(file));
    }

    private void checkWebToolDescriptor(WebToolDescriptor desc) throws Exception{
//        assertNotNull("Descriptor is null", desc);
//        assertNotNull("Author is null", desc.getAuthor());
//        assertNotNull("BasUrl is null", desc.getBaseUrl());
//        assertNotNull("Description is null", desc.getDescription());
//        assertNotNull("Id is null", desc.getId());
//        assertNotNull("name is null",desc.getName());
        checkFieldsNotNull(desc);
        for(FileParameter param: desc.getFileParameters()){
            checkFieldsNotNull(param);
            for(GSDataFormat format: param.getFormats()){
                checkFieldsNotNull(format);
            }
        }
    }

    @Test
    public void testGetUCSCLaunchURL() throws Exception {

        List<WebToolDescriptor> webTools = ATMUtils.getWebTools();

        // Build a map of name -> tool
        Map<String, WebToolDescriptor> toolMap = new HashMap();
        for (WebToolDescriptor wt : webTools) {
            toolMap.put(wt.getName(), wt);
        }

        String toolname = ("UCSC Genome Browser");
        WebToolDescriptor ucscDesc = ATMUtils.getWebTool(toolname);
        checkWebToolDescriptor(ucscDesc);

        String url = StringUtils.decodeURL(ATMUtils.getWebtoolLaunchURL(ucscDesc));
        assertNotNull(url);
    }

    private void checkFieldsNotNull(Object obj) throws Exception{
        checkFieldsNotNull(obj, Collections.<String>emptySet());
    }

    /**
     * Check that all fields in obj are not null
     * @param obj
     * @param excludedFieldNames Names of fields we don't check
     */
    private void checkFieldsNotNull(Object obj, Set<String> excludedFieldNames) throws Exception{
        assertNotNull(obj);
        Field[] fields = obj.getClass().getFields();
        for(Field field: fields){
            field.setAccessible(true);
            String fieldName = field.getName();
            if(excludedFieldNames.contains(fieldName)){
                continue;
            }

            Object value = field.get(obj);
            assertNotNull(field.getName() + " is null", value);
        }

    }

    static class GSTestAuthenticator extends Authenticator {

        @Override
        protected PasswordAuthentication getPasswordAuthentication() {
            return new PasswordAuthentication("igvtest", "igvtest".toCharArray());
        }
    }


}
