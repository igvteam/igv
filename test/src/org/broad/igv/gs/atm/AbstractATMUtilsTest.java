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

import org.broad.igv.AbstractHeadlessTest;
import org.broad.igv.gs.GSUtils;
import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.StringUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

import java.lang.reflect.Field;
import java.util.*;

import static junit.framework.Assert.*;

/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: 8/8/11
 * Time: 9:16 PM
 */
@Ignore
public abstract class AbstractATMUtilsTest extends AbstractHeadlessTest{

    @Before
    public void setUp() throws Exception{
        super.setUp();
        //This is pretty dumb. The reason is that initializing HttpUtils sets the authenticator,
        //and we need to overwrite it in initAuth. It's not actually important which method we call,
        //as long as HttpUtils is initialized so it doesn't get initialized later
        HttpUtils.getInstance().resetAuthenticator();
        GSUtils.logout();
        initAuth();
    }

    @After
    public void tearDown() {
        GSUtils.logout();
        HttpUtils.getInstance().resetAuthenticator();
    }

    protected abstract void initAuth();

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

}
