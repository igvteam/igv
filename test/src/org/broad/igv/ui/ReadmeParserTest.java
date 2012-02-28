package org.broad.igv.ui;

import junit.framework.Assert;
import junit.framework.TestCase;
import org.broad.igv.feature.AminoAcidManager;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.net.URL;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * User: jacob
 * Date: 2012/02/27
 */
public class ReadmeParserTest{
    
    private String path = "docs/igvtools_readme.txt";
    private ReadmeParser parser;

    @Before
    public void setUp(){
        parser = new ReadmeParser(path);
    }

    @After
    public void tearDown(){
        parser = null;
    }

    @Test
    public void testGetBadCmd() throws Exception{
        String fake_cmd = "fake_cmd";
        String info = parser.getDocForCommand(fake_cmd);
        assertEquals("Command " + fake_cmd + " not found", info.split("[\\n,\\r]")[0].trim());
    }
    

    @Test
    public void testGetCommands() throws Exception{
    String[] commands = {"count", "sort", "tile", "index", "sort", "version", "toTDF", "formatexp", "gui"};
    for(String cmd: commands){
        tstGetCommand(cmd);
    }

    }
    
    public void tstGetCommand(String cmd) throws Exception{
        String info = parser.getDocForCommand(cmd);
        assertNotNull(info);
        assertFalse(info.contains("Command " + cmd + " not found"));
        assertFalse(info.matches("Command .* not found.*"));
        assertTrue(info.length() > 0);
    }
}
