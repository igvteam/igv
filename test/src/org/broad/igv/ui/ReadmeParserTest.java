package org.broad.igv.ui;

import junit.framework.TestCase;
import org.junit.Test;

/**
 * User: jacob
 * Date: 2012/02/27
 */
public class ReadmeParserTest extends TestCase {
    
    private String path = "scripts/igvtools_readme.txt";
    private ReadmeParser parser;
    
    public void setUp(){
        parser = new ReadmeParser(path);
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
        assertFalse(info.matches("Command .* not found.*"));
        assertTrue(info.length() > 0);
    }
}
