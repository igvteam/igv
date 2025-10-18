package org.broad.igv.mcp;

import org.junit.Test;

import static org.junit.Assert.*;

public class McpToolsTest {

    @Test
    public void getFileTool() {

        McpTools mcpTools = new McpTools();
        McpTools.ToolDescriptor toolDescriptor = mcpTools.getFileTool();
        assertEquals("loadFile", toolDescriptor.tool().name());
        assertNotNull(toolDescriptor.handler());
    }
}