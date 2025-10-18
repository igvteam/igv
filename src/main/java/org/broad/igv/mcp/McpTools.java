package org.broad.igv.mcp;

import io.modelcontextprotocol.server.McpServerFeatures;
import io.modelcontextprotocol.server.McpSyncServerExchange;
import io.modelcontextprotocol.spec.McpSchema;
import org.broad.igv.batch.CommandExecutor;
import org.broad.igv.ui.IGV;

import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;

public class McpTools {

    static CommandExecutor executor = new CommandExecutor(IGV.getInstance());

    // Existing convenience builder (deprecated in 0.14.x but retained for compatibility)
    public McpServerFeatures.SyncToolSpecification loadFile() {
        McpSchema.Tool tool = loadFileTool();
        return new McpServerFeatures.SyncToolSpecification(tool, loadFileHandler());
    }

    // New: return the tool schema/spec without the handler
    public McpSchema.Tool loadFileTool() {
        McpSchema.JsonSchema pathSchema = new McpSchema.JsonSchema(
                "string",
                Map.of("description", "A file path or URL"),
                null,
                null,
                null,
                null);

        McpSchema.JsonSchema inputSchema = new McpSchema.JsonSchema(
                "object",
                Map.of("path", pathSchema),
                List.of("path"),
                null,
                null,
                null);

        return McpSchema.Tool.builder()
                .name("loadFile")
                .description("Loads a file from a given path or URL")
                .inputSchema(inputSchema)
                .build();
    }

    // New: return the handler for the tool
    public BiFunction<McpSyncServerExchange, Map<String, Object>, McpSchema.CallToolResult> loadFileHandler() {
        return (exchange, args) -> {
            String path = (String) args.get("path");
            String result = executor.execute("load " + path);

            return McpSchema.CallToolResult.builder()
                    .content(List.of(new McpSchema.TextContent(result)))
                    .isError(false)
                    .build();
        };
    }
}