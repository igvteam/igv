package org.broad.igv.mcp;

import io.modelcontextprotocol.server.McpSyncServerExchange;
import io.modelcontextprotocol.spec.McpSchema;
import org.broad.igv.batch.CommandExecutor;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;

public class McpTools {

    private static Logger log = LogManager.getLogger(McpTools.class);

    static IGV igv = IGV.hasInstance() ? IGV.getInstance() : null;   // Null case for unit tests

    static CommandExecutor executor = new CommandExecutor(igv);

    // Existing convenience builder (deprecated in 0.14.x but retained for compatibility)
    public ToolDescriptor getFileTool() {
        //return new ToolDescriptor(loadFileTool(), loadFileHandler());
        return createTool(
                "loadFile",
                "Loads a file from a given path or URL",
                List.of(new ToolArgument("path", "A file path or URL"))
        );
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

    static ToolDescriptor createTool(String name, String description, List<ToolArgument> arguments) {

        Map<String, Object> properties = new HashMap<>();

        for (ToolArgument arg : arguments) {
            McpSchema.JsonSchema argSchema = new McpSchema.JsonSchema(
                    "string",
                    Map.of("description", arg.description()),
                    null,
                    null,
                    null,
                    null);
            properties.put(arg.name(), argSchema);
        }

        McpSchema.JsonSchema inputSchema = new McpSchema.JsonSchema(
                "object",
                properties,
                arguments.stream().map(ToolArgument::name).toList(),
                null,
                null,
                null);

        McpSchema.Tool tool = McpSchema.Tool.builder()
                .name(name)
                .description(description)
                .inputSchema(inputSchema)
                .build();

        BiFunction<McpSyncServerExchange, Map<String, Object>, McpSchema.CallToolResult> handler =

                (exchange, args) -> {

                    String argString = args.values().stream().map(o -> o.toString()).collect(java.util.stream.Collectors.joining(" "));
                    log.info("Executing tool: " + name + " with args: " + argString);
                    String result = executor.execute(name + " " + argString);

                    return McpSchema.CallToolResult.builder()
                            .content(List.of(new McpSchema.TextContent(result)))
                            .isError(false)
                            .build();
                };

        return new ToolDescriptor(tool, handler);
    }

    record ToolDescriptor(McpSchema.Tool tool,
                          BiFunction<McpSyncServerExchange, Map<String, Object>, McpSchema.CallToolResult> handler) {
    }

    record ToolArgument(String name, String description) {
    }
}

/*
	public record JsonSchema( // @formatter:off
		@JsonProperty("type") String type,
		@JsonProperty("properties") Map<String, Object> properties,
		@JsonProperty("required") List<String> required,
		@JsonProperty("additionalProperties") Boolean additionalProperties,
		@JsonProperty("$defs") Map<String, Object> defs,
		@JsonProperty("definitions") Map<String, Object> definitions) { // @formatter:on
	}

 */