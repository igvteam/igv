package org.broad.igv.mcp;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import io.modelcontextprotocol.server.McpSyncServerExchange;
import io.modelcontextprotocol.spec.McpSchema;
import org.broad.igv.batch.CommandExecutor;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;
import org.broad.igv.ui.IGV;

import java.io.IOException;
import java.io.InputStream;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

public class McpTools {

    private static final Logger log = LogManager.getLogger(McpTools.class);

    static IGV igv = IGV.hasInstance() ? IGV.getInstance() : null;   // Null case for unit tests

    static CommandExecutor executor = new CommandExecutor(igv);

    // This record is used for deserializing from JSON
    record ToolDefinition(String name, String description, List<ToolArgument> arguments) {
    }

    // Existing convenience builder (deprecated in 0.14.x but retained for compatibility)
    public static ToolDescriptor getFileTool() {
        //return new ToolDescriptor(loadFileTool(), loadFileHandler());
        return createTool(
                "loadFile",
                "Loads a file from a given path or URL",
                List.of(new ToolArgument("path", "A file path or URL", null, false))
        );
    }

    public static List<ToolDescriptor> getTools() {

        try (InputStream is = McpTools.class.getResourceAsStream("/tools.yaml")) {
            if (is == null) {
                log.error("Could not find tools.yaml");
                return Collections.emptyList();
            }
            ObjectMapper mapper = new ObjectMapper(new YAMLFactory());
            List<ToolDefinition> toolDefs = mapper.readValue(is, new TypeReference<>() {
            });

            return toolDefs.stream()
                    .map(def -> createTool(
                            def.name(),
                            def.description(),
                            def.arguments()))
                    .collect(Collectors.toList());

        } catch (IOException e) {
            log.error("Error reading tools.yaml", e);
            return Collections.emptyList();
        }
    }

    static ToolDescriptor createTool(String name, String description, List<ToolArgument> arguments) {

        Map<String, Object> properties = new LinkedHashMap<>();
        if (arguments != null) {
            for (ToolArgument arg : arguments) {
                Map<String, Object> argSchemaMap = new LinkedHashMap<>();
                argSchemaMap.put("type", "string");
                argSchemaMap.put("description", arg.description());

                if (arg.enumValues() != null && !arg.enumValues().isEmpty()) {
                    List<Map<String, Object>> oneOfList = arg.enumValues().stream()
                            .map(ev -> {
                                Map<String, Object> constSchema = new LinkedHashMap<>();
                                constSchema.put("const", ev.value());
                                constSchema.put("description", ev.description());
                                return constSchema;
                            })
                            .collect(Collectors.toList());
                    argSchemaMap.put("oneOf", oneOfList);
                }
                properties.put(arg.name(), argSchemaMap);
            }
        }

        List<String> requiredArgs = new java.util.ArrayList<>();
        if (arguments != null) {
            requiredArgs = arguments.stream()
                    .filter(arg -> !arg.optional())
                    .map(ToolArgument::name)
                    .toList();
        }

        McpSchema.JsonSchema inputSchema = new McpSchema.JsonSchema(
                "object",
                properties,
                requiredArgs,
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

                    String argString = "";
                    if (args != null && !args.isEmpty()) {
                        if (arguments != null) {
                            argString = arguments.stream()
                                    .map(ToolArgument::name)
                                    .filter(args::containsKey)
                                    .map(args::get)
                                    .map(Object::toString)
                                    .collect(Collectors.joining(" "));
                        } else {
                            // Fallback for safety, though arguments should not be null if args are present
                            argString = args.values().stream().map(Object::toString).collect(java.util.stream.Collectors.joining(" "));
                        }
                    }
                    log.info("Executing tool: " + name + " with args: " + argString);
                    String command = (name + " " + argString).trim();
                    try {
                        String result = executor.execute(command);
                        boolean isError = !"OK".equalsIgnoreCase(result);
                        return McpSchema.CallToolResult.builder()
                                .content(List.of(new McpSchema.TextContent(result)))
                                .isError(isError)
                                .build();
                    } catch (Exception e) {
                        String errorMessage = "Error executing tool '" + name + "': " + e.getMessage();
                        log.error(errorMessage, e);
                        throw new RuntimeException(errorMessage, e);
                    }
                };

        return new ToolDescriptor(tool, handler);
    }

    public record ToolDescriptor(McpSchema.Tool tool,
                                 BiFunction<McpSyncServerExchange, Map<String, Object>, McpSchema.CallToolResult> handler) {
    }

    record EnumValue(String value, String description) {
    }

    record ToolArgument(String name, String description, List<EnumValue> enumValues, boolean optional) {
    }

}
