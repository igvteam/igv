package org.broad.igv.mcp;

import com.fasterxml.jackson.databind.ObjectMapper;
import io.modelcontextprotocol.json.jackson.JacksonMcpJsonMapper;
import io.modelcontextprotocol.server.McpServer;
import io.modelcontextprotocol.server.transport.StdioServerTransportProvider;
import org.broad.igv.logging.LogManager;
import org.broad.igv.logging.Logger;


public class IGVMcpServer {

    static Logger log = LogManager.getLogger(IGVMcpServer.class);

    public static void main(String[] args) {
        start();
    }

    public static void start() {

        log.info("Starting IGV MCP Server");

        StdioServerTransportProvider transportProvider = new StdioServerTransportProvider(
                new JacksonMcpJsonMapper(new ObjectMapper())
        );

        var tools = McpTools.getTools();

        var serverSpecification = McpServer.sync(transportProvider)
                .serverInfo("igv-mcp", "1.0.0");

        for(McpTools.ToolDescriptor td : tools) {
            serverSpecification = serverSpecification.tool(td.tool(), td.handler());
        }

        var server = serverSpecification.build();

        Runtime.getRuntime().addShutdownHook(new Thread(() -> {
            server.close();    // Ensure server is closed on JVM shutdown.  Probably not strictly necessary since stdio will close.
            log.info("IGV MCP Server has shut down.");
        }));
    }
}
