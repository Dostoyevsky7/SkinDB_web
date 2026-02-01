package Servlet;

import com.google.gson.Gson;

import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Servlet for managing the Interactive Visualization Suite
 * Handles starting/stopping the Python Dash visualization server
 */
@WebServlet("/visualization")
public class VisualizationServlet extends HttpServlet {

    // Store running visualization processes
    private static final Map<String, Process> runningProcesses = new ConcurrentHashMap<>();
    private static final int VIZ_PORT_BASE = 8050;

    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {

        String datasetId = request.getParameter("dataset");

        if (datasetId == null || datasetId.isEmpty()) {
            response.sendError(HttpServletResponse.SC_BAD_REQUEST, "Dataset ID is required");
            return;
        }

        try {
            // Get dataset path from mapping.json
            String datasetPath = getDatasetPath(datasetId);

            if (datasetPath == null || !Files.exists(Paths.get(datasetPath))) {
                response.sendError(HttpServletResponse.SC_NOT_FOUND, "Dataset file not found");
                return;
            }

            // Check if visualization server is already running for this dataset
            if (!runningProcesses.containsKey(datasetId)) {
                startVisualizationServer(datasetId, datasetPath);
            }

            // Set attributes for JSP
            request.setAttribute("datasetId", datasetId);
            request.setAttribute("vizPort", VIZ_PORT_BASE);
            String scheme = request.getScheme();
            String host = request.getServerName();
            request.setAttribute("vizUrl", scheme + "://" + host + ":" + VIZ_PORT_BASE + "/viz/");

            // Forward to JSP
            request.getRequestDispatcher("visualization.jsp").forward(request, response);

        } catch (Exception e) {
            e.printStackTrace();
            response.sendError(HttpServletResponse.SC_INTERNAL_SERVER_ERROR,
                             "Error starting visualization: " + e.getMessage());
        }
    }

    /**
     * API endpoint to check visualization server status
     */
    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {

        String action = request.getParameter("action");

        if ("status".equals(action)) {
            String datasetId = request.getParameter("dataset");
            boolean isRunning = runningProcesses.containsKey(datasetId);

            Map<String, Object> result = new HashMap<>();
            result.put("running", isRunning);
            result.put("port", VIZ_PORT_BASE);

            Gson gson = new Gson();
            response.setContentType("application/json");
            response.getWriter().write(gson.toJson(result));

        } else if ("stop".equals(action)) {
            String datasetId = request.getParameter("dataset");
            stopVisualizationServer(datasetId);

            Map<String, Object> result = new HashMap<>();
            result.put("stopped", true);

            Gson gson = new Gson();
            response.setContentType("application/json");
            response.getWriter().write(gson.toJson(result));
        }
    }

    /**
     * Start the Python Dash visualization server
     */
    private void startVisualizationServer(String datasetId, String datasetPath) throws IOException {
        String pythonScript = getServletContext().getRealPath("/WEB-INF/visualization_suite.py");

        // Check if Python script exists
        if (!Files.exists(Paths.get(pythonScript))) {
            throw new IOException("Visualization script not found: " + pythonScript);
        }

        // Build command
        ProcessBuilder pb = new ProcessBuilder(
            "python3",
            pythonScript,
            "--dataset", datasetPath,
            "--port", String.valueOf(VIZ_PORT_BASE)
        );

        pb.redirectErrorStream(true);

        // Start process
        Process process = pb.start();

        // Store process
        runningProcesses.put(datasetId, process);

        // Monitor output in separate thread
        new Thread(() -> {
            try (BufferedReader reader = new BufferedReader(
                    new InputStreamReader(process.getInputStream()))) {
                String line;
                while ((line = reader.readLine()) != null) {
                    System.out.println("[Viz " + datasetId + "] " + line);
                }
            } catch (IOException e) {
                System.err.println("Error reading visualization output: " + e.getMessage());
            }
        }).start();

        System.out.println("Started visualization server for dataset: " + datasetId);

        // Wait a moment for server to start
        try {
            Thread.sleep(3000);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
    }

    /**
     * Stop visualization server for a dataset
     */
    private void stopVisualizationServer(String datasetId) {
        Process process = runningProcesses.remove(datasetId);
        if (process != null && process.isAlive()) {
            process.destroy();
            try {
                process.waitFor();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
            System.out.println("Stopped visualization server for dataset: " + datasetId);
        }
    }

    /**
     * Get dataset file path from mapping.json
     */
    private String getDatasetPath(String datasetId) {
        try {
            String mappingPath = getServletContext().getRealPath("/WEB-INF/classes/mapping.json");
            if (!Files.exists(Paths.get(mappingPath))) {
                mappingPath = getServletContext().getRealPath("/WEB-INF/classes/") + "../resources/mapping.json";
            }

            if (Files.exists(Paths.get(mappingPath))) {
                String content = new String(Files.readAllBytes(Paths.get(mappingPath)));
                Gson gson = new Gson();
                @SuppressWarnings("unchecked")
                Map<String, Map<String, String>> mapping = gson.fromJson(content, Map.class);

                if (mapping.containsKey(datasetId)) {
                    return mapping.get(datasetId).get("file_path");
                }
            }
        } catch (Exception e) {
            System.err.println("Error loading dataset mapping: " + e.getMessage());
        }

        return null;
    }

    /**
     * Clean up resources on servlet destruction
     */
    @Override
    public void destroy() {
        System.out.println("Stopping all visualization servers...");
        for (Map.Entry<String, Process> entry : runningProcesses.entrySet()) {
            Process process = entry.getValue();
            if (process.isAlive()) {
                process.destroy();
            }
        }
        runningProcesses.clear();
        super.destroy();
    }
}
