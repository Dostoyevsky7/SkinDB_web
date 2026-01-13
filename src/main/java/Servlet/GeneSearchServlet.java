package Servlet;

import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.*;
import java.io.*;
import java.util.*;
import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;
import org.apache.commons.csv.*;

@WebServlet("/gene-search")
public class GeneSearchServlet extends HttpServlet {

    private static final Map<String, Map<String, String>> mapping;

    static {
        try (InputStreamReader reader = new InputStreamReader(
                GeneSearchServlet.class.getClassLoader().getResourceAsStream("mapping.csv.json"),
                "UTF-8"
        )) {
            mapping = new Gson().fromJson(
                    reader,
                    new TypeToken<Map<String, Map<String, String>>>() {}.getType()
            );
        } catch (Exception e) {
            throw new ExceptionInInitializerError("Failed to load mapping.csv.json: " + e.getMessage());
        }
    }

    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {

        String query = request.getParameter("q");
        if (query == null || query.trim().isEmpty()) {
            response.sendError(HttpServletResponse.SC_BAD_REQUEST, "Missing query parameter: q");
            return;
        }

        query = query.trim().toLowerCase();
        int maxResults = 500; // Limit results to prevent overwhelming response

        List<Map<String, String>> results = new ArrayList<>();

        // Search across all datasets
        for (Map.Entry<String, Map<String, String>> entry : mapping.entrySet()) {
            if (results.size() >= maxResults) break;

            String said = entry.getKey();
            Map<String, String> meta = entry.getValue();
            String csvRelPath = meta.get("csv_path");

            if (csvRelPath == null) continue;

            String realPath = getServletContext().getRealPath("/" + csvRelPath);
            File csvFile = new File(realPath);

            if (!csvFile.exists()) continue;

            try (
                Reader in = new FileReader(csvFile);
                CSVParser parser = CSVFormat.DEFAULT
                        .withFirstRecordAsHeader()
                        .withIgnoreSurroundingSpaces()
                        .withTrim()
                        .parse(in)
            ) {
                // Build column name mapping (case-insensitive)
                Map<String, String> lowerToOrig = new HashMap<>();
                for (String orig : parser.getHeaderMap().keySet()) {
                    lowerToOrig.put(orig.trim().toLowerCase(), orig);
                }

                // Check required columns exist
                if (!lowerToOrig.containsKey("names") ||
                    !lowerToOrig.containsKey("logfoldchanges") ||
                    !lowerToOrig.containsKey("pvals_adj")) {
                    continue;
                }

                String namesCol = lowerToOrig.get("names");
                String fcCol = lowerToOrig.get("logfoldchanges");
                String pvalCol = lowerToOrig.get("pvals_adj");
                String groupCol = lowerToOrig.get("group");

                for (CSVRecord rec : parser) {
                    if (results.size() >= maxResults) break;

                    try {
                        String geneName = rec.get(namesCol);

                        // Case-insensitive partial match
                        if (geneName.toLowerCase().contains(query)) {
                            double logfc = Double.parseDouble(rec.get(fcCol));
                            double pval = Double.parseDouble(rec.get(pvalCol));
                            String group = groupCol != null ? rec.get(groupCol) : "Unknown";

                            Map<String, String> row = new HashMap<>();
                            row.put("gene", geneName);
                            row.put("said", said);
                            row.put("gse", meta.get("GSE"));
                            row.put("gsm", meta.get("GSM"));
                            row.put("logfc", String.format("%.4f", logfc));
                            row.put("pval", String.format("%.2e", pval));
                            row.put("group", group);
                            results.add(row);
                        }
                    } catch (Exception ignore) {
                        // Skip malformed rows
                    }
                }
            } catch (IOException e) {
                // Skip files that can't be read
            }
        }

        // Sort results by p-value (ascending)
        results.sort((a, b) -> {
            try {
                double pvalA = Double.parseDouble(a.get("pval"));
                double pvalB = Double.parseDouble(b.get("pval"));
                return Double.compare(pvalA, pvalB);
            } catch (Exception e) {
                return 0;
            }
        });

        response.setContentType("application/json;charset=UTF-8");
        Map<String, Object> responseData = new HashMap<>();
        responseData.put("query", query);
        responseData.put("count", results.size());
        responseData.put("results", results);
        new Gson().toJson(responseData, response.getWriter());
    }
}
