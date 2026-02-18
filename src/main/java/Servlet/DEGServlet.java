package Servlet;

import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.*;
import java.io.*;
import java.util.*;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;
import org.apache.commons.csv.*;

@WebServlet("/deg")
public class DEGServlet extends HttpServlet {

    private static final Map<String, Map<String, String>> mapping;

    static {
        try (InputStreamReader reader = new InputStreamReader(
                DEGServlet.class.getClassLoader().getResourceAsStream("mapping.csv.json"),
                "UTF-8"
        )) {
            mapping = new Gson().fromJson(
                    reader,
                    new TypeToken<Map<String, Map<String, String>>>() {}.getType()
            );
        } catch (Exception e) {
            throw new ExceptionInInitializerError("加载 mapping.csv.json 失败: " + e.getMessage());
        }
    }

    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {

        // 1) 获取参数
        String said = request.getParameter("said");
        if (said == null || said.isEmpty()) {
            response.sendError(HttpServletResponse.SC_BAD_REQUEST, "缺少参数: said");
            return;
        }

        Map<String, String> meta = mapping.get(said);
        if (meta == null) {
            response.sendError(HttpServletResponse.SC_NOT_FOUND, "未找到 SAID 映射: " + said);
            return;
        }

        double pvalThreshold = parseDouble(request.getParameter("pval"), 0.05);
        double fcThreshold   = parseDouble(request.getParameter("fc"), 1.0);
        String filterGroup   = request.getParameter("group"); // 可为空

        String csvRelPath = meta.get("csv_path");
        String realPath   = getServletContext().getRealPath("/" + csvRelPath);
        File csvFile      = new File(realPath);
        if (!csvFile.exists()) {
            response.sendError(HttpServletResponse.SC_NOT_FOUND, "CSV 文件不存在: " + realPath);
            return;
        }

        List<Map<String,String>> outList = new ArrayList<>();

        try (
                Reader in = new FileReader(csvFile);
                CSVParser parser = CSVFormat.DEFAULT
                        .withFirstRecordAsHeader()
                        .withIgnoreSurroundingSpaces()
                        .withTrim()
                        .parse(in)
        ) {
            // 建立列名映射（忽略大小写）
            Map<String,String> lowerToOrig = new HashMap<>();
            for (String orig : parser.getHeaderMap().keySet()) {
                lowerToOrig.put(orig.trim().toLowerCase(), orig);
            }

            // 必需列
            String[] required = {"names", "logfoldchanges", "pvals_adj", "scores"};
            for (String key : required) {
                if (!lowerToOrig.containsKey(key)) {
                    response.sendError(HttpServletResponse.SC_INTERNAL_SERVER_ERROR,
                            "CSV 表头缺少列: " + key);
                    return;
                }
            }

            String namesCol = lowerToOrig.get("names");
            String fcCol    = lowerToOrig.get("logfoldchanges");
            String pvalCol  = lowerToOrig.get("pvals_adj");
            String scoreCol = lowerToOrig.get("scores");
            String groupCol = lowerToOrig.get("group");  // 可选列

            for (CSVRecord rec : parser) {
                try {
                    String gene  = rec.get(namesCol);
                    double logfc = Double.parseDouble(rec.get(fcCol));
                    double pval  = Double.parseDouble(rec.get(pvalCol));
                    String score = rec.get(scoreCol);
                    String group = groupCol != null ? rec.get(groupCol) : "Unknown";

                    if (pval <= pvalThreshold && logfc >= fcThreshold) {
                        if (filterGroup == null || filterGroup.isEmpty() || group.equals(filterGroup)) {
                            Map<String,String> row = new HashMap<>();
                            row.put("gene",           gene);
                            row.put("logfoldchanges", String.valueOf(logfc));
                            row.put("pvals_adj",      String.valueOf(pval));
                            row.put("scores",         score);
                            row.put("group",          group); // Cell Type
                            outList.add(row);
                        }
                    }
                } catch (NumberFormatException ignore) {
                    // 忽略格式错误的行
                }
            }

        } catch (IOException e) {
            response.sendError(HttpServletResponse.SC_INTERNAL_SERVER_ERROR, "读取 CSV 出错: " + e.getMessage());
            return;
        }

        // 输出 JSON
        response.setContentType("application/json;charset=UTF-8");
        new Gson().toJson(outList, response.getWriter());
    }

    private double parseDouble(String s, double def) {
        try {
            return Double.parseDouble(s);
        } catch (Exception e) {
            return def;
        }
    }
}
