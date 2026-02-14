package Servlet;

import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.*;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;
import org.apache.commons.csv.*;

@WebServlet("/deg")
public class DEGServlet extends HttpServlet {

    /**
     * ✅ 新：Linux 服务器上 mouse DEG 根目录（按你给的真实结构）
     * /root/SkinDB/mouse/mouseGSE/GSE.../GSM.../<GSM>_annotated_all_DEGs.csv
     *
     * 如果 Tomcat 无法读 /root，请改成你放数据的可读目录，例如：
     * private static final String MOUSE_DEG_BASE = "/data/SkinDB/mouse/mouseGSE";
     */
    private static final String MOUSE_DEG_BASE = "/root/SkinDB/mouse/mouseGSE";

    /**
     * 兼容旧版本：SAID -> csv_path 的映射（mapping.csv.json）
     * 注意：本版本不会因为 mapping 文件不存在而启动失败；仅在缺参时回退使用。
     */
    private static final Map<String, Map<String, String>> mapping = loadOptionalMapping();

    private static Map<String, Map<String, String>> loadOptionalMapping() {
        try (InputStream is = DEGServlet.class.getClassLoader().getResourceAsStream("mapping.csv.json")) {
            if (is == null) {
                // mapping 不存在就返回空 Map；不会让 Tomcat 启动失败
                return new HashMap<>();
            }
            try (InputStreamReader reader = new InputStreamReader(is, StandardCharsets.UTF_8)) {
                Map<String, Map<String, String>> m = new Gson().fromJson(
                        reader,
                        new TypeToken<Map<String, Map<String, String>>>() {}.getType()
                );
                return (m != null) ? m : new HashMap<>();
            }
        } catch (Exception e) {
            // 解析失败也不让启动失败
            return new HashMap<>();
        }
    }

    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response)
            throws ServletException, IOException {

        // -------------------------
        // 0) 读取参数
        // -------------------------
        String said = safe(request.getParameter("said"));
        String gse  = safe(request.getParameter("gse"));
        String gsm  = safe(request.getParameter("gsm"));

        double pvalThreshold = parseDouble(request.getParameter("pval"), 0.05);
        double fcThreshold   = parseDouble(request.getParameter("fc"), 1.0);
        String filterGroup   = safe(request.getParameter("group")); // 可为空

        // -------------------------
        // 1) 定位 CSV 文件路径（优先使用 gse+gsm）
        // -------------------------
        File csvFile = resolveCsvFile(gse, gsm, said, request);

        if (csvFile == null) {
            // 参数不够且 mapping 也没有
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            response.setContentType("application/json;charset=UTF-8");
            response.getWriter().write("{\"error\":\"Missing parameters. Provide gse & gsm (recommended), or said with valid mapping.\"}");
            return;
        }

        if (!csvFile.exists()) {
            response.setStatus(HttpServletResponse.SC_NOT_FOUND);
            response.setContentType("application/json;charset=UTF-8");
            response.getWriter().write("{\"error\":\"CSV file not found\",\"path\":\"" + escJson(csvFile.getAbsolutePath()) + "\"}");
            return;
        }

        if (!csvFile.canRead()) {
            response.setStatus(HttpServletResponse.SC_FORBIDDEN);
            response.setContentType("application/json;charset=UTF-8");
            response.getWriter().write("{\"error\":\"CSV file not readable (permission issue)\",\"path\":\"" + escJson(csvFile.getAbsolutePath()) + "\"}");
            return;
        }

        // -------------------------
        // 2) 读取 CSV + 过滤
        // -------------------------
        List<Map<String,String>> outList = new ArrayList<>();

        // 用 commons-csv，稳
        try (
                Reader in = Files.newBufferedReader(csvFile.toPath(), StandardCharsets.UTF_8);
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

            // 兼容列名：names/gene；scores/score；group/Group 等
            String namesCol = pickCol(lowerToOrig, "names", "gene", "Gene", "symbol", "Symbol");
            String fcCol    = pickCol(lowerToOrig, "logfoldchanges", "logfc", "log2fc", "avg_log2fc");
            String pvalCol  = pickCol(lowerToOrig, "pvals_adj", "p_val_adj", "padj", "adj_pval", "pvalue", "p_val");
            String scoreCol = pickCol(lowerToOrig, "scores", "score", "Score");
            String groupCol = pickCol(lowerToOrig, "group", "Group", "celltype", "cell_type");

            // 必需列检查
            if (namesCol == null || fcCol == null || pvalCol == null || scoreCol == null) {
                response.setStatus(HttpServletResponse.SC_INTERNAL_SERVER_ERROR);
                response.setContentType("application/json;charset=UTF-8");
                response.getWriter().write("{\"error\":\"CSV header missing required columns\","
                        + "\"need\":\"(names/gene), logfoldchanges, pvals_adj, scores\","
                        + "\"headers\":\"" + escJson(String.join("|", parser.getHeaderMap().keySet())) + "\","
                        + "\"path\":\"" + escJson(csvFile.getAbsolutePath()) + "\"}");
                return;
            }

            for (CSVRecord rec : parser) {
                try {
                    String gene  = rec.get(namesCol);
                    double logfc = Double.parseDouble(rec.get(fcCol));
                    double pval  = Double.parseDouble(rec.get(pvalCol));
                    String score = rec.get(scoreCol);
                    String group = (groupCol != null) ? rec.get(groupCol) : "Unknown";

                    // 过滤逻辑：pval <= threshold 且 logfc >= threshold（保持你原来逻辑）
                    if (pval <= pvalThreshold && logfc >= fcThreshold) {
                        if (filterGroup == null || filterGroup.isEmpty() || group.equals(filterGroup)) {
                            Map<String,String> row = new HashMap<>();
                            row.put("gene",           gene);
                            row.put("logfoldchanges", String.valueOf(logfc));
                            row.put("pvals_adj",      String.valueOf(pval));
                            row.put("scores",         score);
                            row.put("group",          group);
                            outList.add(row);
                        }
                    }
                } catch (NumberFormatException ignore) {
                    // 忽略格式错误行
                } catch (IllegalArgumentException ignore) {
                    // 缺列等
                }
            }

        } catch (Exception e) {
            response.setStatus(HttpServletResponse.SC_INTERNAL_SERVER_ERROR);
            response.setContentType("application/json;charset=UTF-8");
            response.getWriter().write("{\"error\":\"Read CSV failed\",\"message\":\"" + escJson(e.getMessage()) + "\",\"path\":\"" + escJson(csvFile.getAbsolutePath()) + "\"}");
            return;
        }

        // -------------------------
        // 3) 输出 JSON
        // -------------------------
        response.setContentType("application/json;charset=UTF-8");
        new Gson().toJson(outList, response.getWriter());
    }

    /**
     * 优先：gse + gsm -> /root/SkinDB/mouse/mouseGSE/GSE.../GSM.../<GSM>_annotated_all_DEGs.csv
     * 回退：said -> mapping.csv.json -> getRealPath("/" + csv_path)
     */
    private File resolveCsvFile(String gse, String gsm, String said, HttpServletRequest request) {

        // 1) 新路径优先（推荐）
        if (gse != null && !gse.isEmpty() && gsm != null && !gsm.isEmpty()) {
            String fileName = gsm + "_annotated_all_DEGs.csv";
            Path p = Paths.get(MOUSE_DEG_BASE, gse, gsm, fileName);
            return p.toFile();
        }

        // 2) 回退到旧 mapping（可选）
        if (said != null && !said.isEmpty() && mapping != null && !mapping.isEmpty()) {
            Map<String, String> meta = mapping.get(said);
            if (meta != null) {
                String csvRelPath = meta.get("csv_path");
                if (csvRelPath != null && !csvRelPath.trim().isEmpty()) {
                    String realPath = request.getServletContext().getRealPath("/" + csvRelPath);
                    if (realPath != null) return new File(realPath);
                }
            }
        }

        return null;
    }

    private static String pickCol(Map<String,String> lowerToOrig, String... keys) {
        // lowerToOrig key 是小写
        for (String k : keys) {
            if (k == null) continue;
            String kk = k.trim().toLowerCase();
            if (lowerToOrig.containsKey(kk)) return lowerToOrig.get(kk);
        }
        // 再做一次忽略大小写遍历兜底
        for (String k : keys) {
            if (k == null) continue;
            for (Map.Entry<String,String> e : lowerToOrig.entrySet()) {
                if (e.getKey() != null && e.getKey().equalsIgnoreCase(k.trim())) {
                    return e.getValue();
                }
            }
        }
        return null;
    }

    private static String safe(String s) {
        if (s == null) return null;
        s = s.trim();
        return s.isEmpty() ? null : s;
    }

    private double parseDouble(String s, double def) {
        try {
            if (s == null) return def;
            return Double.parseDouble(s);
        } catch (Exception e) {
            return def;
        }
    }

    private static String escJson(String s) {
        if (s == null) return "";
        return s.replace("\\", "\\\\").replace("\"", "\\\"");
    }
}
