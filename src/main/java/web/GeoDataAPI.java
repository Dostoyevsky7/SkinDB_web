package web;

import com.fasterxml.jackson.databind.ObjectMapper;

import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Deque;
import java.util.stream.Stream;

@WebServlet("/api/geo")
public class GeoDataAPI extends HttpServlet {
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException {
        // 读取日志文件最后100条记录
        List<Map<String, String>> data = new ArrayList<>();
        Path logPath = resolveLogPath(request);
        if (logPath != null && Files.exists(logPath)) {
            Deque<String> tail = new ArrayDeque<>(100);
            try (Stream<String> lines = Files.lines(logPath)) {
                lines.forEach(line -> {
                    if (tail.size() == 100) {
                        tail.removeFirst();
                    }
                    tail.addLast(line);
                });
            }

            for (String line : tail) {
                Map<String, String> entry = parseLine(line);
                if (entry != null) {
                    data.add(entry);
                }
            }
        }

        response.setContentType("application/json");
        new ObjectMapper().writeValue(response.getWriter(), data);
    }

    private Path resolveLogPath(HttpServletRequest request) {
        String configured = request.getServletContext().getInitParameter("geoLogPath");
        if (configured != null && !configured.trim().isEmpty()) {
            return Paths.get(configured);
        }
        String realPath = request.getServletContext().getRealPath("/WEB-INF/access_geo.log");
        if (realPath != null && !realPath.trim().isEmpty()) {
            return Paths.get(realPath);
        }
        String tmp = System.getProperty("java.io.tmpdir");
        if (tmp != null && !tmp.trim().isEmpty()) {
            return Paths.get(tmp, "access_geo.log");
        }
        return null;
    }

    private Map<String, String> parseLine(String line) {
        if (line == null || line.trim().isEmpty()) {
            return null;
        }
        String ip = null;
        String location = "Unknown";

        if (line.contains("\t")) {
            String[] parts = line.split("\t", 4);
            if (parts.length >= 3) {
                ip = parts[1].replaceAll("[\\[\\]]", "");
                location = parts[2];
            }
        } else {
            int l = line.indexOf('[');
            int r = line.indexOf(']', l + 1);
            if (l >= 0 && r > l) {
                ip = line.substring(l + 1, r);
                String remainder = line.substring(r + 1).trim();
                if (!remainder.isEmpty()) {
                    location = remainder;
                }
            }
        }

        if (ip == null || ip.isEmpty()) {
            return null;
        }
        Map<String, String> entry = new HashMap<>();
        entry.put("ip", ip);
        entry.put("location", location);
        return entry;
    }
}
