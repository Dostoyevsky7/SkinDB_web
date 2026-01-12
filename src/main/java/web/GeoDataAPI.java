package web;

import Servlet.GeoIPServlet;
import com.fasterxml.jackson.databind.ObjectMapper;

import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import static Servlet.GeoIPServlet.LOG_PATH;

@WebServlet("/api/geo")
public class GeoDataAPI extends HttpServlet {
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws IOException {
        // 读取日志文件最后100条记录
        List<Map<String, String>> data = new ArrayList<>();
        try (Stream<String> lines = Files.lines(Paths.get(LOG_PATH))) {
            lines.skip(Math.max(0, Files.size(Paths.get(LOG_PATH)) - 100))
                    .forEach(line -> {
                        // 解析日志格式：日期 [IP] 位置 UserAgent
                        String[] parts = line.split(" ", 4);
                        Map<String, String> entry = new HashMap<>();
                        entry.put("ip", parts[1].replaceAll("[\\[\\]]", ""));
                        entry.put("location", parts[2]);
                        data.add(entry);
                    });
        }

        response.setContentType("application/json");
        new ObjectMapper().writeValue(response.getWriter(), data);
    }
}
