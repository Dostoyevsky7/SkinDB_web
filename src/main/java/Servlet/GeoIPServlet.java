package Servlet;

import javax.servlet.*;
import javax.servlet.http.*;
import javax.servlet.annotation.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import com.maxmind.geoip2.DatabaseReader;
import com.maxmind.geoip2.model.CityResponse;

import java.net.InetAddress;
import java.time.Instant;

@WebServlet("/track")
public class GeoIPServlet extends HttpServlet {
    private String logPath;
    private DatabaseReader geoDBReader;

    @Override
    public void init() throws ServletException {
        try {
            // 从WEB-INF目录加载数据库
            String dbPath = getServletContext().getRealPath("/WEB-INF/GeoLite2-City.mmdb");
            geoDBReader = new DatabaseReader.Builder(new File(dbPath)).build();
            String configuredPath = getServletContext().getInitParameter("geoLogPath");
            if (configuredPath != null && !configuredPath.trim().isEmpty()) {
                logPath = configuredPath;
            } else {
                String realPath = getServletContext().getRealPath("/WEB-INF/access_geo.log");
                if (realPath != null && !realPath.trim().isEmpty()) {
                    logPath = realPath;
                } else {
                    logPath = System.getProperty("java.io.tmpdir")
                            + File.separator + "access_geo.log";
                }
            }
        } catch (IOException e) {
            throw new ServletException("Failed to load GeoIP database", e);
        }
    }

    protected void doGet(HttpServletRequest request, HttpServletResponse response) {
        // 获取真实IP
        String ip = getClientIP(request);

        // 获取地理位置
        String location = "Unknown";
        try {
            CityResponse city = geoDBReader.city(InetAddress.getByName(ip));
            location = String.format("%s, %s",
                    city.getCity().getName(),
                    city.getCountry().getName());
        } catch (Exception e) {
            // 记录查询失败
        }

        // 记录日志
        String userAgent = request.getHeader("User-Agent");
        if (userAgent == null) {
            userAgent = "-";
        }
        String logEntry = String.format("%s\t[%s]\t%s\t%s%n",
                Instant.now().toString(),
                ip,
                location,
                userAgent);

        writeLog(logEntry);
    }

    private synchronized void writeLog(String content) {
        try (FileWriter fw = new FileWriter(logPath, true)) {
            fw.write(content);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private String getClientIP(HttpServletRequest request) {
        String ip = request.getHeader("X-Forwarded-For");
        if (ip == null || ip.isEmpty() || "unknown".equalsIgnoreCase(ip)) {
            ip = request.getRemoteAddr();
        }
        // 处理IPv6映射地址
        return ip.replaceAll("^::ffff:(\\d+\\.\\d+\\.\\d+\\.\\d+)$", "$1");
    }
}
