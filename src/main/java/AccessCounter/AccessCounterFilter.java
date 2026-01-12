package AccessCounter;

import javax.servlet.*;
import javax.servlet.annotation.WebFilter;
import java.io.*;
import java.time.LocalDate;
import java.util.Properties;
import java.util.concurrent.atomic.AtomicInteger;

import javax.servlet.http.Cookie;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

@WebFilter("/*") // 拦截所有请求
public class AccessCounterFilter implements Filter {
    private String counterFilePath;
    private FilterConfig filterConfig;

    @Override
    public void init(FilterConfig filterConfig) throws ServletException {
        this.filterConfig = filterConfig;
        ServletContext context = filterConfig.getServletContext();

        // 从 web.xml 配置中获取保存文件的路径（可以在 web.xml 中配置 <init-param>），否则使用默认路径
        counterFilePath = filterConfig.getInitParameter("counterFile");
        if (counterFilePath == null || counterFilePath.trim().isEmpty()) {
            // 注意：这里的默认路径需要根据你的服务器环境设置，确保有写权限
            counterFilePath = "D:/counterFilter/counter.properties";
        }
        context.setAttribute("counterFile", counterFilePath);

        // 尝试从文件中读取已保存的计数数据
        int total = 0;
        int daily = 0;
        String savedDate = null;
        Properties props = new Properties();
        File file = new File(counterFilePath);
        if (file.exists()) {
            try (FileInputStream fis = new FileInputStream(file)) {
                props.load(fis);
                total = Integer.parseInt(props.getProperty("totalCount", "0"));
                daily = Integer.parseInt(props.getProperty("dailyCount", "0"));
                savedDate = props.getProperty("currentDate");
            } catch (IOException | NumberFormatException e) {
                e.printStackTrace();
            }
        }
        LocalDate today = LocalDate.now();
        // 如果保存的日期与今天不一致，则重置每日计数
        if (savedDate == null || !savedDate.equals(today.toString())) {
            daily = 0;
        }
        context.setAttribute("totalCount", new AtomicInteger(total));
        context.setAttribute("dailyCount", new AtomicInteger(daily));
        context.setAttribute("currentDate", today);
        System.out.println("init: totalCount=" + total + ", dailyCount=" + daily);
    }

    @Override
    public void doFilter(ServletRequest request, ServletResponse response, FilterChain chain)
            throws IOException, ServletException {
        HttpServletRequest req = (HttpServletRequest) request;
        HttpServletResponse res = (HttpServletResponse) response;
        ServletContext context = filterConfig.getServletContext();
        AtomicInteger totalCount = (AtomicInteger) context.getAttribute("totalCount");
        AtomicInteger dailyCount = (AtomicInteger) context.getAttribute("dailyCount");
        LocalDate currentDate = (LocalDate) context.getAttribute("currentDate");

        LocalDate today = LocalDate.now();

        String uri = req.getRequestURI();
        if (uri.endsWith(".css") || uri.endsWith(".js") || uri.endsWith(".png") || uri.endsWith(".jpg") || uri.endsWith(".gif") || uri.contains("favicon")) {
            chain.doFilter(request, response);
            return;
        }
        // 如果日期变了，则重置dailyCount并更新当前日期
        if (!today.equals(currentDate)) {
            dailyCount.set(0);
            context.setAttribute("currentDate", today);
        }
        Cookie[] cookies = req.getCookies();
        Cookie countCookie = null;
        if (cookies != null) {
            for (Cookie c : cookies) {
                if ("count_cookie".equals(c.getName())) {
                    countCookie = c;
                    break;
                }
            }
        }

        // 如果不存在则认为是新访问，进行计数，并创建Cookie
        if (countCookie == null) {
            totalCount.incrementAndGet();
            dailyCount.incrementAndGet();
            System.out.println("doFilter: totalCount=" + totalCount.get() + ", dailyCount=" + dailyCount.get());

            // 创建 Cookie（例如设置值为 "1" 表示已计数过）
            countCookie = new Cookie("count_cookie", "1");
            countCookie.setPath(req.getContextPath());
            // 不设置 maxAge 则为 session Cookie（浏览器关闭后失效）
            res.addCookie(countCookie);
        }
        // 将请求继续传递到下一个过滤器或目标资源
        chain.doFilter(request, response);
    }

    @Override
    public void destroy() {
        // 在 Filter 销毁时将当前计数数据写入文件保存
        ServletContext context = filterConfig.getServletContext();
        AtomicInteger totalCount = (AtomicInteger) context.getAttribute("totalCount");
        AtomicInteger dailyCount = (AtomicInteger) context.getAttribute("dailyCount");
        LocalDate currentDate = (LocalDate) context.getAttribute("currentDate");
        String counterFile = (String) context.getAttribute("counterFile");

        if (totalCount.get() != 0 && dailyCount.get() != 0){
            Properties props = new Properties();
            props.setProperty("totalCount", String.valueOf(totalCount.get()));
            props.setProperty("dailyCount", String.valueOf(dailyCount.get()));
            props.setProperty("currentDate", currentDate.toString());
            props.setProperty("destory?", "yes");
            System.out.println("destroy: totalCount=" + totalCount.get() + ", dailyCount=" + dailyCount.get());

            try (FileOutputStream fos = new FileOutputStream(counterFile)) {
                props.store(fos, "Access Counter Properties");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


    }
}