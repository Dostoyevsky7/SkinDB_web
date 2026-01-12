package Servlet;

import Entity.Dataset;
import Utils.DatasetLoader;
import javax.servlet.http.HttpServlet;
import javax.servlet.annotation.WebServlet;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;


import java.io.IOException;
import java.util.List;

@WebServlet("/SearchServlet")
public class SearchServlet extends HttpServlet {
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        try {
            String excelPath = getServletContext().getRealPath("/WEB-INF/IntegrateTable.xlsx");
            List<Dataset> dataList = DatasetLoader.load(excelPath);
            request.setAttribute("datasetList", dataList);
            request.getRequestDispatcher("search.jsp").forward(request, response);
        } catch (Exception e) {
            throw new ServletException("加载数据出错", e);
        }
    }
}
