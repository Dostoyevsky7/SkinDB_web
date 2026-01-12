package Servlet;

import Utils.CSVUtils;
import com.google.gson.Gson;
import com.opencsv.exceptions.CsvValidationException;

import javax.servlet.*;
import javax.servlet.http.*;
import javax.servlet.annotation.*;
import java.io.IOException;
import java.util.List;

@WebServlet(name = "deg_cell_type_change", value = "/deg_cell_type_change")
public class deg_cell_type_change extends HttpServlet {
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        //获得cellType
        String cellType = request.getParameter("cellType");
        List<List<String>> result;
        System.out.println(cellType);
        /**
         * @param filePath 之后改成从服务器获取文件名
         */
        try {
            result = CSVUtils.readAndFilterCSV("D:/project data/DEG/GSM3758115_DEGs.csv", cellType);
        } catch (CsvValidationException e) {
            System.out.println("can not open the csv file");
            throw new RuntimeException(e);
        }
        Gson gson = new Gson();
        String json = gson.toJson(result);
        response.setContentType("application/json;charset=UTF-8");
        response.getWriter().write(json);
        response.getWriter().flush();
    }

    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

    }
}
