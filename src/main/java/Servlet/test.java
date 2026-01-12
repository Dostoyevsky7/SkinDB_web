package Servlet;

import javax.servlet.*;
import javax.servlet.http.*;
import javax.servlet.annotation.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;

import Utils.CSVUtils;
import com.opencsv.exceptions.CsvValidationException;
import org.json.JSONArray;

@WebServlet(name = "test", value = "/test")
public class test extends HttpServlet {
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        System.out.println("working");
        List<List<String>> csvData;
        try {
            csvData = CSVUtils.readAndFilterCSV("D:/project data/DEG/GSM3758115_DEGs.csv", "0");
        } catch (CsvValidationException e) {
            throw new RuntimeException(e);
        }
        System.out.println(csvData.get(0));
        System.out.println(csvData.get(1));
        System.out.println(csvData.get(2));
        System.out.println(csvData.get(3));

        String py_interpreter = "D:/Python/python.exe";
        String umap_scriptPath = "C:/Users/Colorful404/PycharmProjects/Single-cell test/.venv/Scripts/Main.py";
        String cellType_scriptPath = "C:/Users/Colorful404/PycharmProjects/Single-cell test/.venv/Scripts/cell_types_get.py";

        String path = "C:/Users/Colorful404/PycharmProjects/Single-cell test/.venv/Scripts/merged_anndata_clustered.h5ad";
        ProcessBuilder pb2 = new ProcessBuilder(py_interpreter, cellType_scriptPath, path);
        System.out.println(pb2);
        Process process = pb2.start();
        //添加对错误流的读取
        BufferedReader errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
        StringBuilder errorOutput = new StringBuilder();
        String errorLine;
        while ((errorLine = errorReader.readLine()) != null) {
            errorOutput.append(errorLine).append("\n");
        }
        System.out.println("Error Output: " + errorOutput.toString());
        //获取py的输出
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        StringBuilder output = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            output.append(line);
            System.out.println("out.append");
        }

        // 等待 Python 脚本执行完毕
        try {
            process.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // 输出结果是 JSON 格式的字符串，解析它
        System.out.println(output);
        String resultString = output.toString();
        JSONArray resultArray = new JSONArray(resultString);


    }

    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

    }
}
