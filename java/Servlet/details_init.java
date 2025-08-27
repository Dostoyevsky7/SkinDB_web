package Servlet;

import javax.servlet.*;
import javax.servlet.http.*;
import javax.servlet.annotation.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONObject;

import Utils.CSVUtils;

@WebServlet(name = "details_init", value = "/details_init")
public class details_init extends HttpServlet {
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        // 当有HTTP GET请求到达/getIP路径时触发

        //UMAP的构建
        String py_interpreter = "D:/Python/python.exe";
        String umap_scriptPath = "C:/Users/Colorful404/PycharmProjects/Single-cell test/.venv/Scripts/Main.py";
        String cellType_scriptPath = "C:/Users/Colorful404/PycharmProjects/Single-cell test/.venv/Scripts/cell_types_get.py";


        new Thread(() -> {
            try {
                ProcessBuilder pb1 = new ProcessBuilder(py_interpreter, umap_scriptPath);
                pb1.directory(new File("C:/Users/Colorful404/PycharmProjects/Single-cell test/"));
                pb1.redirectErrorStream(true);  // 合并标准输出和错误输出
                Process process = pb1.start();

                // 异步读取输出，防止缓冲区阻塞
                BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
                String line;
                while ((line = in.readLine()) != null) {
                    System.out.println(line);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }).start();


        //cell type的读取
        //之后改成从客户端获取文件数据
        String anndata_path = "C:/Users/Colorful404/PycharmProjects/Single-cell test/.venv/Scripts/merged_anndata_clustered.h5ad";
        ProcessBuilder pb2 = new ProcessBuilder(py_interpreter, cellType_scriptPath,anndata_path);
        Process process = pb2.start();
        //获取py的输出
        BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
        StringBuilder output = new StringBuilder();
        String line;
        while ((line = reader.readLine()) != null) {
            output.append(line);
        }

        // 等待 Python 脚本执行完毕
        try {
            process.waitFor();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // 输出结果是 JSON 格式的字符串，解析它
        String resultString = output.toString();
        JSONArray resultArray = new JSONArray(resultString);
        System.out.println(resultArray);

        request.setAttribute("cell_type_array", resultArray);
        request.getRequestDispatcher("/details.jsp").forward(request,response);

        //DEG的构建
        //读取csv
    }

    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

    }
}