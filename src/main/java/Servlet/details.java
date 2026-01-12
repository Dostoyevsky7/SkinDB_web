package Servlet;

import javax.servlet.*;
import javax.servlet.http.*;
import javax.servlet.annotation.*;
import java.io.IOException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import javax.servlet.ServletException;
import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

@WebServlet(name = "details", value = "/details")
public class details extends HttpServlet {
    @Override
    protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

    }

    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        String cell_type = request.getParameter("browse");
        String py_interpreter = "D:/Python/python.exe";
        String scriptPath = "C:/Users/Colorful404/PycharmProjects/Single-cell test/.venv/Scripts/Main.py";

        ProcessBuilder pb = new ProcessBuilder(py_interpreter, scriptPath, cell_type, "gene_id/genename");
        pb.start();

    }
}
