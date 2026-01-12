package Servlet;

import javax.servlet.*;
import javax.servlet.http.*;
import javax.servlet.annotation.*;
import java.io.*;
import java.net.URLEncoder; // Import for URL encoding
import java.util.UUID;
import com.google.gson.Gson; // You will need the Gson library for JSON handling

@WebServlet("/integrate")
public class IntegrateServlet extends HttpServlet {
    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
        // 1. Get the list of selected SAIDs from the request
        String[] saids = request.getParameterValues("saids[]");

        if (saids == null || saids.length < 2) {
            response.setContentType("application/json");
            response.setCharacterEncoding("UTF-8");
            HttpJsonResponse jsonResponse = new HttpJsonResponse();
            jsonResponse.error = "Please select at least two datasets to integrate.";
            response.getWriter().write(new Gson().toJson(jsonResponse));
            return;
        }

        // Construct the URL for the Dash application
        // Assuming your Dash app runs on the same host but potentially a different port,
        // or a completely different domain if deployed separately.
        // For simplicity, let's assume it's accessible via /dash on the same server.
        // You might need to adjust the port (8050) if your Dash app runs on a different one.
        StringBuilder dashAppUrl = new StringBuilder("http://localhost:8050/dash?"); // Or your server's actual IP/domain
        dashAppUrl.append("saids=");
        for (int i = 0; i < saids.length; i++) {
            dashAppUrl.append(URLEncoder.encode(saids[i], "UTF-8"));
            if (i < saids.length - 1) {
                dashAppUrl.append(",");
            }
        }

        // 5. Send the response back to the frontend
        response.setContentType("application/json");
        response.setCharacterEncoding("UTF-8");

        HttpJsonResponse jsonResponse = new HttpJsonResponse();
        // Instead of returning an image URL, return the URL to the Dash app
        jsonResponse.redirectUrl = dashAppUrl.toString(); // Use a new field for redirection
        jsonResponse.message = "Successfully prepared Dash app URL.";

        response.getWriter().write(new Gson().toJson(jsonResponse));
    }

    // Helper class for JSON response
    private static class HttpJsonResponse {
        String redirectUrl; // New field for the Dash app URL
        String error;
        String message; // Optional: for success messages
    }
}