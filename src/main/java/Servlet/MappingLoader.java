package Servlet;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

import java.io.InputStreamReader;
import java.util.Map;

public class MappingLoader {

    private static Map<String, Map<String, String>> mapping;

    static {
        try {
            InputStreamReader reader = new InputStreamReader(
                    MappingLoader.class.getClassLoader().getResourceAsStream("mapping.json"),
                    "UTF-8"
            );
            mapping = new Gson().fromJson(reader,
                    new TypeToken<Map<String, Map<String, String>>>() {}.getType());
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("Failed to load mapping.json", e);
        }
    }

    public static Map<String, Map<String, String>> load() {
        return mapping;
    }
}
