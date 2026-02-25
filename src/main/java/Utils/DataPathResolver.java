package Utils;

import javax.servlet.ServletContext;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

public final class DataPathResolver {
    private static final String CONTEXT_DATA_ROOT = "skindbDataRoot";
    private static final String PROPERTY_DATA_ROOT = "skindb.data.root";
    private static final String ENV_DATA_ROOT = "SKINDB_DATA_ROOT";

    private static final String CONTEXT_PYTHON = "skindbPython";
    private static final String PROPERTY_PYTHON = "skindb.python";
    private static final String ENV_PYTHON = "SKINDB_PYTHON";

    private static final String[] DEFAULT_DATA_ROOTS = {"/opt/SkinDB", "/root/SkinDB"};
    private static final String[] DEFAULT_PYTHON_COMMANDS = {
            "/opt/miniconda3/envs/scrna/bin/python",
            "/root/miniconda3/envs/scrna/bin/python",
            "python3"
    };

    private DataPathResolver() {}

    public static String resolveDataRoot(ServletContext context) {
        List<String> candidates = getDataRootCandidates(context);
        for (String candidate : candidates) {
            File dir = new File(candidate);
            if (dir.isDirectory() && dir.canRead()) {
                return dir.getAbsolutePath();
            }
        }
        return candidates.get(0);
    }

    public static List<String> getCandidateFilePaths(ServletContext context, String relativePath) {
        List<String> result = new ArrayList<String>();
        for (String root : getDataRootCandidates(context)) {
            result.add(new File(root, normalizeRelativePath(relativePath)).getAbsolutePath());
        }
        return result;
    }

    public static File resolveReadableFile(ServletContext context, String relativePath) {
        String normalizedRelativePath = normalizeRelativePath(relativePath);
        List<String> dataRoots = getDataRootCandidates(context);
        for (String root : dataRoots) {
            File file = new File(root, normalizedRelativePath);
            if (file.exists() && file.canRead()) {
                return file;
            }
        }
        return new File(dataRoots.get(0), normalizedRelativePath);
    }

    public static String resolvePythonCommand(ServletContext context) {
        LinkedHashSet<String> candidates = new LinkedHashSet<String>();
        addIfPresent(candidates, context == null ? null : context.getInitParameter(CONTEXT_PYTHON));
        addIfPresent(candidates, System.getProperty(PROPERTY_PYTHON));
        addIfPresent(candidates, System.getenv(ENV_PYTHON));
        for (String fallback : DEFAULT_PYTHON_COMMANDS) {
            addIfPresent(candidates, fallback);
        }

        for (String command : candidates) {
            if ("python3".equals(command)) {
                continue;
            }
            File executable = new File(command);
            if (executable.exists() && executable.canExecute()) {
                return executable.getAbsolutePath();
            }
        }
        return "python3";
    }

    private static List<String> getDataRootCandidates(ServletContext context) {
        Set<String> candidates = new LinkedHashSet<String>();
        addIfPresent(candidates, context == null ? null : context.getInitParameter(CONTEXT_DATA_ROOT));
        addIfPresent(candidates, System.getProperty(PROPERTY_DATA_ROOT));
        addIfPresent(candidates, System.getenv(ENV_DATA_ROOT));
        for (String fallback : DEFAULT_DATA_ROOTS) {
            addIfPresent(candidates, fallback);
        }
        return new ArrayList<String>(candidates);
    }

    private static String normalizeRelativePath(String relativePath) {
        if (relativePath == null) {
            return "";
        }
        return relativePath.startsWith("/") ? relativePath.substring(1) : relativePath;
    }

    private static void addIfPresent(Set<String> values, String value) {
        if (value == null) {
            return;
        }
        String trimmed = value.trim();
        if (!trimmed.isEmpty()) {
            values.add(trimmed);
        }
    }
}
