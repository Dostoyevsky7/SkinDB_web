package Utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DataFilter {
    /**
     * 对 CSV 数据进行过滤
     * @param csvData 解析后的 CSV 数据，每一行数据为 Map，键为列名
     * @param group 用户选择的组号（m），字符串形式，可转为整数
     * @param logfcStr 滑条对应的 logFC 阈值，字符串形式（例如 "0.5"）
     * @param pAdjStr 滑条对应的 adjusted p-value 阈值，字符串形式（例如 "0.05"）
     * @param gene 用户选择的基因名，如选择“All”则不过滤基因名
     * @return 返回过滤后的数据列表，每个 Map 只包含基因名称、adjusted p-value 和 log fold change
     */
    public static List<Map<String, Object>> filter(List<Map<String, String>> csvData, String group, String logfcStr, String pAdjStr, String gene) {
        int m;
        double logfcThreshold;
        double pAdjThreshold;

        try {
            m = Integer.parseInt(group);
            logfcThreshold = Double.parseDouble(logfcStr);
            pAdjThreshold = Double.parseDouble(pAdjStr);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("无法解析 group、logfc 或 p_adj 参数，请检查输入格式", e);
        }

        // 构建需要读取的列名
        String nameKey = m + "_names";
        String pAdjKey = m + "_pvals_adj";
        String logfcKey = m + "_logfoldchanges";

        List<Map<String, Object>> filteredList = new ArrayList<>();

        for (Map<String, String> row : csvData) {
            // 获取当前行中指定组的数据
            String geneName = row.get(nameKey);
            String pAdjValueStr = row.get(pAdjKey);
            String logfcValueStr = row.get(logfcKey);

            // 数据可能缺失或格式不对，跳过这行
            if (geneName == null || pAdjValueStr == null || logfcValueStr == null) {
                continue;
            }

            // 如果用户指定了基因名且不为“All”，则进行基因名匹配（忽略大小写）
            if (!"All".equalsIgnoreCase(gene) && !geneName.equalsIgnoreCase(gene)) {
                continue;
            }

            double pAdjValue;
            double logfcValue;
            try {
                pAdjValue = Double.parseDouble(pAdjValueStr);
                logfcValue = Double.parseDouble(logfcValueStr);
            } catch (NumberFormatException e) {
                // 格式错误时跳过当前行
                continue;
            }

            // 筛选条件：
            // 1. p_vals_adj 小于用户输入的 pAdjThreshold
            // 2. logfoldchanges 大于用户输入的 logfcThreshold
            if (pAdjValue >= pAdjThreshold && logfcValue <= logfcThreshold) {
                continue;
            }


            // 如果满足条件，则构造过滤后的记录
            Map<String, Object> filteredRow = new HashMap<>();
            filteredRow.put("geneName", geneName);
            filteredRow.put("pAdj", pAdjValue);
            filteredRow.put("logfc", logfcValue);

            filteredList.add(filteredRow);
        }

        return filteredList;
    }
}
