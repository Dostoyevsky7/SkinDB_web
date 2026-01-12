package Utils;

import com.opencsv.CSVReader;
import com.opencsv.exceptions.CsvValidationException;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;

import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.lang.reflect.Array;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class CSVUtils {
    public static List<List<String>> readAndFilterCSV(String filePath, String cellType) throws IOException, CsvValidationException {
        // 存储最终筛选后的数据
        List<List<String>> result = new ArrayList<>();

        // 使用CSVReader读取文件
        CSVReader reader = new CSVReader(new FileReader(filePath));

        // 读取表头
        String[] st_header = reader.readNext();
        List<String> header = Arrays.asList(st_header);
        // 找到目标表头
        int list_position = 0;
        for (int i = 0; i < header.size(); i++){
            if (header.get(i).contains(cellType)){
                list_position = i;
                break;
            }
        }

        // 读取数据并按条件筛选
        String[] line = new String[0];

        // 遍历表头后的每一行
        while ((line = reader.readNext()) != null){
            // 提取从第k项到第k+3项的数据（确保不超出数组长度）
            if (line.length >= list_position + 4){
                List<String> extractedData = new ArrayList<>();
                extractedData.add(cellType);
                for (int i = list_position; i < list_position + 4;i++){
                    if (i != list_position + 1){
                        extractedData.add(line[i]);
                    }
                }
                result.add(extractedData);
            }
        }

        // 关闭reader
        reader.close();

        return result;
    }
}
