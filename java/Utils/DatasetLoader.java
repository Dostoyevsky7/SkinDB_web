package Utils;

import Entity.Dataset;
import org.apache.poi.ss.usermodel.*;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.List;

public class DatasetLoader {
    public static List<Dataset> load(String excelPath) throws Exception {
        List<Dataset> list = new ArrayList<>();
        FileInputStream fis = new FileInputStream(excelPath);
        Workbook wb = WorkbookFactory.create(fis);
        Sheet sheet = wb.getSheetAt(0);

        String currentGSE = "", currentTitle = "", currentSummary = "";

        for (int i = 1; i <= sheet.getLastRowNum(); i++) {
            Row row = sheet.getRow(i);
            if (row == null) continue;

            Dataset d = new Dataset();
            d.setSaid(getCell(row, 0));

            String gse = getCell(row, 1);
            if (!gse.isEmpty()) currentGSE = gse;
            d.setGse(currentGSE);

            d.setGsm(getCell(row, 9));

            String title = getCell(row, 4);
            if (!title.isEmpty()) currentTitle = title;
            d.setTitle(currentTitle);

            String summary = getCell(row, 5);
            if (!summary.isEmpty()) currentSummary = summary;
            d.setSummary(currentSummary);

            d.setFile(d.getSaid() + ".h5ad");
            list.add(d);
        }
        wb.close();
        return list;
    }

    private static String getCell(Row row, int index) {
        try {
            return row.getCell(index).toString().trim();
        } catch (Exception e) {
            return "";
        }
    }
}
