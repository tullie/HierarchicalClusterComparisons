import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;

public class Utils {
  static List<double[]> readDataset(String datasetFile) throws IOException {
    List<double[]> dataset = new ArrayList<>();
    BufferedReader fileIn = new BufferedReader(new FileReader(datasetFile));
    String line = null;
    while ((line = fileIn.readLine()) != null) {
      String[] lineSplit = line.replaceFirst("^\\s+", "").split("\\s+");
      double[] values = new double[lineSplit.length];
      for (int i = 0; i < values.length; ++i) {
        values[i] = Double.parseDouble(lineSplit[i]);
      }
      dataset.add(values);
    }
    fileIn.close();
    return dataset;
  }
}
