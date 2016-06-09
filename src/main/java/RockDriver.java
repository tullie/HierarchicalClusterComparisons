import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import iweb2.ch4.clustering.rock.ROCKAlgorithm;
import iweb2.ch4.model.DataPoint;
import iweb2.ch4.clustering.hierarchical.Dendrogram;

class RockDriver {
  static final double THRESHOLD_VALUE = 0.2;

  public static void main(String args[]) throws IOException {
    if (args.length != 2) {
      System.out.println("Dataset filename and cluster count required");
      return;
    }

    System.out.println("Running ROCK driver...");

    List<double[]> dataset = Utils.readDataset(args[0]);
    DataPoint[] points = new DataPoint[dataset.size()];
    for (int i = 0; i < dataset.size(); ++i) {
      points[i] = new DataPoint(Integer.toString(i), dataset.get(i));
    }

    ROCKAlgorithm rockClusterer =
        new ROCKAlgorithm(points, Integer.parseInt(args[1]), THRESHOLD_VALUE);
    Dendrogram dendrogram = rockClusterer.cluster();
    dendrogram.printAll();
  }
}
