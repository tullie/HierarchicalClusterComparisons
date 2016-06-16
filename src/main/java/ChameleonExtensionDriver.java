
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import chameleonextension.Chameleon;
import chameleonextension.Node;

class ChameleonExtensionDriver {
  public static void main(String args[]) throws IOException {
    if (args.length != 2) {
      System.out.println("Dataset filename and cluster count required");
      return;
    }

    System.out.println("Running Chameleon extension driver...");

    List<double[]> dataset = Utils.readDataset(args[0]);
    Chameleon chameleonExt = new Chameleon();
    List<List<Node>> clusters =
        chameleonExt.cluster(dataset, Integer.parseInt(args[1]));

    PrintWriter writer =
        new PrintWriter(args[0] + ".chameleon_ext_result", "UTF-8");
    int clusterNumber = 0;
    for (List<Node> cluster : clusters) {
      for (Node node : cluster) {
        String entryLine = "";
        for (double dataPoint : node.values) {
          entryLine += dataPoint + " ";
        }
        writer.println(entryLine + clusterNumber);
      }
      clusterNumber++;
    }
    writer.close();
  }
}
