import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import chameleon.Chameleon;

class ChameleonDriver {
  public static void main(String args[]) throws IOException {
    if (args.length != 2) {
      System.out.println("Dataset filename and partition count required");
      return;
    }

    System.out.println("Running Chameleon driver...");

    List<double[]> dataset = Utils.readDataset(args[0]);
    Chameleon.cluster(dataset, Integer.parseInt(args[1]));
  }
}
