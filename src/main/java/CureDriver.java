import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import cure.Cure;

class CureDriver {
  static final int MIN_REPRESENTIVE_COUNT = 4;
  static final double SHRINK_FACTOR = 0.5;
  static final double REQUIRED_REPRESENTATION_PROBABILITY = 0.1;
  static final int REDUCING_FACTOR_FOR_EACH_PARTITION = 2;
  static final int NUMBER_OF_PARTITIONS = 2;

  public static void main(String args[]) throws IOException {
    if (args.length != 2) {
      System.out.println("Dataset filename and partition count required");
      return;
    }

    System.out.println("Running CURE driver...");

    List<double[]> dataset = Utils.readDataset(args[0]);
    Cure cureClusterer = new Cure(dataset, Integer.parseInt(args[1]));

    cureClusterer.setShrinkFactor(SHRINK_FACTOR);

    cureClusterer.setRequiredRepresentationProbablity(
        REQUIRED_REPRESENTATION_PROBABILITY);

    cureClusterer.setNumberOfPartitions(NUMBER_OF_PARTITIONS);

    cureClusterer.setReducingFactorForEachPartition(
        REDUCING_FACTOR_FOR_EACH_PARTITION);

    ArrayList clusters = cureClusterer.cluster();

    // cureClusterer.showClusters(clusters);
  }
}
