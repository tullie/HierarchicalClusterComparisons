import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

import edu.gatech.gtisc.jbirch.cftree.CFTree;

class BirchDriver {
  static final int MAX_NODE_ENTRIES = 100;
  static final double DIST_THRESHOLD = 1;
  static final int DIST_FUNCTION = CFTree.D0_DIST;
  static final boolean APPLY_MERGE_REFINEMENT = true;
  static final int MEMORY_LIMIT = 1024; // MB.
  static final boolean AUTO_REBUILD = true;
  static final int MEMORY_LIMIT_PERIODIC_CHECK = 10000; // iterations.

  public static void main(String args[]) throws IOException {
    if (args.length != 1) {
      System.out.println("Dataset filename required");
      return;
    }

    System.out.println("Running BIRCH driver...");

    CFTree birchTree = new CFTree(MAX_NODE_ENTRIES, DIST_THRESHOLD,
                                  DIST_FUNCTION, APPLY_MERGE_REFINEMENT);
    birchTree.setAutomaticRebuild(AUTO_REBUILD);
    birchTree.setMemoryLimitMB(MEMORY_LIMIT);
    birchTree.setPeriodicMemLimitCheck(MEMORY_LIMIT_PERIODIC_CHECK);

    List<double[]> dataset = Utils.readDataset(args[0]);

    for (int i = 0; i < dataset.size(); ++i) {
      boolean didInsert = birchTree.insertEntry(dataset.get(i));
      if (!didInsert) {
        throw new IOException("Value on line " + i + " cannot be inserted");
      }
    }
    birchTree.finishedInsertingData();

    // Print indexes of all values for each subcluster.
    // TODO: Save these in files.
    List<ArrayList<Integer>> subclusters = birchTree.getSubclusterMembers();
    for (List<Integer> subcluster : subclusters) {
      System.out.println(Arrays.toString(subcluster.toArray(new Integer[0])));
    }
    System.out.println("Total CF-Nodes: " + birchTree.countNodes());
    System.out.println("Total CF-Entries: " + birchTree.countEntries());
    System.out.println("Total CF-Leaves: " + birchTree.countLeafEntries());
  }
}
