package chameleon;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.LinkedList;

import edu.wlu.cs.levy.CG.KDTree;
import edu.wlu.cs.levy.CG.KDException;

public class Chameleon {
  static final int K_NEAREST_NEIGHBOURS = 3;

  public static class Node {
    int index;
    double[] values;
    List<Integer> knnIndexes;

    Node(int index, double[] values, List<Integer> knnIndexes) {
      this.index = index;
      this.values = values;
      this.knnIndexes = knnIndexes;
    }
  }

  public static void cluster(List<double[]> dataset, int clusterCount) {
    if (dataset.isEmpty() || clusterCount <= 0) {
      System.out.println("Chameleon cluster given invalid parameters");
      return;
    }

    Node[] graph = constructGraph(dataset);

    // Bisect graph
    /*
    // bisection = pref.get(BISECTION, "Kernighan-Lin");
    bisection = pref.get(BISECTION, "Fiduccia-Mattheyses");
    Bisection bisectionAlg =
    BisectionFactory.getInstance().getProvider(bisection);
    if (bisectionAlg instanceof FiducciaMattheyses) {
        FiducciaMattheyses fm = (FiducciaMattheyses) bisectionAlg;
        fm.setIterationLimit(pref.getInt(FiducciaMattheyses.ITERATIONS, 20));
    }

    partitioning = pref.get(PARTITIONING, "Recursive bisection");
    Partitioning partitioningAlg =
    PartitioningFactory.getInstance().getProvider(partitioning);
    partitioningAlg.setBisection(bisectionAlg);
    ArrayList<LinkedList<Node>> partitioningResult =
    partitioningAlg.partition(maxPartitionSize, g, pref);

    String merger = pref.get(MERGER, "pair merger");
    Merger m = MergerFactory.getInstance().getProvider(merger);
    m.setDistanceMeasure(knn.getDistanceMeasure());

    MergeEvaluationFactory mef = MergeEvaluationFactory.getInstance();
    if (m instanceof PairMerger) {
      similarityMeasure = pref.get(SIM_MEASURE, BBK1.name);
      MergeEvaluation me = mef.getProvider(similarityMeasure);
      ((PairMerger)m).setMergeEvaluation(me);
    } else if (m instanceof PairMergerMO) {
      PairMergerMO mo = (PairMergerMO)m;
      mo.clearObjectives();
      mo.addObjective(mef.getProvider(pref.get(OBJECTIVE_1)));
      mo.addObjective(mef.getProvider(pref.get(OBJECTIVE_2)));
    }
    noise = m.initialize(partitioningResult, g, bisectionAlg, pref, noise);
    HierarchicalResult result = m.getHierarchy(dataset, pref);
    result.setNoise(noise);
    return result;
    */
  }

  private static Node[] constructGraph(List<double[]> dataset) {
    int dimensions = dataset.get(0).length;
    KDTree kdTree = new KDTree(dimensions);
    for (int i = 0; i < dataset.size(); ++i) {
      try {
        kdTree.insert(dataset.get(i), i);
      } catch (KDException e) {
        String msg = "Exception when inserting into KDTree: " + e.getMessage();
        System.out.println(msg);
        continue;
      }
    }

    Node[] graph = new Node[dataset.size()];
    for (int i = 0; i < dataset.size(); ++i) {
      List<Integer> knnIndexes;
      try {
        knnIndexes = kdTree.nearest(dataset.get(i), K_NEAREST_NEIGHBOURS + 1);
        knnIndexes.remove(knnIndexes.size() - 1); // Remove self index.
      } catch (Exception e) {
        String msg = "Exception when finding knn: " + e.getMessage();
        System.out.println(msg);
        continue;
      }
      graph[i] = new Node(i, dataset.get(i), knnIndexes);
    }

    return graph;
  }
}
