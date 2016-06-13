package chameleon;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Comparator;

import edu.wlu.cs.levy.CG.KDTree;
import edu.wlu.cs.levy.CG.KDException;

public class Chameleon {
  static final int K_NEAREST_NEIGHBOURS = 3;
  static final int MIN_PARTITION_SIZE = 10; // 1 - 5% of dataset size.

  public static void cluster(List<double[]> dataset, int clusterCount) {
    if (dataset.isEmpty() || clusterCount <= 0) {
      System.out.println("Chameleon cluster given invalid parameters");
      return;
    }

    List<Node> graph = constructGraph(dataset, K_NEAREST_NEIGHBOURS);
    List<List<Node>> subClusters = bisectUntilSize(graph, MIN_PARTITION_SIZE);

    int i = 1;
    for (List<Node> subCluster : subClusters) {
      System.out.println("Cluster " + i);
      i++;

      for (Node node : subCluster) {
        System.out.print(node.index + ", ");
      }
    }

    /*
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

  private static List<Node> constructGraph(List<double[]> dataset, int knn) {
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

    List<Node> graph = new ArrayList<>(dataset.size());
    for (int i = 0; i < dataset.size(); ++i) {
      List<Integer> knnIndexes;
      try {
        knnIndexes = kdTree.nearest(dataset.get(i), knn + 1);
        knnIndexes.remove(knnIndexes.size() - 1); // Remove self index.
      } catch (Exception e) {
        String msg = "Exception when finding knn: " + e.getMessage();
        System.out.println(msg);
        continue;
      }
      graph.add(new Node(i, dataset.get(i), knnIndexes));
    }

    return graph;
  }

  private static List<List<Node>> bisectUntilSize(List<Node> graph,
                                                  int cluster_size) {
    PriorityQueue<List<Node>> clusterSizeHeap =
        new PriorityQueue<>(graph.size(), new Comparator<List<Node>>() {
          @Override
          public int compare(List<Node> lhs, List<Node> rhs) {
            if (lhs.size() == rhs.size())
              return 0;
            if (lhs.size() > rhs.size())
              return 1;
            else
              return -1;
          }
        });

    List<Node> largestCluster = graph;
    while (largestCluster.size() > MIN_PARTITION_SIZE) {
      HMetisInterface hmetis = new HMetisInterface();
      List<List<Node>> subClusters = hmetis.runMetisOnGraph(largestCluster, 2);
      subClusters = correctIndexes(subClusters);
      for (List<Node> subCluster : subClusters) {
        clusterSizeHeap.add(subCluster);
      }
      largestCluster = clusterSizeHeap.poll();
      System.out.println(largestCluster.size());
    }
    return new ArrayList<>(clusterSizeHeap);
  }

  private static List<List<Node>> correctIndexes(List<List<Node>> subClusters) {
    for (List<Node> subCluster : subClusters) {
      Map<Integer, Integer> oldIndexToNew = new HashMap<>();
      for (int i = 0; i < subCluster.size(); ++i) {
        oldIndexToNew.put(subCluster.get(i).index, i);
        subCluster.get(i).index = i;
      }

      for (Node node : subCluster) {
        ListIterator<Integer> neighborIt = node.neighbors.listIterator();
        while (neighborIt.hasNext()) {
          Integer neighbor = neighborIt.next();
          if (oldIndexToNew.containsKey(neighbor)) {
            neighborIt.set(oldIndexToNew.get(neighbor));
          } else {
            neighborIt.remove();
          }
        }
      }
    }
    return subClusters;
  }
}
