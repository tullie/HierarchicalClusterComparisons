package chameleon;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.util.HashSet;
import java.util.Set;
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
  static final double ALPHA = 1;

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

    List<List<Node>> clusters = mergeSubClusters(subClusters, clusterCount);

    int j = 1;
    for (List<Node> subCluster : clusters) {
      System.out.println("Cluster " + j);
      j++;

      for (Node node : subCluster) {
        System.out.print(node.index + ", ");
      }
    }
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

      for (List<Node> subCluster : subClusters) {
        subCluster = rebuildNodeIndexes(subCluster);
        clusterSizeHeap.add(subCluster);
      }
      largestCluster = clusterSizeHeap.poll();
    }
    return new ArrayList<>(clusterSizeHeap);
  }

  private static List<Node> rebuildNodeIndexes(List<Node> cluster) {
    Map<Integer, Integer> oldIndexToNew = new HashMap<>();
    for (int i = 0; i < cluster.size(); ++i) {
      oldIndexToNew.put(cluster.get(i).originalIndex, i);
      cluster.get(i).index = i;
    }

    for (Node node : cluster) {
      node.neighbors.clear();
      for (Integer originalNeighbor : node.originalNeighbors) {
        if (oldIndexToNew.containsKey(originalNeighbor)) {
          node.neighbors.add(oldIndexToNew.get(originalNeighbor));
        }
      }
    }

    return cluster;
  }

  private static List<List<Node>> mergeSubClusters(List<List<Node>> subClusters,
                                                   int clusterCount) {
    PriorityQueue<double[]> clusterDistQueue =
        new PriorityQueue<>(subClusters.size(), new Comparator<double[]>() {
          @Override
          public int compare(double[] lhs, double[] rhs) {
            return Double.compare(lhs[0], rhs[0]);
          }
        });

    Map<Integer, List<double[]>> pqValueByCluster = new HashMap<>();
    for (int i = 0; i < subClusters.size(); ++i) {
      for (int j = i + 1; j < subClusters.size(); ++j) {
        double[] distAndIndexes = new double[3];
        distAndIndexes[0] = distance(subClusters.get(i), subClusters.get(j));
        distAndIndexes[1] = (double)i;
        distAndIndexes[2] = (double)j;
        clusterDistQueue.add(distAndIndexes);
        addToListMap(pqValueByCluster, i, distAndIndexes);
        addToListMap(pqValueByCluster, j, distAndIndexes);
      }
    }

    Set<Integer> mergedIndexes = new HashSet<>();
    int currentClusterCount = subClusters.size();
    while (currentClusterCount > clusterCount) {
      currentClusterCount--;
      double[] distAndIndexes = clusterDistQueue.poll();
      int lhsIndex = (int)distAndIndexes[1];
      int rhsIndex = (int)distAndIndexes[2];
      mergedIndexes.add(lhsIndex);
      mergedIndexes.add(rhsIndex);

      for (double[] pqValue : pqValueByCluster.get(lhsIndex)) {
        pqValueByCluster.remove(pqValue);
      }
      for (double[] pqValue : pqValueByCluster.get(rhsIndex)) {
        pqValueByCluster.remove(pqValue);
      }

      List<Node> lhsCluster = subClusters.get(lhsIndex);
      List<Node> rhsCluster = subClusters.get(rhsIndex);

      List<Node> mergedCluster = mergeClusters(lhsCluster, rhsCluster);

      for (int i = 0; i < subClusters.size(); ++i) {
        double[] mergedDistAndIndexes = new double[3];
        mergedDistAndIndexes[0] = distance(mergedCluster, subClusters.get(i));
        mergedDistAndIndexes[1] = (double)i;
        mergedDistAndIndexes[2] = (double)subClusters.size();
        clusterDistQueue.add(distAndIndexes);
        addToListMap(pqValueByCluster, i, distAndIndexes);
        addToListMap(pqValueByCluster, subClusters.size(), distAndIndexes);
      }

      subClusters.add(mergedCluster);
    }

    List<List<Node>> resultClusters =  new ArrayList<>();
    for (int i = 0; i < subClusters.size(); ++i) {
      if (!mergedIndexes.contains(i)) {
        resultClusters.add(subClusters.get(i));
      }
    }

    return resultClusters;
  }

  private static double distance(List<Node> lhsCluster, List<Node> rhsCluster) {
    double connectivity = relativeInterConnectivity(lhsCluster, rhsCluster);
    double closeness = relativeCloseness(lhsCluster, rhsCluster);
    return connectivity * Math.pow(closeness, ALPHA);
  }

  private static int relativeInterConnectivity(List<Node> lhsCluster,
                                               List<Node> rhsCluster) {
    List<Node> mergedCluster = mergeClusters(lhsCluster, rhsCluster);
    int mergedEdgeCut = edgeCutWeight(mergedCluster);
    int lhsEdgeCut = edgeCutWeight(lhsCluster);
    int rhsEdgeCut = edgeCutWeight(rhsCluster);
    return mergedEdgeCut / ((lhsEdgeCut + rhsEdgeCut) / 2);
  }

  private static int relativeCloseness(List<Node> lhsCluster,
                                       List<Node> rhsCluster) {
    return 0;
  }

  private static int edgeCutWeight(List<Node> cluster) { return 0; }

  private static List<Node> mergeClusters(List<Node> lhsCluster,
                                          List<Node> rhsCluster) {
    lhsCluster.addAll(rhsCluster);
    return rebuildNodeIndexes(lhsCluster);
  }

  private static <E, K> void addToListMap(Map<E, List<K>> map, E key, K value) {
    if (!map.containsKey(key)) {
      map.put(key, new ArrayList<>());
    }

    List<K> valueList = map.get(key);
    valueList.add(value);
    map.put(key, valueList);
  }
}
