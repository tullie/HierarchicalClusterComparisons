package chameleonextension;

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
import java.util.Random;

import edu.wlu.cs.levy.CG.KDTree;
import edu.wlu.cs.levy.CG.KDException;

public class Chameleon {
  static final int K_NEAREST_NEIGHBOURS = 3;
  static final double ALPHA = 1;
  static final double REQUIRED_REPRESENTATION_PROBABLITY = 0.1;

  HMetisInterface hmetis = new HMetisInterface();
  Set<Integer> sampledIndexes = new HashSet<>();

  public List<List<Node>> cluster(List<double[]> dataset, int clusterCount) {
    if (dataset.isEmpty() || clusterCount <= 0) {
      System.out.println("Chameleon cluster given invalid parameters");
      return null;
    }

    int samplePoints = calculateSampleSize(clusterCount, dataset.size());
    System.out.println("Found sample size: " + samplePoints);
    List<double[]> sampledDataset = selectRandomPoints(samplePoints, dataset);

    List<Node> graph = constructGraph(dataset, K_NEAREST_NEIGHBOURS);

    System.out.println("Bisecting....");
    int minPartitonSize = Math.max((int)(dataset.size() * 0.025), 1);
    List<List<Node>> subClusters = bisectUntilSize(graph, minPartitonSize);

    System.out.println("Merging....");
    List<List<Node>> clusters = mergeSubClusters(subClusters, clusterCount);

    return labelRemainingDataPoints(clusters, dataset);
  }

  private int calculateSampleSize(int clusterCount, int datasetSize) {
    // Based off Chernoff Bounds mentioned in CURE algorithm.
    return (
        int)((0.5 * datasetSize) +
             (clusterCount *
              Math.log10(1 / REQUIRED_REPRESENTATION_PROBABLITY)) +
             (clusterCount *
              Math.sqrt(
                  Math.pow(Math.log10(1 / REQUIRED_REPRESENTATION_PROBABLITY),
                           2) +
                  (datasetSize / clusterCount) *
                      Math.log10(1 / REQUIRED_REPRESENTATION_PROBABLITY))));
  }

  private List<double[]> selectRandomPoints(int sampleSize,
                                            List<double[]> dataset) {
    List<double[]> datasetSample = new ArrayList<>();
    Random random = new Random();
    for (int i = 0; i < sampleSize; i++) {
      int index = random.nextInt(dataset.size());
      if (sampledIndexes.contains(index)) {
        i--;
        continue;
      }
      double[] dataPoints = dataset.get(index);
      datasetSample.add(dataPoints);
      sampledIndexes.add(index);
    }
    return datasetSample;
  }

  private List<List<Node>> labelRemainingDataPoints(List<List<Node>> clusters,
                                                    List<double[]> dataset) {
    Map<Integer, Integer> clusterFromNodeIndex = new HashMap<>();
    KDTree kdTree = new KDTree(clusters.get(0).get(0).values.length);
    for (int i = 0; i < clusters.size(); ++i) {
      for (Node node : clusters.get(i)) {
        clusterFromNodeIndex.put(node.originalIndex, i);
        try {
          kdTree.insert(node.values, node.originalIndex);
        } catch (KDException e) {
          String msg =
              "Exception when inserting into KDTree: " + e.getMessage();
          System.out.println(msg);
          continue;
        }
      }
    }

    for (int i = 0; i < dataset.size(); ++i) {
      if (sampledIndexes.contains(i)) {
        continue;
      }

      try {
        Integer nearestIndex = (Integer)kdTree.nearest(dataset.get(i));
        int clusterToAssign = clusterFromNodeIndex.get(nearestIndex);
        int newIndex = clusters.get(clusterToAssign).size();
        Node assignedNode = new Node(newIndex, dataset.get(i), null);
        clusters.get(clusterToAssign).add(assignedNode);
      } catch (Exception e) {
        String msg = "Exception when finding knn: " + e.getMessage();
        System.out.println(msg);
      }
    }

    return clusters;
  }

  private List<Node> constructGraph(List<double[]> dataset, int knn) {
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
      List<Node.Edge> neighbors = new ArrayList<>(knnIndexes.size());
      for (int neighborIndex : knnIndexes) {
        double dist =
            euclideanDistance(dataset.get(i), dataset.get(neighborIndex));
        neighbors.add(new Node.Edge(neighborIndex, dist));
      }
      graph.add(new Node(i, dataset.get(i), neighbors));
    }

    // Add directed edge back from neighbor nodes
    for (Node node : graph) {
      for (Node.Edge edge : node.neighbors) {
        Node neighbor = graph.get(edge.neighborIndex);

        boolean selfReferenceFound = false;
        for (Node.Edge neighborEdge : neighbor.neighbors) {
          if (neighborEdge.neighborIndex == node.index) {
            selfReferenceFound = true;
            break;
          }
        }

        if (!selfReferenceFound) {
          Node.Edge backEdge = new Node.Edge(node.index, edge.weight);
          neighbor.originalNeighbors.add(backEdge);
          neighbor.neighbors.add(backEdge);
        }
      }
    }

    return graph;
  }

  private double euclideanDistance(double[] lhs, double[] rhs) {
    double innerSum = 0;
    for (int i = 0; i < lhs.length; ++i) {
      innerSum += Math.pow(lhs[i] + rhs[i], 2);
    }
    return Math.sqrt(innerSum);
  }

  private List<List<Node>> bisectUntilSize(List<Node> graph, int cluster_size) {
    PriorityQueue<List<Node>> clusterSizeHeap =
        new PriorityQueue<>(graph.size(), new Comparator<List<Node>>() {
          @Override
          public int compare(List<Node> lhs, List<Node> rhs) {
            if (lhs.size() == rhs.size())
              return 0;
            if (lhs.size() > rhs.size())
              return -1;
            else
              return 1;
          }
        });

    List<Node> largestCluster = graph;
    while (largestCluster.size() > cluster_size) {
      if (largestCluster.size() <= 1) {
        break;
      }
      List<List<Node>> subClusters = hmetis.runMetisOnGraph(largestCluster, 2);

      System.out.println("Heap size: " + clusterSizeHeap.size());
      boolean emptyCluster = false;
      for (List<Node> subCluster : subClusters) {
        if (subCluster.size() == 0) {
          emptyCluster = true;
          continue;
        }
        subCluster = rebuildNodeIndexes(subCluster);
        clusterSizeHeap.add(subCluster);
      }

      if (emptyCluster) {
        break;
      }

      largestCluster = clusterSizeHeap.poll();
    }
    clusterSizeHeap.add(largestCluster);
    return new ArrayList<>(clusterSizeHeap);
  }

  private List<Node> rebuildNodeIndexes(List<Node> cluster) {
    Map<Integer, Integer> oldIndexToNew = new HashMap<>();
    for (int i = 0; i < cluster.size(); ++i) {
      oldIndexToNew.put(cluster.get(i).originalIndex, i);
      cluster.get(i).index = i;
    }

    for (Node node : cluster) {
      node.neighbors.clear();
      for (Node.Edge edge : node.originalNeighbors) {
        if (oldIndexToNew.containsKey(edge.neighborIndex)) {
          int newIndex = oldIndexToNew.get(edge.neighborIndex);
          node.neighbors.add(new Node.Edge(newIndex, edge.weight));
        }
      }
    }

    return cluster;
  }

  private List<List<Node>> mergeSubClusters(List<List<Node>> subClusters,
                                            int clusterCount) {
    PriorityQueue<double[]> clusterDistQueue =
        new PriorityQueue<>(subClusters.size(), new Comparator<double[]>() {
          @Override
          public int compare(double[] lhs, double[] rhs) {
            return Double.compare(rhs[0], lhs[0]);
          }
        });

    for (int i = 0; i < subClusters.size(); ++i) {
      for (int j = i + 1; j < subClusters.size(); ++j) {
        double[] distAndIndexes = new double[3];
        distAndIndexes[0] = distance(subClusters.get(i), subClusters.get(j));
        distAndIndexes[1] = (double)i;
        distAndIndexes[2] = (double)j;
        clusterDistQueue.add(distAndIndexes);
      }
    }

    Set<Integer> mergedIndexes = new HashSet<>();
    int currentClusterCount = subClusters.size();
    while (currentClusterCount > clusterCount) {
      double[] distAndIndexes = clusterDistQueue.poll();
      int lhsIndex = (int)distAndIndexes[1];
      int rhsIndex = (int)distAndIndexes[2];
      if (mergedIndexes.contains(lhsIndex) ||
          mergedIndexes.contains(rhsIndex)) {
        continue;
      }

      currentClusterCount--;
      System.out.println("Current cluster count: " + currentClusterCount);
      mergedIndexes.add(lhsIndex);
      mergedIndexes.add(rhsIndex);
      List<Node> lhsCluster = subClusters.get(lhsIndex);
      List<Node> rhsCluster = subClusters.get(rhsIndex);
      lhsCluster.addAll(rhsCluster);
      List<Node> mergedCluster = rebuildNodeIndexes(lhsCluster);

      for (int i = 0; i < subClusters.size(); ++i) {
        if (!mergedIndexes.contains(i)) {
          double[] mergedDistAndIndexes = new double[3];
          mergedDistAndIndexes[0] = distance(mergedCluster, subClusters.get(i));
          mergedDistAndIndexes[1] = (double)i;
          mergedDistAndIndexes[2] = (double)subClusters.size();
          clusterDistQueue.add(mergedDistAndIndexes);
        }
      }

      subClusters.add(mergedCluster);
    }

    List<List<Node>> resultClusters = new ArrayList<>();
    for (int i = 0; i < subClusters.size(); ++i) {
      if (!mergedIndexes.contains(i)) {
        resultClusters.add(subClusters.get(i));
      }
    }

    return resultClusters;
  }

  private double distance(List<Node> lhsCluster, List<Node> rhsCluster) {
    List<Node.Edge> lhsCrossedEdges = new ArrayList<>();
    if (lhsCluster.size() > 1) {
      List<List<Node>> lhsSplit = hmetis.runMetisOnGraph(lhsCluster, 2);
      lhsCrossedEdges = getCrossClusterEdges(lhsSplit.get(0), lhsSplit.get(0));
    }

    List<Node.Edge> rhsCrossedEdges = new ArrayList<>();
    if (rhsCluster.size() > 1) {
      List<List<Node>> rhsSplit = hmetis.runMetisOnGraph(rhsCluster, 2);
      rhsCrossedEdges = getCrossClusterEdges(rhsSplit.get(0), rhsSplit.get(0));
    }

    List<Node.Edge> crossedEdges = getCrossClusterEdges(lhsCluster, rhsCluster);
    List<Node.Edge> lhsInnerEdges =
        getCrossClusterEdges(lhsCluster, lhsCluster);
    List<Node.Edge> rhsInnerEdges =
        getCrossClusterEdges(rhsCluster, rhsCluster);

    double connectivity = relativeInterConnectivity(
        lhsCrossedEdges, rhsCrossedEdges, crossedEdges);
    System.out.println("Connectivity: " + connectivity);

    double closeness =
        relativeCloseness(lhsCrossedEdges, rhsCrossedEdges, crossedEdges,
                          lhsInnerEdges, rhsInnerEdges);
    System.out.println("Closeness: " + closeness);

    return closeness * Math.pow(closeness, ALPHA);
  }

  private double relativeInterConnectivity(List<Node.Edge> lhsCrossedEdges,
                                           List<Node.Edge> rhsCrossedEdges,
                                           List<Node.Edge> allCrossedEdges) {
    double mergedEdgeCutSum = totalWeightOfEdges(allCrossedEdges);
    double lhsEdgeCutSum = totalWeightOfEdges(lhsCrossedEdges);
    double rhsEdgeCutSum = totalWeightOfEdges(rhsCrossedEdges);

    if (lhsEdgeCutSum == 0 && rhsEdgeCutSum == 0) {
      return 0;
    } else {
      double result = mergedEdgeCutSum / ((lhsEdgeCutSum + rhsEdgeCutSum) / 2);
      return result;
    }
  }

  private double relativeCloseness(List<Node.Edge> lhsCrossedEdges,
                                   List<Node.Edge> rhsCrossedEdges,
                                   List<Node.Edge> allCrossedEdges,
                                   List<Node.Edge> lhsInnerEdges,
                                   List<Node.Edge> rhsInnerEdges) {
    double lhsEdgeCutAverage = averageWeightOfEdges(lhsCrossedEdges);
    double lhsWeightAverage = averageWeightOfEdges(lhsInnerEdges);
    double rhsEdgeCutAverage = averageWeightOfEdges(rhsCrossedEdges);
    double rhsWeightAverage = averageWeightOfEdges(rhsInnerEdges);
    double crossEdgeCut = averageWeightOfEdges(allCrossedEdges);
    double lhsRhsWeight = lhsWeightAverage + rhsWeightAverage;
    if (lhsRhsWeight == 0)
      return 0;
    double leftTerm = lhsWeightAverage / lhsRhsWeight * lhsEdgeCutAverage;
    double rightTerm = rhsWeightAverage / lhsRhsWeight * rhsEdgeCutAverage;
    if (leftTerm + rightTerm == 0)
      return 0;
    double result = crossEdgeCut / (leftTerm + rightTerm);
    return result;
  }

  private List<Node.Edge> getCrossClusterEdges(List<Node> lhsCluster,
                                               List<Node> rhsCluster) {
    Set<Integer> indexesInCluster = new HashSet<>();
    for (Node node : lhsCluster) {
      indexesInCluster.add(node.originalIndex);
    }

    List<Node.Edge> crossedEdges = new ArrayList<>();
    for (Node node : rhsCluster) {
      for (Node.Edge edge : node.originalNeighbors) {
        if (indexesInCluster.contains(edge.neighborIndex)) {
          crossedEdges.add(edge);
        }
      }
    }

    return crossedEdges;
  }

  private double averageWeightOfEdges(List<Node.Edge> edges) {
    if (edges.size() == 0) {
      return 0;
    }
    return totalWeightOfEdges(edges) / edges.size();
  }

  private static double totalWeightOfEdges(List<Node.Edge> edges) {
    double edgeWeightSum = 0;
    for (Node.Edge edge : edges) {
      edgeWeightSum += edge.weight;
    }
    return edgeWeightSum;
  }
}
