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
  static final int K_NEAREST_NEIGHBOURS = 7;
  static final double ALPHA = 2;
  static final int MIN_REPRESENTATION_COUNT = 200;
  static final int MIN_PARTITION_COUNT = 0;

  static final boolean RANDOM_SAMPLE = true;
  static final double REQUIRED_REPRESENTATION_PROBABLITY = 0.1;

  HMetisInterface hmetis = new HMetisInterface();
  Set<Integer> sampledIndexes = new HashSet<>();

  public List<List<Node>> cluster(List<double[]> dataset, int clusterCount) {
    if (dataset.isEmpty() || clusterCount <= 0) {
      System.out.println("Chameleon cluster given invalid parameters");
      return null;
    }

    List<Node> graph;
    if (RANDOM_SAMPLE) {
      int samplePoints = calculateSampleSize(clusterCount, dataset.size());
      System.out.println("Found sample size: " + samplePoints);
      List<double[]> sampleDataset = selectRandomPoints(samplePoints, dataset);
      graph = constructGraph(sampleDataset, K_NEAREST_NEIGHBOURS);
    } else {
      graph = constructGraph(dataset, K_NEAREST_NEIGHBOURS);
    }

    System.out.println("Bisecting....");
    int minPartitionSize = MIN_PARTITION_COUNT;
    if (minPartitionSize == 0) {
      minPartitionSize = Math.max((int)(dataset.size() * 0.025), 1);
    }
    List<List<Node>> bisectedGraph = bisectUntilSize(graph, minPartitionSize);

    System.out.println("Merging....");
    List<Cluster> subClusters = new ArrayList<>(bisectedGraph.size());
    for (List<Node> subCluster : bisectedGraph) {
      subClusters.add(new Cluster(subCluster));
    }
    List<Cluster> clusters = mergeSubClusters(subClusters, clusterCount);

    List<List<Node>> clusterNodes = new ArrayList<>(clusters.size());
    for (Cluster cluster : clusters) {
      clusterNodes.add(cluster.nodes);
    }

    if (RANDOM_SAMPLE) {
      return labelRemainingDataPoints(clusterNodes, dataset);
    } else {
      return clusterNodes;
    }
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

  private List<List<Node>> labelRemainingDataPoints(List<List<Node>> clusterNodes,
                                                    List<double[]> dataset) {
    Map<Integer, Integer> clusterFromNodeIndex = new HashMap<>();
    KDTree kdTree = new KDTree(clusterNodes.get(0).get(0).values.length);
    for (int i = 0; i < clusterNodes.size(); ++i) {
      for (Node node : clusterNodes.get(i)) {
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
        int newIndex = clusterNodes.get(clusterToAssign).size();
        Node assignedNode = new Node(newIndex, dataset.get(i), null);
        clusterNodes.get(clusterToAssign).add(assignedNode);
      } catch (KDException e) {
        String msg = "Exception when finding nearest: " + e.getMessage();
        System.out.println(msg);
      }
    }

    return clusterNodes;
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

  private List<Cluster> mergeSubClusters(List<Cluster> subClusters,
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
        distAndIndexes[0] = distance(subClusters.get(i).reprPoints,
                                     subClusters.get(j).reprPoints);
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
      Cluster lhsCluster = subClusters.get(lhsIndex);
      Cluster rhsCluster = subClusters.get(rhsIndex);
      Cluster mergedCluster = mergeClusters(lhsCluster, rhsCluster);

      for (int i = 0; i < subClusters.size(); ++i) {
        if (!mergedIndexes.contains(i)) {
          double[] mergedDistAndIndexes = new double[3];
          mergedDistAndIndexes[0] =
              distance(mergedCluster.reprPoints, subClusters.get(i).reprPoints);
          mergedDistAndIndexes[1] = (double)i;
          mergedDistAndIndexes[2] = (double)subClusters.size();
          clusterDistQueue.add(mergedDistAndIndexes);
        }
      }

      subClusters.add(mergedCluster);
    }

    List<Cluster> resultClusters = new ArrayList<>();
    for (int i = 0; i < subClusters.size(); ++i) {
      if (!mergedIndexes.contains(i)) {
        resultClusters.add(subClusters.get(i));
      }
    }

    return resultClusters;
  }

  private Cluster mergeClusters(Cluster lhsCluster, Cluster rhsCluster) {
    Cluster mergedCluster = lhsCluster;
    mergedCluster.reprPoints.addAll(rhsCluster.reprPoints);
    mergedCluster.nodes.addAll(rhsCluster.nodes);
    rebuildNodeIndexes(mergedCluster.reprPoints);
    rebuildNodeIndexes(mergedCluster.nodes);
    /*
    double[] mean = calculateClusterMean(mergedCluster.reprPoints);
    List<Node> reprNodes = new ArrayList();
    for (int i = 0; i < MIN_REPRESENTATION_COUNT && i <
    mergedCluster.reprPoints.size(); ++i) {
      double maxDist = 0;
      double minDist = 0;
      Node maxNode = null;
      for (int j = 0; j < mergedCluster.reprPoints.size(); ++j) {
        Node node = mergedCluster.reprPoints.get(j);
        if (j == 0) {
          minDist = euclideanDistance(node.values, mean);
        } else {
          minDist = calculateMinDistanceFromSet(node, reprNodes);
        }
        if (minDist >= maxDist) {
          maxDist = minDist;
          maxNode = node;
        }
      }
      reprNodes.add(maxNode);
    }

    mergedCluster.reprPoints = reprNodes;
    */
    return mergedCluster;
  }

  private double calculateMinDistanceFromSet(Node node, List<Node> set) {
    double minDistance = Double.MIN_VALUE;
    for (Node setNode : set) {
      if (setNode.originalIndex == node.originalIndex) {
        continue;
      }

      double distance = euclideanDistance(node.values, setNode.values);
      if (minDistance > distance) {
        minDistance = distance;
      }
    }

    if (minDistance == Double.MIN_VALUE) {
      return 0;
    } else
      return minDistance;
  }

  private double[] calculateClusterMean(List<Node> cluster) {
    int dimension = cluster.get(0).values.length;
    double[] mean = new double[dimension];
    for (Node node : cluster) {
      for (int i = 0; i < mean.length; ++i) {
        mean[i] += node.values[i];
      }
    }

    for (int i = 0; i < mean.length; ++i) {
      mean[i] /= cluster.size();
    }

    return mean;
  }

  private double distance(List<Node> lhsCluster, List<Node> rhsCluster) {
    List<Node.Edge> lhsCrossedEdges = new ArrayList<>();
    if (lhsCluster.size() > 1) {
      List<List<Node>> lhsSplit = hmetis.runMetisOnGraph(lhsCluster, 2);
      lhsCrossedEdges = getCrossClusterEdges(lhsSplit.get(0), lhsSplit.get(1));
    }

    List<Node.Edge> rhsCrossedEdges = new ArrayList<>();
    if (rhsCluster.size() > 1) {
      List<List<Node>> rhsSplit = hmetis.runMetisOnGraph(rhsCluster, 2);
      rhsCrossedEdges = getCrossClusterEdges(rhsSplit.get(0), rhsSplit.get(1));
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

    if (lhsEdgeCutSum + rhsEdgeCutSum == 0) {
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
