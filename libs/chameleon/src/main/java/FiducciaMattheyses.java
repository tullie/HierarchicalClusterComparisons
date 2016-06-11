package chameleon;

import java.util.List;
import java.util.ArrayList;

public class FiducciaMattheyses {
  static final int ITERATION_LIMIT = 20;

  int maxDegree = 0;
  List<Node> graph;
  Integer[] differenceBuckets;
  Integer[] swapHistory;
  Integer[] swapHistoryCost;
  Integer[] nodeNext;
  Integer[] nodePrevious;
  Integer[] nodeDifference;
  boolean[] clusterForNode;

  public FiducciaMattheyses(List<Node> graph) { this.graph = graph; }

  public List<List<Node>> bisect() {

    // Find max Node.
    for (Node node : graph) {
      maxDegree = Math.max(node.neighbors.size(), maxDegree);
    }

    differenceBuckets = new Integer[2 * maxDegree + 1];
    nodeNext = new Integer[graph.size()];
    nodePrevious = new Integer[graph.size()];
    nodeDifference = new Integer[graph.size()];
    swapHistory = new Integer[graph.size()];
    swapHistoryCost = new Integer[graph.size()];
    clusterForNode = new boolean[graph.size()];

    // Initial node partition.
    for (int i = graph.size() / 2; i < graph.size(); ++i) {
      clusterForNode[i] = true;
    }

    for (int i = 0; i < ITERATION_LIMIT; ++i) {
      boolean[] usedNodes = new boolean[graph.size()];

      // Prepare difference buckets.
      for (int j = 0; j <= maxDegree; ++j) {
        differenceBuckets[j] = null;
      }

      // Compute differences.
      for (int j = 0; j < graph.size(); ++j) {
        int difference = 0;
        for (int neighbor : graph.get(j).neighbors) {
          if (clusterForNode[j] == clusterForNode[neighbor]) {
            difference--;
          } else {
            difference++;
          }
        }

        // Add to bucket.
        addIntoBucket(difference, j);
        nodeDifference[j] = difference;
      }

      System.out.println("Finished computing differences");

      // Compute costs.
      for (int j = 0; j < graph.size(); ++j) {
        int bestNodeIndex = findBestNodeToSwapFromCluster(j % 2 == 0);

        // Add to history.
        swapHistory[j] = bestNodeIndex;
        swapHistoryCost[j] = nodeDifference[bestNodeIndex];
        usedNodes[bestNodeIndex] = true;

        System.out.println("Best index: " + bestNodeIndex);

        // Update differences.
        for (int neighbor : graph.get(bestNodeIndex).neighbors) {
          if (usedNodes[neighbor]) {
            continue;
          }

          removeFromBucket(neighbor);
          if (clusterForNode[neighbor] == clusterForNode[bestNodeIndex]) {
            nodeDifference[neighbor] += 2;
          } else {
            nodeDifference[neighbor] -= 2;
          }
          addIntoBucket(nodeDifference[neighbor], neighbor);
        }
        removeFromBucket(bestNodeIndex);
      }

      System.out.println("Finished computing costs");

      // Swap up to best index.
      int maxDifference = 0;
      int differenceSum = 0;
      int maxDifferenceIndex = -1;
      for (int j = 0; j < graph.size(); j++) {
        differenceSum += swapHistoryCost[j];

        // Keep the bisection balanced, only swap pairs
        if (differenceSum > maxDifference && j % 2 == 0) {
          maxDifference = differenceSum;
          maxDifferenceIndex = j;
        }
      }

      if (maxDifferenceIndex == -1) {
        break;
      }

      for (int j = 0; j <= maxDifferenceIndex; j++) {
        int swapHist = swapHistory[j];
        clusterForNode[swapHist] = !clusterForNode[swapHist];
      }
    }

    // Construct clusters from label.
    List<Node> leftCluster = new ArrayList<>(graph.size() / 2);
    List<Node> rightCluster = new ArrayList<>(graph.size() / 2);
    for (int i = 0; i < graph.size(); ++i) {
      if (clusterForNode[i]) {
        leftCluster.add(graph.get(i));
      } else {
        rightCluster.add(graph.get(i));
      }

    }
    List<List<Node>> subClusters = new ArrayList<>();
    subClusters.add(leftCluster);
    subClusters.add(rightCluster);
    return subClusters;
  }

  private int findBestNodeToSwapFromCluster(boolean cluster) {
    for (int j = 0; j < graph.size(); ++j) {
      int k = maxDegree * 2;
      while (k >= 0) {
        Integer nodeIndex = differenceBuckets[k];
        while (nodeIndex != null) {
          if (clusterForNode[nodeIndex] == cluster) {
            return nodeIndex;
          }
          nodeIndex = nodeNext[nodeIndex];
        }
        k--;
      }
    }
    return findBestNodeToSwapFromCluster(!cluster);
  }

  private void removeFromBucket(int node) {
    if (nodePrevious[node] != null) {
      nodeNext[nodePrevious[node]] = nodeNext[node];
    } else {
      differenceBuckets[nodeDifference[node] + maxDegree] = nodeNext[node];
    }

    if (nodeNext[node] != null) {
      nodePrevious[nodeNext[node]] = nodePrevious[node];
    }
  }

  private void addIntoBucket(int difference, int node) {
    int bucketPos = difference + maxDegree;
    System.out.println("BucketPos: " + bucketPos);
    if (differenceBuckets[bucketPos] == null) {
      differenceBuckets[bucketPos] = node;
      nodeNext[node] = null;
    } else {
      nodePrevious[differenceBuckets[bucketPos]] = node;
      nodeNext[node] = differenceBuckets[bucketPos];
      differenceBuckets[bucketPos] = node;
    }
    nodePrevious[node] = null;
  }
}
