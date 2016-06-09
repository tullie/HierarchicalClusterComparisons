package cure;

import java.io.*;
import java.util.*;

public class Cure {
  private int totalNumberOfPoints;
  private int numberOfClusters;
  private int minRepresentativeCount;
  private double shrinkFactor;
  private double requiredRepresentationProbablity;
  private int numberOfPartitions;
  private int reducingFactorForEachPartition;

  private Point[] dataPoints;
  private ArrayList outliers;
  private HashMap dataPointsMap;

  private static int currentRepAdditionCount;
  private Hashtable integerTable;

  public Cure(List<double[]> dataset, int clusterCount) {
    initializeParameters(dataset, clusterCount);
  }

  public void setShrinkFactor(double shrinkFactor) {
    this.shrinkFactor = shrinkFactor;
  }

  public void
  setRequiredRepresentationProbablity(double requiredRepresentationProbablity) {
    this.requiredRepresentationProbablity = requiredRepresentationProbablity;
  }

  public void setNumberOfPartitions(int numberOfPartitions) {
    this.numberOfPartitions = numberOfPartitions;
  }

  public void
  setReducingFactorForEachPartition(int reducingFactorForEachPartition) {
    this.reducingFactorForEachPartition = reducingFactorForEachPartition;
  }

  public ArrayList cluster() {
    KDClusterNode KDNode = new KDClusterNode(totalNumberOfPoints);
    KDNode.pointsInCluster = totalNumberOfPoints;
    KDClusterNode.noOfDimension = 2;
    for (int j = 0; j < KDNode.pointsInCluster; j++) {
      KDNode.kdpoints[j] = new PointsNDim(KDClusterNode.noOfDimension);
    }

    for (int j = 0; j < KDNode.pointsInCluster; j++) {
      PointsNDim p = (PointsNDim)KDNode.kdpoints[j];

      p.coord[0] = dataPoints[j].x;
      p.coord[1] = dataPoints[j].y;
    }
    KDNode.parent = null;
    KDNode.left = null;
    KDNode.right = null;

    KDClusterNode.normalizePoints(KDNode);
    KDNode.calculateMedianMinMax();
    KDNode.calculateDensity();
    if (KDNode.pointsInCluster > 10) {
      KDNode.isLeaf = false;
      KDClusterNode.buildTree(KDNode, 0);
    } else {
      KDNode.isLeaf = true;
    }
    KDClusterNode.denormalizePoints(KDNode);
    KDClusterNode.clustercount = 0;
    KDClusterNode.inOrder(KDNode);
    numberOfClusters = KDClusterNode.clustercount;
    long beginTime = System.currentTimeMillis();

    int sampleSize = calculateSampleSize();
    if (sampleSize >= totalNumberOfPoints) {
      sampleSize = totalNumberOfPoints;
    }

    ArrayList randomPointSet = selectRandomPoints(sampleSize);
    ArrayList[] partitionedPointSet = partitionPointSet(randomPointSet);
    ArrayList subClusters = clusterSubPartitions(partitionedPointSet);

    if (reducingFactorForEachPartition >= 10) {
      eliminateOutliersFirstStage(subClusters, 1);
    } else {
      eliminateOutliersFirstStage(subClusters, 0);
    }
    ArrayList clusters = clusterAll(subClusters);
    clusters = labelRemainingDataPoints(clusters);

    long time = System.currentTimeMillis() - beginTime;

    System.out.println("The Algorithm took " + time +
                       " milliseconds to complete.");

    return clusters;
  }


  private void initializeParameters(List<double[]> dataset, int clusterCount) {
    totalNumberOfPoints = dataset.size();
    numberOfClusters = clusterCount;
    minRepresentativeCount = 4;
    shrinkFactor = 0.5;
    requiredRepresentationProbablity = 0.1;
    numberOfPartitions = clusterCount;
    reducingFactorForEachPartition = 2;
    dataPoints = new Point[totalNumberOfPoints];
    dataPointsMap = new HashMap();
    currentRepAdditionCount = totalNumberOfPoints;
    integerTable = new Hashtable();
    outliers = new ArrayList();

    for (int i = 0; i < dataset.size(); ++i) {
      if (dataset.get(i).length != 2) {
        System.out.println("Element " + i + " has incorrect dimension");
        continue;
      }
      double x = dataset.get(i)[0];
      double y = dataset.get(i)[1];
      dataPoints[i] = new Point(x, y, i);
      dataPointsMap.put(i, dataPoints[i]);
    }
  }

  /**
   * Calculates the Sample Size based on Chernoff Bounds Mentioned in the CURE
   * Algorithm
   */
  private int calculateSampleSize() {
    return (int)((0.5 * totalNumberOfPoints) +
                 (numberOfClusters *
                  Math.log10(1 / requiredRepresentationProbablity)) +
                 (numberOfClusters *
                  Math.sqrt(
                      Math.pow(Math.log10(1 / requiredRepresentationProbablity),
                               2) +
                      (totalNumberOfPoints / numberOfClusters) *
                          Math.log10(1 / requiredRepresentationProbablity))));
  }

  /**
   * Select random points from the data set
   */
  private ArrayList selectRandomPoints(int sampleSize) {
    ArrayList randomPointSet = new ArrayList();
    Random random = new Random();
    for (int i = 0; i < sampleSize; i++) {
      int index = random.nextInt(totalNumberOfPoints);
      if (integerTable.containsKey(index)) {
        i--;
        continue;
      } else {
        Point point = dataPoints[index];
        randomPointSet.add(point);
        integerTable.put(index, "");
      }
    }
    return randomPointSet;
  }

  /**
   * Partition the sampled data points to p partitions (p = numberOfPartitions)
   */
  private ArrayList[] partitionPointSet(ArrayList pointSet) {
    ArrayList partitionedSet[] = new ArrayList[numberOfPartitions];
    Iterator iter = pointSet.iterator();
    for (int i = 0; i < numberOfPartitions - 1; i++) {
      partitionedSet[i] = new ArrayList();
      int pointIndex = 0;
      while (pointIndex < pointSet.size() / numberOfPartitions) {
        partitionedSet[i].add(iter.next());
        pointIndex++;
      }
    }
    partitionedSet[numberOfPartitions - 1] = new ArrayList();
    while (iter.hasNext()) {
      partitionedSet[numberOfPartitions - 1].add(iter.next());
    }
    return partitionedSet;
  }

  /**
   * Cluster each partitioned set to n/pq clusters
   */
  private ArrayList clusterSubPartitions(ArrayList partitionedSet[]) {
    ArrayList clusters = new ArrayList();
    int numberOfClusterInEachPartition =
        totalNumberOfPoints /
        (numberOfPartitions * reducingFactorForEachPartition);
    for (int i = 0; i < partitionedSet.length; i++) {
      ClusterSet clusterSet =
          new ClusterSet(partitionedSet[i], numberOfClusterInEachPartition,
                         minRepresentativeCount, shrinkFactor, dataPointsMap);
      Cluster[] subClusters = clusterSet.getAllClusters();
      for (int j = 0; j < subClusters.length; j++) {
        clusters.add(subClusters[j]);
      }
    }
    return clusters;
  }

  /**
   * Eliminates outliers after pre-clustering
   */
  private void eliminateOutliersFirstStage(ArrayList clusters,
                                           int outlierEligibilityCount) {
    Iterator iter = clusters.iterator();
    ArrayList clustersForRemoval = new ArrayList();
    while (iter.hasNext()) {
      Cluster cluster = (Cluster)iter.next();
      if (cluster.getClusterSize() <= outlierEligibilityCount) {
        updateOutlierSet(cluster);
        clustersForRemoval.add(cluster);
      }
    }
    while (!clustersForRemoval.isEmpty()) {
      Cluster c = (Cluster)clustersForRemoval.remove(0);
      clusters.remove(c);
    }
  }

  /**
   * Cluster all remaining clusters. Merge all clusters using CURE's
   * hierarchical clustering algorithm till specified number of clusters
   * remain.
   */
  private ArrayList clusterAll(ArrayList clusters) {
    ClusterSet clusterSet =
        new ClusterSet(clusters, numberOfClusters, minRepresentativeCount,
                       shrinkFactor, dataPointsMap, true);
    return clusterSet.mergeClusters();
  }

  /**
   * Assign all remaining data points which were not part of the sampled data
   * set to set of clusters formed
   */
  private ArrayList labelRemainingDataPoints(ArrayList clusters) {

    for (int index = 0; index < dataPoints.length; index++) {
      if (integerTable.containsKey(index))
        continue;
      Point p = dataPoints[index];
      double smallestDistance = 1000000;
      int nearestClusterIndex = -1;
      for (int i = 0; i < clusters.size(); i++) {
        ArrayList rep = ((Cluster)clusters.get(i)).rep;
        for (int j = 0; j < rep.size(); j++) {
          double distance = p.calcDistanceFromPoint((Point)rep.get(j));
          if (distance < smallestDistance) {
            smallestDistance = distance;
            nearestClusterIndex = i;
          }
        }
      }
      if (nearestClusterIndex != -1) {
        ((Cluster)clusters.get(nearestClusterIndex)).pointsInCluster.add(p);
      }
    }
    return clusters;
  }

  /**
   * Update the outlier set for the clusters which have been identified as
   * outliers
   */
  private void updateOutlierSet(Cluster cluster) {
    ArrayList outlierPoints = cluster.getPointsInCluster();
    Iterator iter = outlierPoints.iterator();
    while (iter.hasNext()) {
      outliers.add(iter.next());
    }
  }

  private void debug(Exception e) { e.printStackTrace(System.out); }

  /**
   * Gets the current representative count so that the new points added do not
   * conflict with older KD Tree indices
   */
  public static int getCurrentRepCount() { return ++currentRepAdditionCount; }

  public void showClusters(ArrayList clusters) {
    for (int i = 0; i < numberOfClusters; i++) {
      Cluster cluster = (Cluster)clusters.get(i);
      logCluster(cluster, "cluster" + i);
    }
    logOutlier();
    logPlotScript(clusters.size());
  }

  private BufferedWriter getWriterHandle(String filename) {
    BufferedWriter out = null;
    try {
      FileWriter fw = new FileWriter(filename, true);
      out = new BufferedWriter(fw);
    } catch (Exception e) {
      debug(e);
    }
    return out;
  }

  private void closeWriterHandle(BufferedWriter out) {
    try {
      out.flush();
      out.close();
    } catch (Exception e) {
      debug(e);
    }
  }

  private void logCluster(Cluster cluster, String filename) {
    BufferedWriter out = getWriterHandle(filename);
    try {
      out.write("#\tX\tY\n");
      for (int j = 0; j < cluster.pointsInCluster.size(); j++) {
        Point p = (Point)cluster.pointsInCluster.get(j);
        out.write("\t" + p.x + "\t" + p.y + "\n");
      }
    } catch (Exception e) {
      debug(e);
    }
    closeWriterHandle(out);
  }

  private void logOutlier() {
    BufferedWriter out = getWriterHandle("outliers");
    try {
      out.write("#\tX\tY\n");
      for (int j = 0; j < outliers.size(); j++) {
        Point p = (Point)outliers.get(j);
        out.write("\t" + p.x + "\t" + p.y + "\n");
      }
    } catch (Exception e) {
      debug(e);
    }
    closeWriterHandle(out);
  }

  private void logPlotScript(int totalClusters) {
    BufferedWriter out = getWriterHandle("plotcure.txt");
    try {
      setPlotStyle(out);
      out.write("plot");
      for (int i = 0; i < totalClusters; i++) {
        out.write(" \"cluster" + i + "\",");
      }
      out.write(" \"outliers\"");
    } catch (Exception e) {
      debug(e);
    }
    closeWriterHandle(out);
  }

  private void setPlotStyle(BufferedWriter out) {
    try {
      out.write("reset\n");
      out.write("set size ratio 2\n");
      out.write("unset key\n");
      out.write("set title \"CURE\"\n");
    } catch (Exception e) {
      debug(e);
    }
  }
}
