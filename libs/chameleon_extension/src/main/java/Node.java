package chameleonextension;

import java.util.List;
import java.util.ArrayList;

public class Node {
  public static class Edge {
    int neighborIndex;
    double weight;

    public Edge(int neighborIndex, double weight) {
      this.neighborIndex = neighborIndex;
      this.weight = weight;
    }
  }

  public int index;
  public double[] values;
  public List<Edge> neighbors;
  public int originalIndex;
  public List<Edge> originalNeighbors;

  public Node(int index, double[] values, List<Edge> neighbors) {
    this.index = index;
    this.values = values;
    this.originalIndex = index;
    if (neighbors != null) {
      this.neighbors = neighbors;
      this.originalNeighbors = new ArrayList<>(neighbors);
    }
  }

  public Node(Node cpyNode) {
    this.index = cpyNode.index;
    this.values = cpyNode.values;
    this.neighbors = cpyNode.neighbors;
    this.originalIndex = cpyNode.originalIndex;
    this.originalNeighbors = cpyNode.originalNeighbors;
  }
}
