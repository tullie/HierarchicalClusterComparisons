package chameleon;

import java.util.List;

public class Node {
  public int index;
  public double[] values;
  public List<Integer> neighbors;

  Node(int index, double[] values, List<Integer> neighbors) {
    this.index = index;
    this.values = values;
    this.neighbors = neighbors;
  }
}
