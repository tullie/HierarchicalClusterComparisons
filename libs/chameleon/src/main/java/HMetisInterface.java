package chameleon;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.io.IOException;
import java.net.URL;
import java.io.UnsupportedEncodingException;
import java.io.PrintWriter;
import java.util.List;
import java.util.ArrayList;
import java.lang.Runtime;
import java.lang.ProcessBuilder;
import java.util.Enumeration;

public class HMetisInterface {
  static final double UFACTOR = 5.0; // Balance factor (%).

  static int current = 0;

  File metisFile = null;

  public List<List<Node>> runMetisOnGraph(List<Node> graph, int k) {

    // Format graph into input file.
    File inputFile = new File("inputGraph-" + current++);
    String filename = inputFile.getName();
    hMetisExport(graph, inputFile);

    if (metisFile == null) {
      metisFile = getBinary();
    }

    // Run hMetis.
    String space = " ";
    StringBuilder sb = new StringBuilder(metisFile.getAbsolutePath());
    sb.append(space)
        .append(filename)
        .append(space)
        .append(String.valueOf(k))
        .append(space)
        .append(String.valueOf(UFACTOR))
        .append(space);

    try {
      Process p = Runtime.getRuntime().exec(sb.toString());
      // System.out.println(readStdout(p));
      // System.out.println(readStderr(p));
      readStdout(p);
      p.waitFor();
    } catch (IOException | InterruptedException e) {
      e.printStackTrace();
      return null;
    }

    inputFile.delete();

    List<List<Node>> subclusters = importMetrisResult(filename, k, graph);
    return subclusters;
  }

  private List<List<Node>> importMetrisResult(String path, int k,
                                              List<Node> graph) {
    List<List<Node>> result = new ArrayList<>();
    for (int i = 0; i < k; i++) {
      result.add(new ArrayList<Node>());
    }

    File outputFile = new File(path + ".part." + String.valueOf(k));
    try (BufferedReader br = new BufferedReader(new FileReader(outputFile))) {
      String line;
      int i = 0;
      while ((line = br.readLine()) != null) {
        int partitionNumber = Integer.parseInt(line);
        result.get(partitionNumber).add(graph.get(i));
        i++;
      }
      outputFile.delete();
    } catch (IOException e) {
      e.printStackTrace();
    }

    return result;
  }

  private void cleanup() {
    if (metisFile != null && metisFile.exists()) {
      System.out.println("deleting " + metisFile.getAbsolutePath());
      metisFile.deleteOnExit();
    }
  }

  private File getBinary() {
    ClassLoader classLoader = Thread.currentThread().getContextClassLoader();

    try {
      Enumeration<URL> roots = classLoader.getResources(".");
      while (roots.hasMoreElements()) {
        URL url = roots.nextElement();
        URL urlWithHmetisPath =
            new URL(url.getProtocol(), url.getHost(), url.getPort(),
                    url.getFile() + "shmetis", null);
        File f = new File(urlWithHmetisPath.toURI());
        if (f.exists()) {
          Process p =
              Runtime.getRuntime().exec("chmod ugo+x " + f.getAbsolutePath());
          readStdout(p);
          readStderr(p);
          return f;
        }
      }
    } catch (IOException | URISyntaxException e) {
      e.printStackTrace();
    }

    return null;
  }

  private void hMetisExport(List<Node> graph, File target) {
    StringBuilder sb;
    try (PrintWriter writer = new PrintWriter(target, "UTF-8")) {

      sb = new StringBuilder();

      // A hyperedge is formed by node's neighbourhood.
      sb.append(graph.size()).append(" ").append(graph.size()).append("\n");
      String space = " ";
      for (int i = 0; i < graph.size(); i++) {

        // Append self.
        sb.append(i + 1);
        for (Node.Edge edge : graph.get(i).neighbors) {
          sb.append(space).append(edge.neighborIndex + 1);
        }
        sb.append("\n");
      }
      writer.write(sb.toString());
    } catch (UnsupportedEncodingException | FileNotFoundException e) {
      e.printStackTrace();
    }
  }

  public String readStdout(Process p) {
    StringBuilder sb = new StringBuilder();
    try {
      BufferedReader input =
          new BufferedReader(new InputStreamReader(p.getInputStream()));
      String line;
      while ((line = input.readLine()) != null) {
        sb.append(line).append("\n");
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
    return sb.toString();
  }

  public String readStderr(Process p) {
    StringBuilder sb = new StringBuilder();
    try {
      BufferedReader input =
          new BufferedReader(new InputStreamReader(p.getErrorStream()));
      String line;
      while ((line = input.readLine()) != null) {
        sb.append(line).append("\n");
      }
    } catch (IOException e) {
      e.printStackTrace();
    }
    return sb.toString();
  }
}
