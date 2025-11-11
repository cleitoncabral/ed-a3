
import java.io.*;
import java.nio.file.*;
import java.util.*;

class Node {

    public final int index;
    public double g;
    public double h;
    public double f;
    public Node parent;

    public Node(int index) {
        this.index = index;
        this.g = 0;
        this.h = 0;
        this.f = 0;
        this.parent = null;
    }

    public Node(int index, double cost) {
        this.index = index;
        this.g = cost;
        this.h = 0;
        this.f = 0;
        this.parent = null;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null || getClass() != obj.getClass()) {
            return false;
        }
        Node node = (Node) obj;
        return index == node.index;
    }

    @Override
    public int hashCode() {
        return Objects.hash(index);
    }

    @Override
    public String toString() {
        return String.valueOf(index);
    }
}

class SearchResult {

    public List<Node> path;
    public double totalCost;
    public int nodesExpanded;
    public long executionTime;

    public SearchResult() {
        this.path = new ArrayList<>();
        this.totalCost = 0;
        this.nodesExpanded = 0;
        this.executionTime = 0;
    }

    public String pathToString() {
        if (path.isEmpty()) {
            return "";
        }

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < path.size(); i++) {
            sb.append(path.get(i));
            if (i < path.size() - 1) {
                sb.append(" -> ");
            }
        }
        return sb.toString();
    }

    public String pathToCoordinates(int graphSize) {
        if (path.isEmpty()) {
            return "";
        }

        int cols = (int) Math.sqrt(graphSize);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < path.size(); i++) {
            int index = path.get(i).index;
            int row = index / cols;
            int col = index % cols;
            sb.append("(").append(row).append(",").append(col).append(")");
            if (i < path.size() - 1) {
                sb.append(" -> ");
            }
        }
        return sb.toString();
    }

    @Override
    public String toString() {
        return String.format("Caminho: %s, Custo: %.2f, Nós expandidos: %d, Tempo: %.2fms",
                pathToString(), totalCost, nodesExpanded, executionTime / 1_000_000.0);
    }
}

abstract class SearchAlgorithm {

    public abstract SearchResult findPath(int[][] graph, Node start, Node goal);

    public abstract String getAlgorithmName();

    public abstract String getHeuristicName();

    protected boolean isTraversable(int[][] graph, int from, int to) {
        return graph[from][to] > 0;
    }

    protected List<Node> getNeighbors(Node node, int[][] graph) {
        List<Node> neighbors = new ArrayList<>();
        int size = graph.length;

        for (int i = 0; i < size; i++) {
            if (isTraversable(graph, node.index, i)) {
                neighbors.add(new Node(i, graph[node.index][i]));
            }
        }
        return neighbors;
    }

    protected List<Node> reconstructPath(Node goal) {
        List<Node> path = new ArrayList<>();
        Node current = goal;

        while (current != null) {
            path.add(0, current);
            current = current.parent;
        }
        return path;
    }

    protected double calculatePathCost(List<Node> path) {
        if (path.size() < 2) {
            return 0;
        }

        double totalCost = 0;
        for (int i = 0; i < path.size() - 1; i++) {
            Node next = path.get(i + 1);
            totalCost += next.g;
        }
        return totalCost;
    }
}

class Heuristic {

    public static double manhattan(int size, Node a, Node b) {

        int cols = (int) Math.sqrt(size);
        int x1 = a.index / cols;
        int y1 = a.index % cols;
        int x2 = b.index / cols;
        int y2 = b.index % cols;

        return Math.abs(x1 - x2) + Math.abs(y1 - y2);
    }

    public static double euclidean(int size, Node a, Node b) {

        int cols = (int) Math.sqrt(size);
        int x1 = a.index / cols;
        int y1 = a.index % cols;
        int x2 = b.index / cols;
        int y2 = b.index % cols;

        return Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2));
    }

    public static double calculate(String heuristicType, int size, Node a, Node b) {
        return switch (heuristicType.toLowerCase()) {
            case "manhattan" ->
                manhattan(size, a, b);
            case "euclidean", "euclidiana" ->
                euclidean(size, a, b);
            default ->
                0;
        };
    }
}

class BFS extends SearchAlgorithm {

    @Override
    public SearchResult findPath(int[][] graph, Node start, Node goal) {
        SearchResult result = new SearchResult();
        long startTime = System.nanoTime();

        int size = graph.length;
        boolean[] visited = new boolean[size];

        Queue<Node> queue = new LinkedList<>();
        queue.add(start);
        visited[start.index] = true;

        while (!queue.isEmpty()) {
            Node current = queue.poll();
            result.nodesExpanded++;

            if (current.equals(goal)) {
                result.path = reconstructPath(current);
                result.totalCost = calculatePathCost(result.path);
                break;
            }

            for (Node neighbor : getNeighbors(current, graph)) {
                if (!visited[neighbor.index]) {
                    visited[neighbor.index] = true;
                    neighbor.parent = current;
                    queue.add(neighbor);
                }
            }
        }

        result.executionTime = System.nanoTime() - startTime;
        return result;
    }

    @Override
    public String getAlgorithmName() {
        return "BFS";
    }

    @Override
    public String getHeuristicName() {
        return "";
    }
}

class DFS extends SearchAlgorithm {

    @Override
    public SearchResult findPath(int[][] graph, Node start, Node goal) {
        SearchResult result = new SearchResult();
        long startTime = System.nanoTime();

        int size = graph.length;
        boolean[] visited = new boolean[size];

        Stack<Node> stack = new Stack<>();
        stack.push(start);
        visited[start.index] = true;

        while (!stack.isEmpty()) {
            Node current = stack.pop();
            result.nodesExpanded++;

            if (current.equals(goal)) {
                result.path = reconstructPath(current);
                result.totalCost = calculatePathCost(result.path);
                break;
            }

            for (Node neighbor : getNeighbors(current, graph)) {
                if (!visited[neighbor.index]) {
                    visited[neighbor.index] = true;
                    neighbor.parent = current;
                    stack.push(neighbor);
                }
            }
        }

        result.executionTime = System.nanoTime() - startTime;
        return result;
    }

    @Override
    public String getAlgorithmName() {
        return "DFS";
    }

    @Override
    public String getHeuristicName() {
        return "";
    }
}

class Dijkstra extends SearchAlgorithm {

    @Override
    public SearchResult findPath(int[][] graph, Node start, Node goal) {
        SearchResult result = new SearchResult();
        long startTime = System.nanoTime();

        int size = graph.length;
        double[] dist = new double[size];
        boolean[] visited = new boolean[size];

        Arrays.fill(dist, Double.MAX_VALUE);
        dist[start.index] = 0;

        PriorityQueue<Node> queue = new PriorityQueue<>((a, b) -> Double.compare(dist[a.index], dist[b.index]));
        queue.add(start);

        while (!queue.isEmpty()) {
            Node current = queue.poll();
            result.nodesExpanded++;

            if (visited[current.index]) {
                continue;
            }
            visited[current.index] = true;

            if (current.equals(goal)) {
                result.path = reconstructPath(current);
                result.totalCost = dist[goal.index];
                break;
            }

            for (Node neighbor : getNeighbors(current, graph)) {
                if (!visited[neighbor.index]) {
                    double newDist = dist[current.index] + graph[current.index][neighbor.index];
                    if (newDist < dist[neighbor.index]) {
                        dist[neighbor.index] = newDist;
                        neighbor.parent = current;
                        neighbor.g = newDist;
                        queue.add(neighbor);
                    }
                }
            }
        }

        result.executionTime = System.nanoTime() - startTime;
        return result;
    }

    @Override
    public String getAlgorithmName() {
        return "DIJKSTRA";
    }

    @Override
    public String getHeuristicName() {
        return "";
    }
}

class GreedyBestFirstSearch extends SearchAlgorithm {

    private final String heuristicType;

    public GreedyBestFirstSearch(String heuristicType) {
        this.heuristicType = heuristicType;
    }

    @Override
    public SearchResult findPath(int[][] graph, Node start, Node goal) {
        SearchResult result = new SearchResult();
        long startTime = System.nanoTime();

        int size = graph.length;
        boolean[] visited = new boolean[size];

        PriorityQueue<Node> queue = new PriorityQueue<>((a, b) -> Double.compare(a.h, b.h));

        start.h = Heuristic.calculate(heuristicType, size, start, goal);
        queue.add(start);
        visited[start.index] = true;

        while (!queue.isEmpty()) {
            Node current = queue.poll();
            result.nodesExpanded++;

            if (current.equals(goal)) {
                result.path = reconstructPath(current);
                result.totalCost = calculatePathCost(result.path);
                break;
            }

            for (Node neighbor : getNeighbors(current, graph)) {
                if (!visited[neighbor.index]) {
                    visited[neighbor.index] = true;
                    neighbor.parent = current;
                    neighbor.h = Heuristic.calculate(heuristicType, size, neighbor, goal);
                    queue.add(neighbor);
                }
            }
        }

        result.executionTime = System.nanoTime() - startTime;
        return result;
    }

    @Override
    public String getAlgorithmName() {
        return "GREEDY BEST-FIRST-SEARCH";
    }

    @Override
    public String getHeuristicName() {
        return heuristicType.toUpperCase();
    }
}

class AStarSearch extends SearchAlgorithm {

    private final String heuristicType;

    public AStarSearch(String heuristicType) {
        this.heuristicType = heuristicType;
    }

    @Override
    public SearchResult findPath(int[][] graph, Node start, Node goal) {
        SearchResult result = new SearchResult();
        long startTime = System.nanoTime();

        int size = graph.length;
        double[] gScore = new double[size];
        double[] fScore = new double[size];
        boolean[] visited = new boolean[size];

        Arrays.fill(gScore, Double.MAX_VALUE);
        Arrays.fill(fScore, Double.MAX_VALUE);

        gScore[start.index] = 0;
        fScore[start.index] = Heuristic.calculate(heuristicType, size, start, goal);

        PriorityQueue<Node> queue = new PriorityQueue<>((a, b) -> Double.compare(fScore[a.index], fScore[b.index]));
        queue.add(start);

        while (!queue.isEmpty()) {
            Node current = queue.poll();
            result.nodesExpanded++;

            if (current.equals(goal)) {
                result.path = reconstructPath(current);
                result.totalCost = gScore[goal.index];
                break;
            }

            if (visited[current.index]) {
                continue;
            }
            visited[current.index] = true;

            for (Node neighbor : getNeighbors(current, graph)) {
                if (!visited[neighbor.index]) {
                    double tentativeGScore = gScore[current.index] + graph[current.index][neighbor.index];

                    if (tentativeGScore < gScore[neighbor.index]) {
                        neighbor.parent = current;
                        gScore[neighbor.index] = tentativeGScore;
                        fScore[neighbor.index] = tentativeGScore
                                + Heuristic.calculate(heuristicType, size, neighbor, goal);
                        queue.add(neighbor);
                    }
                }
            }
        }

        result.executionTime = System.nanoTime() - startTime;
        return result;
    }

    @Override
    public String getAlgorithmName() {
        return "A*";
    }

    @Override
    public String getHeuristicName() {
        return heuristicType.toUpperCase();
    }
}

public class Main {

    public static void main(String[] args) {
        if (args.length != 3) {
            System.out.println("Uso: java Main <arquivo> \"origem\" \"destino\"");
            System.out.println("Exemplo: java Main teste_4x4.txt \"0,0\" \"2,2\"");
            return;
        }

        String inputFile = args[0];
        String originStr = args[1];
        String destinationStr = args[2];

        try {
            int[][] grafo = searchAndReadFile(inputFile);
            System.out.println("Arquivo lido com sucesso!");
            System.out.println("Tamanho do grafo: " + grafo.length + "x" + grafo[0].length);

            int originIndex = parseCoordinatesToIndex(originStr, grafo.length);
            int destinationIndex = parseCoordinatesToIndex(destinationStr, grafo.length);

            Node origin = new Node(originIndex);
            Node destination = new Node(destinationIndex);

            if (!isTraversable(grafo, origin.index, destination.index) && grafo[origin.index][destination.index] <= 0) {
                System.out.println(
                        "AVISO: Não há conexão direta entre origem e destino, mas algorithms podem encontrar caminho indireto.");
            }

            System.out.println("Origem: nó " + origin.index + " (coordenadas: " + originStr + ")");
            System.out.println("Destino: nó " + destination.index + " (coordenadas: " + destinationStr + ")");

            String outputFolderPath = createOutputFolder(inputFile);

            List<SearchAlgorithm> algorithms = new ArrayList<>();
            algorithms.add(new BFS());
            algorithms.add(new DFS());
            algorithms.add(new Dijkstra());
            algorithms.add(new GreedyBestFirstSearch("manhattan"));
            algorithms.add(new GreedyBestFirstSearch("euclidean"));
            algorithms.add(new AStarSearch("manhattan"));
            algorithms.add(new AStarSearch("euclidean"));

            for (SearchAlgorithm algorithm : algorithms) {
                System.out.println("\nExecutando: " + algorithm.getAlgorithmName()
                        + (algorithm.getHeuristicName().isEmpty() ? "" : " (" + algorithm.getHeuristicName() + ")"));

                SearchResult result = algorithm.findPath(grafo, origin, destination);
                System.out.println("Resultado: " + result);
                if (!result.path.isEmpty()) {
                    System.out.println("Caminho em coordenadas: " + result.pathToCoordinates(grafo.length));
                }

                saveResult(inputFile, algorithm, result, originStr, destinationStr, outputFolderPath, grafo.length);
            }

            System.out.println("\nTodos os algorithms executados e resultados salvos em: " + outputFolderPath);

        } catch (IOException e) {
            System.err.println("Erro ao ler o arquivo: " + e.getMessage());
        } catch (Exception e) {
            System.err.println("Erro inesperado: " + e.getMessage());
        }
    }

    private static String createOutputFolder(String inputFile) {
        File file = new File(inputFile);
        String filename = file.getName();
        String basename = filename.substring(0, filename.lastIndexOf('.'));

        String path = "./output/" + basename + "/";
        File pastaOutput = new File(path);

        if (!pastaOutput.exists()) {
            if (pastaOutput.mkdirs()) {
                System.out.println("Pasta output criada: " + pastaOutput.getAbsolutePath());
            } else {
                System.err.println("AVISO: Não foi possível criar a pasta " + path + ". Salvando no diretório atual.");
                path = "./";
            }
        }

        return path;
    }

    private static void saveResult(String inputFile, SearchAlgorithm algorithm,
            SearchResult result, String origin, String destination, String outputFolder, int graphSize) {

        String filename = new File(inputFile).getName();
        String basename = filename.substring(0, filename.lastIndexOf('.'));
        String extension = getExtensionAlgorithm(algorithm);

        String outputFile = outputFolder + basename + "." + extension;

        try (PrintWriter writer = new PrintWriter(new FileWriter(outputFile))) {
            writer.println("ALGORITMO: " + algorithm.getAlgorithmName());

            String heuristic = algorithm.getHeuristicName();
            writer.println("HEURISTICA: " + (heuristic.isEmpty() ? "" : heuristic));

            writer.println("ORIGEM: " + origin);
            writer.println("DESTINO: " + destination);

            String pathStr = result.path.isEmpty() ? "" : result.pathToString();
            writer.println("CAMINHO (NOS): " + pathStr);

            String pathCoords = result.path.isEmpty() ? "" : result.pathToCoordinates(graphSize);
            writer.println("CAMINHO (COORDENADAS): " + pathCoords);

            writer.println("CUSTO: " + (result.path.isEmpty() ? "" : String.format("%.2f", result.totalCost)));
            writer.println("NOS EXPANDIDOS: " + result.nodesExpanded);
            writer.println("TEMPO (ms): " + String.format("%.2f", result.executionTime / 1_000_000.0));

            System.out.println("Resultado salvo em: " + outputFile);

        } catch (IOException e) {
            System.err.println("Erro ao salvar resultado em " + outputFile + ": " + e.getMessage());
        }
    }

    private static String getExtensionAlgorithm(SearchAlgorithm algorithm) {
        String name = algorithm.getAlgorithmName();
        String heuristic = algorithm.getHeuristicName().toLowerCase();

        return switch (name) {
            case "BFS" ->
                "bfs";
            case "DFS" ->
                "dfs";
            case "DIJKSTRA" ->
                "dijkstra";
            case "GREEDY BEST-FIRST-SEARCH" ->
                "gbs." + heuristic;
            case "A*" ->
                "a." + heuristic;
            default ->
                name.toLowerCase();
        };
    }

    private static boolean isTraversable(int[][] graph, int from, int to) {
        return graph[from][to] > 0;
    }

    private static int parseCoordinatesToIndex(String coordinates, int graphSize) {
        String split = coordinates.replaceAll("[\"()\\s]", "");
        String[] parts = split.split(",");

        if (parts.length != 2) {
            throw new IllegalArgumentException("Formato de coordenadas inválido. Use: x,y");
        }

        try {
            int x = Integer.parseInt(parts[0]);
            int y = Integer.parseInt(parts[1]);

            int cols = (int) Math.sqrt(graphSize);
            if (cols * cols != graphSize) {
                throw new IllegalArgumentException("Grid não é quadrado perfeito: " + graphSize + " nós");
            }

            int index = x * cols + y;

            if (index < 0 || index >= graphSize) {
                throw new IllegalArgumentException("Coordenadas fora do grid: (" + x + "," + y + ")");
            }

            return index;

        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Coordenadas devem ser números inteiros");
        }
    }

    public static int[][] searchAndReadFile(String filename) throws IOException {
        File file = new File(filename);
        if (!file.exists()) {
            throw new FileNotFoundException("Arquivo não encontrado: " + filename);
        }

        List<String> lines = Files.readAllLines(Paths.get(filename));
        lines.removeIf(linha -> linha.trim().isEmpty());

        if (lines.isEmpty()) {
            throw new IOException("Arquivo vazio: " + filename);
        }

        int numNodes = lines.size();
        int[][] graph = new int[numNodes][numNodes];

        for (int i = 0; i < numNodes; i++) {
            String linha = lines.get(i).trim();
            String[] values = linha.split("\\s+");

            if (values.length != numNodes) {
                throw new IOException("Linha " + (i + 1) + " tem " + values.length
                        + " colunas, mas esperava " + numNodes + ". A matriz deve ser quadrada.");
            }

            for (int j = 0; j < numNodes; j++) {
                try {
                    graph[i][j] = Integer.parseInt(values[j]);
                } catch (NumberFormatException e) {
                    throw new IOException("Valor inválido na posição [" + i + "," + j + "]: '" + values[j] + "'");
                }
            }
        }

        return graph;
    }
}
