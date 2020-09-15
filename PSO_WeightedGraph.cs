using System;
using System.Collections.Generic;
using System.Linq;
using FibonacciHeap;
using MathNet.Numerics.LinearAlgebra;

namespace Theory_Definitions
{
    /// <summary>
    /// A simple, directed weighted graph. Weights are assumed positive.
    /// </summary>
    public class PSO_WeightedGraph
    {
        #region Properties
        public int Size { get; }
        public Matrix<double> AdjacencyMatrix { get; }
        #endregion

        #region Constructors
        public PSO_WeightedGraph(int nodes)
        {
            this.Size = nodes;
        }

        public PSO_WeightedGraph(int nodes, ICollection<Tuple<int, int, double>> edges) : this(nodes)
        {
            AdjacencyMatrix = Matrix<double>.Build.Sparse(nodes, nodes);
            foreach (var edge in edges)
            {
                AdjacencyMatrix[edge.Item1, edge.Item2] = edge.Item3;
            }
        }

        private PSO_WeightedGraph(Matrix<double> m) : this(m.ColumnCount)
        {
            AdjacencyMatrix = m;
        }
        #endregion

        public IEnumerable<int> GetNeighbors(int vertex)
        {
            return
                from i in Enumerable.Range(0, Size)
                where AdjacencyMatrix[vertex, i] > 0
                select i;
        }

        /// <summary>
        /// Implementation of uniform-cost search algorithm using a Fibonacci heap data structure to minimize computing times
        /// </summary>
        /// <param name="source">The node on which to start the search</param>
        /// <param name="destination">The node we try to find the shortest path to, starting from source</param>
        public List<int> LeastCostPath(int source, int destination)
        {
            var Predecessor = new Dictionary<int, int>();
            var Distance = new Dictionary<int, double>();
            var Frontier = new FibonacciHeap<int, double>(0);
            Frontier.Insert(new FibonacciHeapNode<int, double>(source, 0));
            var Explored = new List<int>();
            Predecessor.Add(source, -1); //value of -1 indicates this node has no predecessors

            while (true)
            {
                if (Frontier.IsEmpty())
                    throw new Exception("LeastCostPath: Failed to find path between source (" + source + ") and destination (" + destination + ").");

                var minNode = Frontier.RemoveMin();
                if (minNode.Data == destination)
                {
                    List<int> LCP = new List<int> { minNode.Data };
                    int pred = Predecessor[minNode.Data];
                    while (pred != -1)
                    {
                        LCP.Add(pred);
                        pred = Predecessor[pred];
                    }
                    LCP.Reverse();
                    return LCP;
                }

                Explored.Add(minNode.Data);
                foreach (int neighbor in this.GetNeighbors(minNode.Data))
                    if (!Explored.Contains(neighbor))
                    {
                        var neighborCost = minNode.Key + AdjacencyMatrix[minNode.Data, neighbor];
                        Frontier.Insert(new FibonacciHeapNode<int, double>(neighbor, neighborCost));
                        if (Distance.TryGetValue(neighbor, out double cost))
                        {
                            if (neighborCost < cost)
                            {
                                Predecessor[neighbor] = minNode.Data;
                                Distance[neighbor] = neighborCost;
                            }
                        }
                        else
                        {
                            Predecessor.Add(neighbor, minNode.Data);
                            Distance.Add(neighbor, neighborCost);
                        }
                    }
            }
        }

        /// <summary>
        /// Implementation of Yen's algorithm which finds the shortest path between nodes, and then the
        /// K-1 shortest deviations from this path. Paths returned will be simple and loopless
        /// </summary>
        /// <param name="K">The number of shortest paths to find.</param>
        /// <param name="source">The node on which to start the search</param>
        /// <param name="destination">The node we try to find the shortest path to, starting from source</param>
        /// <returns></returns>
        public List<List<int>> KShortestPaths(int K, int source, int destination)
        {
            List<List<int>> ShortestPaths = new List<List<int>>();
            var PotentialPaths = new FibonacciHeap<List<int>, double>(0);
            ShortestPaths.Add(LeastCostPath(source, destination));

            //now find next K-1 shortest paths
            foreach (int k in Enumerable.Range(1, K - 1))
            {
                //The spur node ranges from the first node to the next to last node in the previous k-shortest path.
                int spurNodeCount = ShortestPaths[k - 1].Count - 1;
                foreach (int i in Enumerable.Range(0, spurNodeCount))
                {
                    int spurNode = ShortestPaths[k - 1][i];
                    List<int> rootPath = ShortestPaths[k - 1].GetRange(0, i + 1);
                    PSO_WeightedGraph AlteredGraph = this.Clone();

                    //temporarily remove edges to avoid retracing our steps
                    foreach (List<int> shortPath in ShortestPaths)
                        if (rootPath.SequenceEqual(shortPath.Take(i + 1)))
                            AlteredGraph.AdjacencyMatrix[shortPath[i], shortPath[i + 1]] = 0;

                    //To avoid looping back over a previous path, we disconnect nodes in the root path (except the spur node)
                    //by setting the weights of the edges that connect them to the graph to 0
                    foreach (int x in Enumerable.Range(0, Size).Where(a => a != spurNode & rootPath.Contains(a)))
                    {
                        var v = Vector<double>.Build.Sparse(Size);
                        AlteredGraph.AdjacencyMatrix.SetColumn(x, v);
                        AlteredGraph.AdjacencyMatrix.SetRow(x, v);
                    }

                    //build spur path and connect the spur path to the root
                    List<int> spurPath = new List<int>();
                    //finding the least cost path may fail due to removing the edges above; just ignore and continue
                    try { spurPath = AlteredGraph.LeastCostPath(spurNode, destination); }
                    catch (Exception ex) { break; }

                    List<int> totalPath = rootPath;
                    totalPath.AddRange(spurPath.Where(node => node != spurNode).ToList());
                    PotentialPaths.Insert(new FibonacciHeapNode<List<int>, double>(totalPath, this.PathCost(totalPath)));
                }

                if (PotentialPaths.IsEmpty())
                    break;

                ShortestPaths.Add(PotentialPaths.RemoveMin().Data);
            }

            return ShortestPaths;
        }

        /// <summary>
        /// Finds the cost of traversing the path along the nodes specified by the input parameter
        /// </summary>
        /// <param name="nodes">The nodes of the path</param>
        /// <returns></returns>
        public double PathCost(List<int> nodes)
        {
            double cost = 0;
            foreach (int i in Enumerable.Range(0, nodes.Count - 1))
            {
                var addCost = AdjacencyMatrix[nodes[i], nodes[i + 1]];
                if (addCost > 0)
                    cost += addCost;
                else
                    throw new ArgumentException("Cannot calculate path cost; nodes given do not yield a valid path.");
            }
            return cost;
        }

        #region Miscellaneous
        public PSO_WeightedGraph Clone()
        {
            return new PSO_WeightedGraph(AdjacencyMatrix.Clone());
        }
        #endregion
    }
}
