#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include "edge.cpp"
#include <vector>
#include <utility>
#include <cmath>
#include <map>

//Initially there will be no parallel edges in this implementation since connecting two edges that already exist will result in the smaller weight edge being replaced
int n_of_edges(int n); //function to calculate the supposed number of edges given the number of word tokens
template <class T> 
bool does_exist(T element, std::vector<T>& arr);
bool does_edge_exist(int sourceS, int destinationS, std::vector<Edge>& arr);
Edge * get_Edge_Vector(int source, int destination, std::vector<Edge>& arr);
void delete_Edge(int source, int destination, std::vector<Edge>& arr);
//this is the class that defines a directed graph that will be used to get the minimum spanning tree of the graph and get the dependency tree
class Graph {
private:
	int no_of_vertices; // the number of vertices in the graph
	std::vector<Edge> arr_edge; //initialise an empty array (well, technically a vector) of edges
public:
	Graph(int n); //initialise the graph based on the number of tokens. Note that n is the number of token, so the actual number of edges will be n*(n+1) and the number of vertices is n+1
	Graph(int v, std::vector<Edge>& vect);	//create a graph given the number of vertices and all the edges
	Graph(const Graph& ref); //cctor
	Graph& operator=(const Graph& ref); //copy assignment
	~Graph();
	void Disconnect(int source, int destination); //remove an edge connecting source and destination
	void Connect(int source, int destination, double weight); //make a connection. //IF THE EDGE ALREADY EXISTS, THE WEIGHT WILL BE UPDATED IF PARAMETER WEIGHT IS LOWER THAN THE CURRENT WEIGHT
	void Connect(int source, int destination); //make a connection with default weight. Do nothing if an edge already connects the two vertices or if the source and destination indices are out of bound
	bool isEdge(int source, int destination); //check if there is an edge connecting two vertices
	Edge * get_Edge(int source, int destination); //get the edge connecting the source and the destination, return null_ptr if none exists
	int edgeCount(); //returns the number of edges in the graph
	int vertexCount(); //returns the number of vertices in the graph
	void setEdgeWeight(int source, int destination, double weight);
	void setEdgeWeightDefault(int source, int destination);
	void printGraph();
	Edge * get_min_incoming_edge(int node); //get the minimum incoming edge that ends at a particular node
	std::vector<Edge> get_min_incoming_edge(); //get all the edges in the graph that form the minimum incoming edge
	void remove_edges_ending_in(int node); //remove all edges that end at the particular node
	bool is_edge_ending_in(int destination);
	bool does_cycle_exist();
	bool does_cycle_exist(int root); //detect if there is a cycle, with the first node inspected is the root
	std::pair<int, int> get_First_Cycle(); //get a pair of nodes that forms the first detected cycle in the graph. IMPORTANT: Only call whenever you are sure that the graph does have a cycle (using the previous function)
	std::vector<Edge> get_Edges_Starting_From(int node); //get all the edges that start from a particular node
	Edge * get_Edge_By_Weight(double weight, Edge& alternative_1, Edge& alternative_2); // pick which edge (between alternative 1 and 2) has the same weight as the specified weight
	std::vector<Edge> getMST(int root); //calculate the minimum spanning tree for this directed graph (the negative log of the probability) using the Chu-Liu/Edmonds Algorithm
	std::vector<std::string> getPairs(int root); //get the parsing result as pairs
	//NEED A NEW METHOD TO REMOVE PARALLEL EDGES AND TAKE THE SMALLEST ONE
};

#endif //GRAPH_H