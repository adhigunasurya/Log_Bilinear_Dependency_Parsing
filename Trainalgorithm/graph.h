#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include "edge.cpp"
#include <vector>
#include <utility>
#include <cmath>
#include <map>
#include <Eigen/Dense> //use Eigen for the adjacency matrix

bool detect_cycle(std::vector<Edge>& vector_of_edges);
bool does_reverse_edge_exist(Edge& edge, std::vector<Edge>& vector_of_edges);
std::vector<int> get_cycle(std::vector<Edge>& vector_of_edges); //precondition: this function must only be called only when there is a cycle on the minimum incoming edges
bool contains_edge(int source, int destination, std::vector<Edge>& vector_of_edges);
double get_weight(int source, int destination, std::vector<Edge>& vector_of_edges);
bool does_exist_edge_ending_in(int destination, std::vector<Edge>& vector_of_edges);
int getSource(int destination, std::vector<Edge>& vector_of_edges);
class Graph {
private:
	int no_of_vertices;
	int root;
	Eigen::ArrayXXd adjacency_matrix;
	void remove_edges_ending_in_root();
public:
	Graph(int no_of_vertices, int root); 
	Graph(int no_of_vertices, int root, std::vector<Edge>& Edges); 
	Graph(const Graph& ref); //cctor
	std::vector<Edge> get_min_incoming_edges(); //get the minimum incoming edges for each vertex
	Graph& operator=(const Graph& ref); //copy assignment
	~Graph();
	//important to remember: during Connect, MAKE SURE THAT IT'S NOT CONNECTING A VERTEX TO ITSELF, AND MAKE SURE THAT THE DESTINATION IS NOT THE ROOT
	void Connect(int source, int destination, double weight); //make a connection, IF THE EDGE ALREADY EXISTS, THE WEIGHT WILL BE UPDATED IF PARAMETER WEIGHT IS LOWER THAN THE CURRENT WEIGHT
	void Connect(int source, int destination);
	void Disconnect(int source, int destination); //might be used to disconnect all edges that end at the root
	int get_no_of_vertices();
	Eigen::ArrayXXd get_adjacency_matrix();
	void printGraph(); //print graph to std::cout
	void printGraph(std::ofstream& file); //print graph to an external file
	double get_minimum_incoming_weight(int node);
	std::vector<Edge> get_MST(); //use Chu-Liu/Edmonds algorithm to find the minimum spanning tree of a directed graph
};

#endif