#include "graph.h"
#include <exception> //exception is important since we are allocating memories dynamically where an exception may be thrown
#include <stack>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

Graph::Graph(int no_of_vertices, int root) {
	this->no_of_vertices = no_of_vertices;
	this->root = root;
	adjacency_matrix = Eigen::ArrayXXd::Zero(no_of_vertices, no_of_vertices);
	//initialise the adjacency matrix with default weight
	for (int i = 0; i < adjacency_matrix.rows(); i++) {
		for (int j = 0; j < adjacency_matrix.cols(); j++) {
			adjacency_matrix(i, j) = DEFAULT_WEIGHT;
		}
	}
}

Graph::Graph(int no_of_vertices, int root, std::vector<Edge>& Edges)  {
	this->no_of_vertices = no_of_vertices;
	this->root = root;
	adjacency_matrix = Eigen::ArrayXXd::Zero(no_of_vertices, no_of_vertices);
	for (int i = 0; i < adjacency_matrix.rows(); i++) {
		for (int j = 0; j < adjacency_matrix.cols(); j++) {
			adjacency_matrix(i, j) = DEFAULT_WEIGHT;
		}
	}
	for (std::vector<Edge>::iterator it = Edges.begin(); it != Edges.end(); it++) {
		if (it->getWeight() < adjacency_matrix(it->getSource(), it->getDestination())) {
			adjacency_matrix(it->getSource(), it->getDestination()) = it->getWeight();
		}
	}
}

Graph::Graph(const Graph& ref) {
	this->no_of_vertices = ref.no_of_vertices;
	this->root = ref.root;
	adjacency_matrix = Eigen::ArrayXXd::Zero(no_of_vertices, no_of_vertices);
	for (int i = 0; i < adjacency_matrix.rows(); i++) {
		for (int j = 0; j < adjacency_matrix.cols(); j++) {
			adjacency_matrix(i, j) = ref.adjacency_matrix(i, j);
		}
	}
}

Graph& Graph::operator=(const Graph& ref) {
	this->no_of_vertices = ref.no_of_vertices;
	this->root = ref.root;
	//resize the current adjacency matrix to match the size of the reference
	adjacency_matrix.resize(no_of_vertices, no_of_vertices);
	for (int i = 0; i < adjacency_matrix.rows(); i++) {
		for (int j = 0; j < adjacency_matrix.cols(); j++) {
			adjacency_matrix(i, j) = ref.adjacency_matrix(i, j);
		}
	}
}

Graph::~Graph() {
	//resize the adjacency matrix to become empty
	adjacency_matrix.resize(0, 0);
}

//important to remember: during Connect, MAKE SURE THAT IT'S NOT CONNECTING A VERTEX TO ITSELF, AND MAKE SURE THAT THE DESTINATION IS NOT THE ROOT
//make a connection, IF THE EDGE ALREADY EXISTS, THE WEIGHT WILL BE UPDATED IF PARAMETER WEIGHT IS LOWER THAN THE CURRENT WEIGHT
void Graph::Connect(int source, int destination, double weight) {
	if (source == destination) {
		//do nothing
	} //else if (destination == root) {
		else if (destination == root && weight != DEFAULT_WEIGHT) {
		//do nothing
	} else if (source >= no_of_vertices || destination >= no_of_vertices || source < 0 || destination < 0) {
		//do nothing
	} else {
		if (weight < adjacency_matrix(source, destination)) {
			adjacency_matrix(source, destination) = weight;
		}
	}
}

void Graph::Connect(int source, int destination) {
	if (source == destination) {
		//do nothing
	} else if (source >= no_of_vertices || destination >= no_of_vertices || source < 0 || destination < 0) {
		//do nothing
	} else {
		if (DEFAULT_WEIGHT < adjacency_matrix(source, destination)) {
			adjacency_matrix(source, destination) = DEFAULT_WEIGHT;
		}
	}
}

void Graph::Disconnect(int source, int destination) {
	if (source >= no_of_vertices || destination >= no_of_vertices || source < 0 || destination < 0) {
		//do nothing
	} else {
		adjacency_matrix(source, destination) = DEFAULT_WEIGHT;
	}
}

void Graph::remove_edges_ending_in_root() {
	for (int i = 0; i < no_of_vertices; i++) {
		Disconnect(i, root);
	}
}

int Graph::get_no_of_vertices() {
	return no_of_vertices;
}

Eigen::ArrayXXd Graph::get_adjacency_matrix() {
	return adjacency_matrix;
}

//get the minimum incoming edges for each vertex EXCEPT FOR THE ROOT
std::vector<Edge> Graph::get_min_incoming_edges() {
	std::vector<Edge> result;
	for (int i = 0; i < adjacency_matrix.cols(); i++) {
		double min = DEFAULT_WEIGHT + 1.0;
		int min_idx = -1;
		for (int j = 0; j < adjacency_matrix.rows(); j++) {
			if (i != j) {
				if (adjacency_matrix(j, i) < min) {
					min = adjacency_matrix(j, i);
					min_idx = j;
				}
			}	
		}
		if (min != DEFAULT_WEIGHT) {
			result.push_back(*(new Edge(min_idx, i, min)));
		}
	}
	return result;
}

void Graph::printGraph() {
	std::cout << adjacency_matrix << std::endl;
}

void Graph::printGraph(std::ofstream& file) {
	file << adjacency_matrix << std::endl;
}

bool does_reverse_edge_exist(Edge& edge, std::vector<Edge>& vector_of_edges) {
	for (int i = 0; i < vector_of_edges.size(); i++) {
		if (vector_of_edges.at(i).getSource() == edge.getDestination() && vector_of_edges.at(i).getDestination() == edge.getSource()) {
			return true;
		}
	}
	return false;
}

bool detect_cycle(std::vector<Edge>& vector_of_edges) {
	for (int i = 0; i < vector_of_edges.size(); i++) {
		if (does_reverse_edge_exist(vector_of_edges[i], vector_of_edges)) {
			return true;
		}
	}
	return false;
}

bool does_exist_edge_ending_in(int destination, std::vector<Edge>& vector_of_edges) {
	for (std::vector<Edge>::iterator it = vector_of_edges.begin(); it != vector_of_edges.end(); it++) {
		if (it->getDestination() == destination) {
			return true;
		}
	}
	return false;
}

//precondition: this function must only be called only when there is a cycle on the graph
//if the result's elements are both -1, then something is wrong
std::vector<int> get_cycle(std::vector<Edge>& vector_of_edges) {
	std::vector<int> result(2, -1);
	for (std::vector<Edge>::iterator it = vector_of_edges.begin(); it != vector_of_edges.end(); it++) {
		if (does_reverse_edge_exist(*it, vector_of_edges)) {
			result[0] = it->getSource();
			result[1] = it->getDestination();
			break;
		}
	}
	return result;
}

bool contains_edge(int source, int destination, std::vector<Edge>& vector_of_edges) {
	for (std::vector<Edge>::iterator it = vector_of_edges.begin(); it != vector_of_edges.end(); it++) {
		if (it->getSource() == source && it->getDestination() == destination) {
			return true;
		}
	}
	return false;
}

double get_weight(int source, int destination, std::vector<Edge>& vector_of_edges) {
	for (std::vector<Edge>::iterator it = vector_of_edges.begin(); it != vector_of_edges.end(); it++) {
		if (it->getSource() == source && it->getDestination() == destination) {
			return it->getWeight();
		}
	}
	return DEFAULT_WEIGHT;
}

int getSource(int destination, std::vector<Edge>& vector_of_edges) {
	for (std::vector<Edge>::iterator it = vector_of_edges.begin(); it != vector_of_edges.end(); it++) {
		if (it->getDestination() == destination) {
			return it->getSource();
		}
	}
}

double Graph::get_minimum_incoming_weight(int node) {
	double min = DEFAULT_WEIGHT+1;
	for (int i = 0; i < adjacency_matrix.rows(); i++) {
		if (i != node) {
			if (adjacency_matrix(i, node) < min) {
				min = adjacency_matrix(i, node);
			}
		}
	}
	return min;
}

std::vector<Edge> Graph::get_MST() {
	this->remove_edges_ending_in_root();
	std::vector<Edge> result;
	std::vector<Edge> curr_min_incoming_edges = this->get_min_incoming_edges();
	if (!detect_cycle(curr_min_incoming_edges)) {
		result = curr_min_incoming_edges;
	} else {
		std::vector<int> cycle = get_cycle(curr_min_incoming_edges);
		std::map<int, int> mapping;
		int j = 0;
		int root_idx = -1;
		//map the indices of the current graph to new indices where the vertices that form a cycle are condensed to become one vertex
		for (int i = 0; i < no_of_vertices; i++) {
			if (i != cycle[0] && i != cycle[1]) {
				if (i == root) {
					root_idx = j;
				}
				mapping[i] = j;
				j++;
			} else { //map the vertices that form a cycle to the last vertex of the condensed graph
				mapping[i] = no_of_vertices - 2;
			}
		}
		if (root_idx == -1) {
			std::cout << "Something wrong as the root is not detected" << std::endl;
		} else {
			Graph recursive_graph(no_of_vertices-1, root_idx);
			for (int i = 0; i < adjacency_matrix.rows(); i++) {
				for (int j = 0; j < adjacency_matrix.cols(); j++) {
					if (i != j) {
						if (i != cycle[0] && i != cycle[1] && j != cycle[0] && j != cycle[1]) {
							recursive_graph.Connect(mapping.at(i), mapping.at(j), adjacency_matrix(i, j));
						} else if ((i == cycle[0] || i == cycle[1]) && j != cycle[0] && j != cycle[1]) {
							recursive_graph.Connect(mapping.at(i), mapping.at(j), adjacency_matrix(i, j));
						} else if (i != cycle[0] && i != cycle[1] && (j == cycle[0] || j == cycle[1])) {
							recursive_graph.Connect(mapping.at(i), mapping.at(j), adjacency_matrix(i, j) - get_minimum_incoming_weight(j));
						} else if ((i == cycle[0] && j == cycle[1]) || (i == cycle[1] && j == cycle[0])){
							//do nothing, the vertices are not included as they are condensed
						} else {
							std::cout << "something is likely wrong" << std::endl;
						}
					}
				}
			}
			std::vector<Edge> recursive_min_incoming_edges = recursive_graph.get_MST();
			//include both vertices connected by the cycle
			result.push_back(*(new Edge(cycle[0], cycle[1], adjacency_matrix(cycle[0], cycle[1]))));
			result.push_back(*(new Edge(cycle[1], cycle[0], adjacency_matrix(cycle[1], cycle[0]))));
			bool done = false;
			for (int i = 0; i < adjacency_matrix.rows(); i++) {
				for (int j = 0; j < adjacency_matrix.cols(); j++) {
					if (i != j) {
						//do nothing if both vertices are in the cycle
						if ((i == cycle[0] && j == cycle[1]) || (i == cycle[1] && j == cycle[0])) {
							//do nothing
							//result.push_back(*(new Edge(i, j, adjacency_matrix(i, j))));
						} //case 2: both edges are not in the cycle
						else if (contains_edge(mapping.at(i), mapping.at(j), recursive_min_incoming_edges) && i != cycle[0] && i != cycle[1] && j != cycle[0] && j != cycle[1]) {
							result.push_back(*(new Edge(i, j, adjacency_matrix(i, j))));
						} //case 3: the source is in the cycle but the destination is not
						else if (contains_edge(mapping.at(i), mapping.at(j), recursive_min_incoming_edges) && (i == cycle[0] || i == cycle[1]) && j != cycle[0] && j != cycle[1]) {
							if (!does_exist_edge_ending_in(j, result) && adjacency_matrix(i, j) == get_weight(mapping.at(i), mapping.at(j), recursive_min_incoming_edges)) {
								result.push_back(*(new Edge(i, j, adjacency_matrix(i, j))));
							}
						} //case 4: the source is not in the cycle but the destination is 
						else if (contains_edge(mapping.at(i), mapping.at(j), recursive_min_incoming_edges) && (j == cycle[0] || j == cycle[1]) && i != cycle[0] && i != cycle[1]) {
							if (!done && (adjacency_matrix(i, j) - get_minimum_incoming_weight(j)) == get_weight(mapping.at(i), mapping.at(j), recursive_min_incoming_edges)) {
								done = true;
								for (int k = 0; k < result.size(); k++) {
									if (result.at(k).getDestination() == j) {
										result.erase(result.begin() + k);
										break;
									}
								}
								result.push_back(*(new Edge(i, j, adjacency_matrix(i, j))));
							}
						} 
					}
				}
			}
		} 
	}
	return result;
}


