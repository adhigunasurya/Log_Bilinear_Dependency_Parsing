#include "graph.h"
#include <exception> //exception is important since we are allocating memories dynamically where an exception may be thrown
#include <stack>

//function to calculate resulting number of edges
template <class T> 
bool does_exist(T element, std::vector<T> arr) {
	for (typename std::vector<T>::iterator it = arr.begin() ; it != arr.end(); it++) {
		if (*it == element) {
			return true;
		}
	}
	return false;
}

int n_of_edges(int n) { //return the number of edges that will be in the directed graph
	return ((n+1) * n);
}

int Graph::edgeCount() {
	return no_of_edges;
}

int Graph::vertexCount() {
	return no_of_vertices;
}

//sole constructor
Graph::Graph(int n) { //IMPORTANT: n is the number of tokens, thus the size of the graph will be n+1
	no_of_edges = n_of_edges(n); 
	no_of_vertices = n + 1; //n+1 since we add the dummy "root" token
	try {
		//create edges from every vertex to every other vertex, all edges are assigned default weight (remember we are making a strongly connected graph)
		for (int i = 0; i < no_of_vertices; i++) {
			for (int j = 0; j < no_of_vertices; j++) {
				if (i != j) { //avoid creating edges from a vertex to itself
					arr_edge.push_back(*(new Edge(i, j)));
				}
			}
		}
	} catch (std::exception& e) {
		std::cout << "Error detected during vector of edges initialisation. Exception description: " << e.what() << std::endl;
	
	}
}

Graph::Graph(const Graph& ref) { //cctor
	//set the size
	this -> no_of_edges = ref.no_of_edges;
	this -> no_of_vertices = ref.no_of_vertices;
	(this -> arr_edge).clear();
	arr_edge = *(new std::vector<Edge> (ref.arr_edge)); //copy the vector of edges that belong to ref
}

Graph& Graph::operator=(const Graph& ref) { //copy assignment
	this -> no_of_edges = ref.no_of_edges;
	this -> no_of_vertices = ref.no_of_vertices;
	//delete all the elements from the object's array of edges
	(this -> arr_edge).clear();
	//initialise the array of edges using copy constructor with ref.arr_edge as its sole parameter
	this -> arr_edge = *(new std::vector<Edge> (ref.arr_edge));
}

Graph::~Graph() { //destructor
	(this -> arr_edge).clear();
}

void Graph::Disconnect(int sourceS, int destinationS) { //inefficient operation as the whole vector would then have to be repositioned if we delete anything other than the last element
	if (sourceS >= no_of_vertices || destinationS >= no_of_vertices) {
		//do nothing
		std::cout << "source or destination indices are out of bound" << std::endl;
	} else {
		int i = 0;
		for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
			if ((it -> getSource()) == sourceS && (it -> getDestination()) == destinationS) {
				arr_edge.erase(arr_edge.begin() + i); //erase the element
				no_of_edges--;
				break;
			} else {
				i++;
			}
		}
	}	
}

bool Graph::isEdge(int sourceS, int destinationS) {
	if (sourceS >= no_of_vertices || destinationS >= no_of_vertices) {
		//do nothing
		std::cout << "source or destination indices are out of bound" << std::endl;
		return false; //there is no edge connecting them as the indices are invalid
	}
	for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
		if ((it -> getSource()) == sourceS && (it -> getDestination()) == destinationS) {
			return true; // there is an edge
		} 
	}
	return false; // there is no edge connecting them
}

void Graph::setEdgeWeight(int sourceS, int destinationS, double weightS) {
	if (sourceS >= no_of_vertices || destinationS >= no_of_vertices) {
		//do nothing
		std::cout << "source or destination indices are out of bound" << std::endl;
	} else {
		for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
			if ((it -> getSource()) == sourceS && (it -> getDestination()) == destinationS) {
				it -> setWeight(weightS); //if found, update the weight
				break;
			} 
		}
	}
	//otherwise, do nothing
}

void Graph::setEdgeWeightDefault(int sourceS, int destinationS) {
	if (sourceS >= no_of_vertices || destinationS >= no_of_vertices) {
		//do nothing
		std::cout << "source or destination indices are out of bound" << std::endl;
	} else {
		for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
			if ((it -> getSource()) == sourceS && (it -> getDestination()) == destinationS) {
				it -> setWeight(DEFAULT_WEIGHT); //if found, update the weight
				break;
			} 
		}
	//otherwise, do nothing
	}
}

void Graph::printGraph() {
	if (arr_edge.empty() == true) {
		std::cout << "The graph is empty" << std::endl;
	} else {
		int i = 1;
		for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
			std::cout << "Edge " << i << " source: " << it -> getSource() << " destination: " << it -> getDestination() << " weight: " << it -> getWeight() << std::endl;
			i++;
		}
		std::cout << std::endl;
	}
	
}

void Graph::Connect(int source, int destination, double weight) { //IF THE EDGE ALREADY EXISTS, THE WEIGHT WILL BE UPDATED IF PARAMETER WEIGHT IS LOWER THAN THE CURRENT WEIGHT
	if (source >= no_of_vertices || destination >= no_of_vertices) {
		//do nothing
		std::cout << "source or destination indices are out of bound" << std::endl;
	} else {

		double curweight = 0;
		bool is_exist = false;
		for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
			if ((it -> getSource()) == source && (it -> getDestination()) == destination) {
				is_exist = true; //an edge already connects the two vertices. Check the weight of the previous one, replace if the added weight is smaller
				curweight = it -> getWeight();
				if (weight < curweight) {
					it -> setWeight(weight);
				}
				break;
			}
		}
		if (is_exist == false) { //if no edge is connecting the two vertices
			arr_edge.push_back(*(new Edge(source, destination, weight)));
			no_of_edges++;
		} 
	}
}

void Graph::Connect(int source, int destination) {
	if (source >= no_of_vertices || destination >= no_of_vertices) {
		//do nothing
		std::cout << "source or destination indices are out of bound" << std::endl;
	} else {
		bool is_exist = false;
		for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
			if (it -> getSource() == source && it -> getDestination() == destination) {
				is_exist = true; //an edge already connects the two vertices. Abort and do nothing
				break;
			}
		}
		if (is_exist == false) { //if no edge is connecting the two vertices
			arr_edge.push_back(*(new Edge(source, destination)));
			no_of_edges++;
		} else { 
			//if there is already an edge connecting the two vertices, do nothing 
		}
	}
}

Edge * Graph::get_min_incoming_edge(int node) { //node is the destination. return the minimum incoming edge to "node", or null is no edge ends at the node
	int min = 999.5; //dummy value that will certainly be replaced
	Edge * result = nullptr;
	for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
		if ((it -> getDestination() == node) && (it -> getWeight() < min)) {
			result = &(*it);
			min = it -> getWeight();
		}
	}
	return result;
}

bool Graph::is_edge_ending_in(int destination) {
	for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
		if (it -> getDestination() == destination) {
			return true;
		}
	}
	return false;
}

void Graph::remove_edges_ending_in(int node) {
	std::vector<Edge>::iterator it = arr_edge.begin();
	while (it != arr_edge.end()) {
		if (it -> getDestination () == node) {
			it = arr_edge.erase(it);
		} else {
			it++;
		}
	}
}

std::vector<Edge> Graph::get_Edges_Starting_From(int node) {
	std::vector<Edge> result;
	for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
		if (it -> getSource() == node) {
			result.push_back(*it);
		}
	}
}

//use DFS to detect cycle
bool Graph::does_cycle_exist() {
	bool result = false;
	std::vector<int> visited; //vector of vertices that have been visited
	std::stack<int> frontier; //the frontier: the pointers to edges that need to be visited
	std::vector<Edge> actions;
	//now, get the first edge that we wish to explore
	Edge * first = &arr_edge.front();
	visited.push_back(first -> getSource()); //mark the starting node as already visited
	frontier.push(first -> getSource()); //put the first node into the frontier
	while (!frontier.empty()) {
		if (visited.size() == no_of_vertices) {
			break; //if we've visited all the vertices, thus we've iterated over the whole tree
		}
		for (std::vector<Edge>::iterator e = )
		//actions = get_Edges_Starting_From(frontier.top()); //get all the possible actions from the current top of the stack
		std::cout << "Size of actions: " << actions.size() << std::endl;
		frontier.pop(); //remove the top of the stack
		for (std::vector<Edge>::iterator e = actions.begin(); e != actions.end(); e++) {
			if (does_exist<int> (e -> getDestination(), visited)) {
				result = true;
				break;
			}
			frontier.push(e -> getDestination());
			visited.push_back(e -> getDestination());
		}
		if (result == true) {
			break;
		}
		actions.clear();
	}
	return result;
}


//remove all the edges that have "node" as their destinations
/*
void Graph::remove_edges_ending_in(int node) {
	for (std::vector<Edge>::iterator it = arr_edge.begin(); it != arr_edge.end();) {
		if (it -> getDestination() == node) {
			arr_edge.erase(it);
		} else {
			it++;
		}
	}
} */



