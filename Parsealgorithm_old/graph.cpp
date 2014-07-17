#include "graph.h"
#include <exception> //exception is important since we are allocating memories dynamically where an exception may be thrown
#include <stack>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


//function to calculate resulting number of edges
template <class T> 
bool does_exist(T element, std::vector<T>& arr) {
	for (typename std::vector<T>::iterator it = arr.begin() ; it != arr.end(); it++) {
		if (*it == element) {
			return true;
		}
	}
	return false;
}

bool does_edge_exist(int source, int destination, std::vector<Edge>& arr) {
	for (std::vector<Edge>::iterator it = arr.begin(); it != arr.end(); it++) {
		if (it->getSource() == source && it->getDestination() == destination) {
			return true;
		}
	}
	return false;
}

//get the (unique) edge that starts at source and ends at destination out of a vector of Edges called arr
Edge * get_Edge_Vector (int source, int destination, std::vector<Edge>& arr) {
	for (std::vector<Edge>::iterator it = arr.begin(); it != arr.end(); it++) {
		if (it->getSource() == source && it->getDestination() == destination) {
			return &(*it);
		}
	}
	return nullptr;
}

void delete_Edge(int source, int destination, std::vector<Edge>& arr) {
	for (std::vector<Edge>::iterator it = arr.begin(); it != arr.end(); it++) {
		if (it->getSource() == source && it->getDestination() == destination) {
			arr.erase(it);
		}
	}
}


int n_of_edges(int n) { //return the number of edges that will be in the directed graph
	return ((n+1) * n);
}

int Graph::edgeCount() {
	return arr_edge.size();
}

int Graph::vertexCount() {
	return no_of_vertices;
}

//sole constructor
Graph::Graph(int n) { //IMPORTANT: n is the number of tokens, thus the size of the graph (number of vertices in the graph) will be n+1
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

//recall that n_v is the number of vertices (not the number of tokens in a sentence such as in the previous case)
Graph::Graph(int n_v, std::vector<Edge>& vect) {
	no_of_vertices = n_v;
	for (std::vector<Edge>::iterator it = vect.begin(); it != vect.end(); it++) {
		arr_edge.push_back(*it); //deep copy
	}
}

Graph::Graph(const Graph& ref) { //cctor
	//set the size
	this -> no_of_vertices = ref.no_of_vertices;
	//(this -> arr_edge).clear();
	arr_edge = *(new std::vector<Edge> (ref.arr_edge)); //copy the vector of edges that belong to ref (deep copy)
}

Graph& Graph::operator=(const Graph& ref) { //copy assignment
	this -> no_of_vertices = ref.no_of_vertices;
	//delete all the elements from the object's array of edges
	(this -> arr_edge).clear();
	//initialise the array of edges using copy constructor with ref.arr_edge as its sole parameter (deep copy)
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

Edge * Graph::get_Edge(int sourceS, int destinationS) { //get the edge connecting the source and the destination, return null_ptr if none exists {
	if (sourceS >= no_of_vertices || destinationS >= no_of_vertices) {
		//do nothing if the source and destination parameters are beyond the valid index
		std::cout << "source or destination indices are out of bound" << std::endl;
		return nullptr; //there is no edge connecting them as the indices are invalid
	}
	for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
		if ((it -> getSource()) == sourceS && (it -> getDestination()) == destinationS) {
			return &(*it); // return the edge connecting the source and the destination
		} 
	}
	return nullptr; 
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
	} else if (!arr_edge.empty()){
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
		if (is_exist == false && source != destination) { //if no edge is connecting the two vertices and the edge is not connecting a vertex to itself
			arr_edge.push_back(*(new Edge(source, destination, weight)));
		} 
	} else if (arr_edge.empty()){
		if (source != destination) {
			arr_edge.push_back(*(new Edge(source, destination, weight)));
		}
	}
}

void Graph::Connect(int source, int destination) {
	if (source >= no_of_vertices || destination >= no_of_vertices) {
		//do nothing
		std::cout << "source or destination indices are out of bound" << std::endl;
	} else if (!arr_edge.empty()) {
		bool is_exist = false;
		for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
			if (it -> getSource() == source && it -> getDestination() == destination) {
				is_exist = true; //an edge already connects the two vertices. Abort and do nothing
				break;
			}
		}
		if (is_exist == false && source != destination) { //if no edge is connecting the two vertices and the edge is not connecting a vertex to itself
			arr_edge.push_back(*(new Edge(source, destination)));
		} else { 
			//if there is already an edge connecting the two vertices, do nothing 
		}
	} else if (arr_edge.empty()) {
		if (source != destination) {
			arr_edge.push_back(*(new Edge(source, destination)));
		}
	}
}

Edge * Graph::get_min_incoming_edge(int node) { //node is the destination. return the minimum incoming edge to "node", or null is no edge ends at the node
	int min = 99999.5; //dummy value that will certainly be replaced
	Edge * result = nullptr;
	for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
		if ((it -> getDestination() == node) && (it -> getWeight() < min)) {
			result = &(*it);
			min = it -> getWeight();
		}
	}
	return result;
}

std::vector<Edge> Graph::get_min_incoming_edge() {
	std::vector<int> min_temp(no_of_vertices, 9999.0);
	std::vector<Edge*> edge_temp(no_of_vertices, nullptr); //each vertex would have at most one, except the root which has none
	std::vector<Edge> result;
	for (std::vector<Edge>::iterator it = arr_edge.begin() ; it != arr_edge.end(); it++) {
			if (it->getWeight() < min_temp.at(it->getDestination())) {
				min_temp[it->getDestination()] = it->getWeight();
				edge_temp[it->getDestination()] = &(*it);
			}
	}
	for (int n = 0; n < no_of_vertices; n++) {
		if (edge_temp[n] != nullptr) {
			result.push_back(*(edge_temp[n]));
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
	return result;
}

//use DFS to detect cycle
bool Graph::does_cycle_exist() {
	bool result = false;
	std::vector<int> visited; //vector of vertices that have been visited
	std::stack<int> frontier; //the frontier: the pointers to edges that need to be visited
	//std::vector<Edge> actions;
	std::vector<Edge> actions;
	//now, get the first edge that we wish to explore (arbitrarily)
	Edge * first = &arr_edge.front();
	frontier.push(first -> getSource()); //put the first node into the frontier
	while (!frontier.empty()) {
		if (visited.size() == no_of_vertices) {
			break; //if we've visited all the vertices, thus we've iterated over the whole tree
		}
		//First, put the current node being examined into the list of visited nodes
		visited.push_back(frontier.top());
		//get all possible "actions" (obtain all edges that start from the current node being examined)
		actions = get_Edges_Starting_From(frontier.top()); //get all the possible actions from the current top of the stack
		//we thus no longer need the top of the stack
		frontier.pop(); 
		// iterate over all possible actions, add the "destination node" of each edge into the frontier to be examined later
		for (std::vector<Edge>::iterator e = actions.begin(); e != actions.end(); e++) {
			//if the destination node has already been visited before, then we detect a cycle and terminate the whole loop to save time
			if (does_exist<int> (e -> getDestination(), visited)) {
				result = true;
				break;
			}
			frontier.push(e -> getDestination());
		}
		if (result == true) {
			break;
		}
		actions.clear();
	}
	return result;
}

bool Graph::does_cycle_exist(int root) {
	bool result = false;
	std::vector<int> visited; //vector of vertices that have been visited
	std::stack<int> frontier; //the frontier: the pointers to edges that need to be visited
	//std::vector<Edge> actions;
	std::vector<Edge> actions;
	//now, get the first edge that we wish to explore, in this case the first one to be inspected is the root
	frontier.push(root); //put the root into the frontier
	while (!frontier.empty()) {
		if (visited.size() == no_of_vertices) {
			break; //if we've visited all the vertices, thus we've iterated over the whole tree
		}
		//First, put the current node being examined into the list of visited nodes
		visited.push_back(frontier.top());
		//get all possible "actions" (obtain all edges that start from the current node being examined)
		actions = get_Edges_Starting_From(frontier.top()); //get all the possible actions from the current top of the stack
		//we thus no longer need the top of the stack
		frontier.pop(); 
		// iterate over all possible actions, add the "destination node" of each edge into the frontier to be examined later
		for (std::vector<Edge>::iterator e = actions.begin(); e != actions.end(); e++) {
			//if the destination node has already been visited before, then we detect a cycle and terminate the whole loop to save time
			if (does_exist<int> (e -> getDestination(), visited)) {
				result = true;
				break;
			}
			frontier.push(e -> getDestination());
		}
		if (result == true) {
			break;
		}
		actions.clear(); //delete the current actions to save memory
	}
	return result;
}

//get the first cycle. IMPORTANT: Only call whenever you are sure that the graph does have a cycle (using the previous function)
std::pair<int, int> Graph::get_First_Cycle() {
	bool cycle_identified = false; //indicates whether a cycle has been identified. Used to break from the loop once we've identified a cycle. 
	int prev_node = -1; //dummy previous node
	std::pair<int, int> result_pair;
	std::vector<int> visited; //vector of vertices that have been visited
	std::stack<int> frontier; //the frontier: the pointers to edges that need to be visited
	//std::vector<Edge> actions;
	std::vector<Edge> actions;
	//now, get the first edge that we wish to explore (arbitrarily)
	Edge * first = &arr_edge.front();
	frontier.push(first -> getSource()); //put the first node into the frontier
	while (!frontier.empty()) {
		if (visited.size() == no_of_vertices) {
			break; //if we've visited all the vertices, thus we've iterated over the whole tree
		}
		//First, put the current node being examined into the list of visited nodes
		visited.push_back(frontier.top());
		//get all possible "actions" (obtain all edges that start from the current node being examined)
		actions = get_Edges_Starting_From(frontier.top()); //get all the possible actions from the current top of the frontier
		//we thus no longer need the top of the stack, but keep its value with variable prev_node
		prev_node = frontier.top();
		frontier.pop(); 
		// iterate over all possible actions, add the "destination node" of each edge into the frontier to be examined later
		for (std::vector<Edge>::iterator e = actions.begin(); e != actions.end(); e++) {
			//if the destination node has already been visited before, then we detect a cycle and terminate the whole loop to save time
			if (does_exist<int> (e -> getDestination(), visited)) {
				cycle_identified = true;
				result_pair = std::make_pair (prev_node, e -> getDestination());
				break;
			}
			frontier.push(e -> getDestination());
		}
		if (cycle_identified == true) { //we only need one pair of nodes that form a cycle
			break;
		}
		actions.clear();
	}
	return result_pair;
}


Edge * Graph::get_Edge_By_Weight(double weight, Edge& alternative_1, Edge& alternative_2) { // pick which edge (between alternative 1 and 2) has the same weight as the specified weight
	Edge * result;
	if (alternative_1.getWeight() == weight) {
		result = &alternative_1;
	} else if (alternative_2.getWeight() == weight) {
		result = &alternative_2;
	} else {
		return nullptr;
	}
}

//NOTE: by calling this function, the graph is modified (all connections to the root are disabled)
std::vector<Edge> Graph::getMST(int root) {
	Graph copy_graph(*this);
	std::vector<Edge> result;
	std::vector<Edge> lowest_incoming_edges;
	//first, remove all parallel edges with a single one (unimplemented)
	//second, remove all edges ending at the root 
	copy_graph.remove_edges_ending_in(root);
	//third, get a single lowest incoming edge for each vertex other than the root
	lowest_incoming_edges = copy_graph.get_min_incoming_edge();
	//delete the lowest incoming edge that ends at the root
	//Edge * temporary_edge = get_min_incoming_edge(root);
	/*for (int i = 0; i < copy_graph.no_of_vertices; i++) {
		if (i != root) {
			if (copy_graph.get_min_incoming_edge(i) != nullptr) {
				lowest_incoming_edges.push_back(*(copy_graph.get_min_incoming_edge(i)));
			}
			
		}	
	} */
	//fourth, construct a graph consisting of all the vertices and only the lowest incoming edges for all vertices except the root
	Graph lowest_edges_graph(no_of_vertices, lowest_incoming_edges);
	//lowest_edges_graph.printGraph();
	//if the graph with only the lowest incoming edges have no cycle, then the process terminates and the lowest incoming edges is the result
	//check for cycles
	bool temp = false;
	for (int n=0; n < no_of_vertices; n++) {
		if (lowest_edges_graph.does_cycle_exist(n)) {
			temp = true;
			break;
		}
	}

	//there is no cycle, then result is simply those lowest incoming edges
	if (temp == false) {
		result = lowest_incoming_edges;
		lowest_incoming_edges.clear();
	}
	/*if (!lowest_edges_graph.does_cycle_exist()) {
		result = lowest_incoming_edges;
		lowest_incoming_edges.clear(); */
	else { //otherwise, a cycle exists. Identify the cycle
		std::pair<int, int> cycle = lowest_edges_graph.get_First_Cycle();
		//define a new directed graph where the two vertices that form a cycle is "condensed" into one vertex
		/*Recreate the same graph, BUT with the two vertices that form a cycle as declared above "condensed together" in the internal representation of the graph 
		Recall that the parameter for the constructor of the graph is the number of tokens, which is number of vertices - 2
		This is because the new graph would then have precisely n-1 vertices, which is exactly what we want (condensed)*/
		Graph temp_graph(no_of_vertices - 2); 
		//clear all the edges from temp_graph to start from the beginning
		temp_graph.arr_edge.clear();
		//create a mapping for the vertices (such as vertex 1 in the original graph is equivalent with vertex 2 in the temp graph)
		//IMPLEMENT MAPPING
		std::map<int, int> mapping;
		int j = 0;
		for (int i = 0; i < copy_graph.no_of_vertices; i++) {
			if (i != std::get<0> (cycle) && i != std::get<1> (cycle)) { //if the current vertex is not one of the nodes that form a cycle
				mapping[i] = j;
				j++;
			} else {
				mapping[i] = temp_graph.no_of_vertices - 1; //reserve the last vertex of temp_graph for the "condensed vertex"
			}
		}
		//iterate over all the edges in the original graph, add the edges according to the rules for the Edmond's algorithm page in Wikipedia
		int new_idx1 = 0;
		int new_idx2 = 0;
		double temp_weight = 0.0;
		for (std::vector<Edge>::iterator it = copy_graph.arr_edge.begin(); it != copy_graph.arr_edge.end(); it++) {
			new_idx1 = mapping.find(it->getSource())->second;
			new_idx2 = mapping.find(it->getDestination())->second;
			//case 1: both vertices connected by the particular edge are not in the cycle
			if (it->getSource() != std::get<0>(cycle) && it->getSource() != std::get<1>(cycle) && it->getDestination() != std::get<1>(cycle) && it->getDestination() != std::get<0>(cycle)) {
				temp_graph.Connect(new_idx1, new_idx2, it->getWeight()); //add the edge, including its weight, (and thus both vertices) into the condensed graph
				//case 2: the source vertex belongs to the cycle but the target vertex doesn't
			} else if ((it->getSource() == std::get<0>(cycle) || it->getSource() == std::get<1>(cycle)) && it->getDestination() != std::get<0>(cycle) && it->getDestination() != std::get<1>(cycle)) { 
				temp_graph.Connect(new_idx1, new_idx2, it->getWeight());
				//case 3: the source vertex does not belong to the cycle but the target vertex does
			} else if (it->getSource() != std::get<0>(cycle) && it->getSource() != std::get<1>(cycle) && (it->getDestination() == std::get<0>(cycle) || it->getDestination() == std::get<1>(cycle))) {
				//as defined in Wikipedia, the weight is then: the weight of the edge - the weight of the lowest-weighting edge that ends at the destination vertex
				temp_weight = it->getWeight() - copy_graph.get_min_incoming_edge(it->getDestination())->getWeight();
				temp_graph.Connect(new_idx1, new_idx2, temp_weight);
			}
		}
			//perform a recursive call to find the minimum spanning tree of temp_graph
			std::vector<Edge> temp_MST = temp_graph.getMST(root);
			//First, include all the edges in the cycle:
			Edge * temp_edge = copy_graph.get_Edge(std::get<0>(cycle), std::get<1>(cycle));
			if (temp_edge != nullptr) {
				result.push_back(*temp_edge);
			}
			temp_edge = copy_graph.get_Edge(std::get<1>(cycle), std::get<0>(cycle));
			if (temp_edge != nullptr) {
				result.push_back(*temp_edge);
			}
			//second, mark all the "common edges" between copy_graph and temp_graph
			int temp_idx1;
			int temp_idx2;
			Edge * tmp_edge;
			for (std::vector<Edge>::iterator it = copy_graph.arr_edge.begin(); it != copy_graph.arr_edge.end(); it++) {
				temp_idx1 = mapping.find(it->getSource())->second;
				temp_idx2 = mapping.find(it->getDestination())->second;
				if (does_edge_exist(temp_idx1, temp_idx2, temp_MST)) { //this is a common edge between the minimum incoming edges in temp_MST and the original graph
					// get the particular edge that is being examined
					tmp_edge = get_Edge_Vector(temp_idx1, temp_idx2, temp_MST);
					//case 1: both temp_idx1 and temp_idx2 are not the "condensed" node
					if (temp_idx1 != (temp_graph.no_of_vertices - 1) && temp_idx2 != (temp_graph.no_of_vertices - 1)) {
						//then simply push the edge into the result
						result.push_back(*it);
						//case 2: the source is the condensed node but the destination is not
					} else if (temp_idx1 == (temp_graph.no_of_vertices - 1) && temp_idx2 != (temp_graph.no_of_vertices - 1)) {
						if (it->getWeight() == tmp_edge->getWeight()) { 
							result.push_back(*it);
						}
						//case 3: the source is not the condensed node but the destination is the condesed node
					} else if (temp_idx1 != (temp_graph.no_of_vertices - 1) && temp_idx2 == (temp_graph.no_of_vertices - 1)) {
						if (tmp_edge->getWeight() == (it->getWeight() - copy_graph.get_min_incoming_edge(it->getDestination())->getWeight())) {
							result.push_back(*it);
							delete_Edge(copy_graph.get_min_incoming_edge(it->getDestination())->getSource(), it->getDestination(), result);
						}
					}
				}
			}
		lowest_incoming_edges.clear();
	}
	return result;
}

std::vector<std::string> Graph::getPairs(int root) {
	std::vector<Edge> temp = getMST(root);
	std::vector<std::string> result;
	for (std::vector<Edge>::iterator it = temp.begin(); it != temp.end(); it++) {
		std::ostringstream s1;
		s1 << it->getSource();
		std::ostringstream s2;
		s2 << it->getDestination();
		result.push_back(s1.str() + " " + s2.str());
	}
	return result;
}



