#include "graph.cpp"

int main() {
	Graph G(3, 0);
	G.Connect(0, 1, 100);
	G.Connect(0, 2, 90);
	G.Connect(1, 0, 70);
	G.Connect(1, 2, 20);
	G.Connect(2, 0, 20);
	G.Connect(2, 1, 30);
	G.printGraph();
	std::vector<Edge> min_incoming_edges = G.get_min_incoming_edges();
	for (std::vector<Edge>::iterator it = min_incoming_edges.begin(); it != min_incoming_edges.end(); it++) {
		std::cout << "This is an edge from " << it->getSource() << " to " << it->getDestination() << std::endl;
	}
	if (detect_cycle(min_incoming_edges)) {
		std::cout << "There is a cycle in the min incoming edges" << std::endl;
		std::vector<int> the_cycle = get_cycle(min_incoming_edges);
		std::cout << "The cycle is " << the_cycle[0] << " and " << the_cycle[1] << std::endl;
	} else { 
		std::cout << "There isn't a cycle in the min incoming edges" << std::endl;
	}
	std::cout << "Applying MST algorithm" << std::endl;
	std::vector<Edge> MST = G.get_MST();
	for (std::vector<Edge>::iterator it = MST.begin(); it != MST.end(); it++) {
		std::cout << "This is an edge from " << it->getSource() << " to " << it->getDestination() << std::endl;
	}
	return 0;
}