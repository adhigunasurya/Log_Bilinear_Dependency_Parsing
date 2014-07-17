#include "graph.cpp"
#include <ctime>

int main() {
	/* test to check if cycle exists
	Graph gr2(1);
	gr2.printGraph();
	gr2.Disconnect(1, 0);
	if (gr2.does_cycle_exist()) {
		std::cout << "A cycle exists" << std::endl; 
		//print out the pair of vertices that cause the cycle
		std::cout << "The pair is " << std::get<0> (gr2.get_First_Cycle()) << " " << std::get<1> (gr2.get_First_Cycle()) << std::endl;
	} else {
		std::cout << "No cycle exists" << std::endl;
	} */

	/* test the find MST function, without cycle 
	Graph gr3(2);
	gr3.printGraph();
	gr3.Connect(0, 1, -log(0.2));
	gr3.Connect(0, 2, -log(0.5));
	gr3.Connect(1, 0, -log(0.05));
	gr3.Connect(1, 2, -log(0.1));
	gr3.Connect(2, 1, -log(0.4));
	gr3.Connect(2, 0, -log(0.04));
	gr3.printGraph();
	std::vector<Edge> result = gr3.getMST(0);
	for (std::vector<Edge>::iterator it = result.begin(); it != result.end(); it++) {
		it -> printEdge();
	} */ 
	Graph gr3(50); //BIG EXAMPLE
	gr3.Connect(0, 3, -log(20));
	gr3.Connect(3, 2, -log(18));
	gr3.Connect(2, 1, -log(15));
	gr3.Connect(3, 5, -log(19));
	gr3.Connect(5, 4, -log(20));
	gr3.Connect(3, 6, -log(21));
	gr3.Connect(3, 9, -log(23));
	gr3.Connect(9, 8, -log(20));
	gr3.Connect(8, 7, -log(20));
	gr3.Connect(9, 10, -log(10));
	//gr3.printGraph();
	clock_t t1, t2;
	t1 = clock();
	std::vector<std::string> result = gr3.getPairs(0);
	for (std::vector<std::string>::iterator it = result.begin(); it != result.end(); it++) {
		std::cout << *it << std::endl;
	} 
	t2 = clock();
	double diff = ((double) t2 - (double) t1);
	//std::vector<Edge> result = gr3.getMST(0);
	
	/*std::cout<<"total program running time: "<<diff/CLOCKS_PER_SEC<<std::endl;
	for (std::vector<Edge>::iterator it = result.begin(); it != result.end(); it++) {
		it -> printEdge();
	} */

//small example
	/*Graph gr3(3);
	gr3.Connect(0, 1, -log(9));
	gr3.Connect(0, 2, -log(10));
	gr3.Connect(0, 3, -log(9));
	gr3.Connect(1, 2, -log(20));
	gr3.Connect(1, 3, -log(3));
	gr3.Connect(1, 0, -log(0.2));
	gr3.Connect(2, 1, -log(30));
	gr3.Connect(2, 3, -log(30));
	gr3.Connect(2, 0, -log(0.01));
	gr3.Connect(3, 1, -log(11));
	gr3.Connect(3, 2, -log(0.00001));
	gr3.Connect(3, 0, -log(0.0928));
	gr3.printGraph();
	
	clock_t t1, t2;
	t1 = clock();
	std::vector<Edge> result = gr3.getMST(0);
	t2 = clock();
	double diff = ((double) t2 - (double) t1);
	std::cout<<"total program running time: "<<diff<<std::endl;
	for (std::vector<Edge>::iterator it = result.begin(); it != result.end(); it++) {
		it -> printEdge();
	}
	gr3.printGraph();
 */
	/* test to check cycles with 4 vertices
	Graph gr1(3); //test constructor
	gr1.printGraph();
	//gr1.Disconnect(0, 2);
	gr1.Disconnect(0, 3);
	gr1.Disconnect(2, 1);
	gr1.Disconnect(1, 0);
	gr1.Disconnect(2, 0);
	gr1.Disconnect(1, 3);
	gr1.Disconnect(3, 0);
	gr1.Disconnect(3, 1);
	gr1.Disconnect(3, 2);
	gr1.printGraph();
	
	if (gr1.does_cycle_exist()) {
		std::cout << "A cycle exists" << std::endl; 
		//print out the pair of vertices that cause the cycle
		std::cout << "The pair is " << std::get<0> (gr1.get_First_Cycle()) << " " << std::get<1> (gr1.get_First_Cycle()) << std::endl;
	} else {
		std::cout << "No cycle exists" << std::endl;
	}  */
	/* test to remove the edges
	gr1.Disconnect(0, 2);
	gr1.Disconnect(2, 0);
	gr1.printGraph();
	std::cout << "number of vertex is: " << gr1.vertexCount() << " and number of edge is: " << gr1.edgeCount() << std::endl;
	gr1.Connect(0, 2, -2.5);
	gr1.printGraph();	
	gr1.Connect(2, 0, -3.7);
	gr1.Connect(0, 1, -10.7);
	gr1.Connect(0, 2, -2.0);
	gr1.Connect(1, 0, -8.9);
	gr1.Connect(1, 2, -6.4);
	gr1.Connect(2, 1, -4.77);
	gr1.Connect(3, 1, -4.77);
	gr1.printGraph();
	gr1.remove_edges_ending_in(0);
	gr1.printGraph();
	gr1.remove_edges_ending_in(1);
	gr1.printGraph();
	gr1.remove_edges_ending_in(2);
	gr1.printGraph();
	*/
	return 0;
}