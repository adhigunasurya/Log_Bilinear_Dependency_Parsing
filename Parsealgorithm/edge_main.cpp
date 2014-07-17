#include "edge.cpp"

int main() {
	Edge * Edge1 = new Edge(1, 7, 9.0);
	Edge Edge2(2, 5, 7.5);
	Edge * Edge3 = new Edge(2, 8);
	Edge Edge4(*Edge1);
	std::cout << "Edge 1: " << Edge1->getSource() << " " << Edge1->getDestination() << " " << Edge1->getWeight() << std::endl;
	std::cout << "Edge 2: " << Edge2.getSource() << " " << Edge2.getDestination() << " " << Edge2.getWeight() << std::endl;
	std::cout << "Edge 3: " << Edge3->getSource() << " " << Edge3->getDestination() << " " << Edge3->getWeight() << std::endl;
	std::cout << "Edge 4: " << Edge4.getSource() << " " << Edge4.getDestination() << " " << Edge4.getWeight() << std::endl;
	Edge4 = *Edge3;
	std::cout << "Edge 4: " << Edge4.getSource() << " " << Edge4.getDestination() << " " << Edge4.getWeight() << std::endl;
	Edge1->setWeight(5.728123);
	Edge1->setDestination(10);
	Edge1->setSource(8);
	std::cout << "Edge 1: " << Edge1->getSource() << " " << Edge1->getDestination() << " " << Edge1->getWeight() << std::endl;
	return 0;
}