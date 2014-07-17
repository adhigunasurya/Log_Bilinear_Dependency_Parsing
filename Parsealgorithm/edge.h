#ifndef EDGE_H
#define EDGE_H

#ifndef DEFAULT_WEIGHT
#define DEFAULT_WEIGHT 99999.5 //set the default weight for edges connecting unconnected vertices
#endif

#include <iostream>

//avoid using namespace std, especially in headers, since all implementations will be using that

class Edge {
private:
	int source = 0;
	int destination = 0;
	double weight = 0.0;

public:
	//declare the constructors
	Edge() {source = 0; destination = 0; weight = 0.0;} //default constructor, do nothing
	Edge(int sourceS, int destinationS, double weightS); //constructor, set everything the edge needs to store
	Edge(int sourceS, int destinationS); //constructor, called to build originally non-existing edges (assign a default weight), to make the graph strongly connected
	//declare the copy constructor
	Edge(const Edge& ref);
	//declare the copy assignment
	Edge& operator= (const Edge& ref) {	this->source = ref.source;
	this->destination = ref.destination;
	this->weight = ref.weight;}
	//declare all the getters
	double getWeight(); //get the weight of the edge
	int getSource(); //get the source of the edge
	int getDestination(); //get the destination of the edge
	//declare all the setters
	void setWeight(double weightS);
	void setSource(int sourceS);
	void setDestination(int destinationS);
	void printEdge();
	inline bool operator==(const Edge& comparison);
};

#endif //EDGE_H