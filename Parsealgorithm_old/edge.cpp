#include <iostream>

#include "edge.h"

//default constructor
Edge::Edge (int sourceS, int destinationS, double weightS) {
	this->source = sourceS;
	this->destination = destinationS;
	this->weight = weightS;
}

//constructor with missing weight, assign value of default weight
Edge::Edge(int sourceS, int destinationS) {
	this->source = sourceS;
	this->destination = destinationS;
	this->weight = DEFAULT_WEIGHT;
}

Edge::Edge(const Edge& ref) {
	this->source = ref.source;
	this->destination = ref.destination;
	this->weight = ref.weight;
}


double Edge::getWeight() {
	return this->weight;
}

int Edge::getSource() {
	return this->source;
}

int Edge::getDestination() {
	return this->destination;
}

void Edge::setWeight(const double weightS) {
	this->weight = weightS;
}

void Edge::setSource(int sourceS) {
	this->source = sourceS;
}

void Edge::setDestination(int destinationS) {
	this->destination = destinationS;
}

void Edge::printEdge() {
	std::cout << "This is an edge connecting vertex " << source << " to vertex " << destination << " with weight " << weight << std::endl;
}

inline bool Edge::operator==(const Edge& comparison) {
	if (source == comparison.source && destination == comparison.destination && weight == comparison.weight) {
		return true;
	} else {
		return false;
	}
}