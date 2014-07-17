#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

template <class T>
bool contains(T element, std::vector<T> vect) {
	for (typename std::vector<T>::iterator it = vect.begin(); it != vect.end(); it++) {
		if (element == *it) {
			return true;
		}
	}
	return false;
}

int main() {
	std::ifstream pairs("pairs.txt");
	std::ofstream output("only_postag_bias.txt");
	if (!pairs.is_open() || !output.is_open()) {
		std::cout << "The files cannot open" << std::endl;
	} else {
		std::map<std::string, int> mapping;
		std::string line;
		while (getline(pairs, line)) {
			if (line.length() > 1) {
				std::vector<std::string> temp = split(line, ' ');
				if (mapping.find(temp[0]) == mapping.end()) {
					mapping[temp[0]] = 1;
				} else {
					mapping[temp[0]] = mapping[temp[0]] + 1;
				}
			}
		}
		for (std::map<std::string, int>::iterator it = mapping.begin(); it != mapping.end(); it++) {
			output << it->first << " " << it->second << std::endl;
		}
	}

	return 0;
}