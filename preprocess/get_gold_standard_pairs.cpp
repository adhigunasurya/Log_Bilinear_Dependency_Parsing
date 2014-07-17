#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
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

template <typename T>
T StringToNumber ( const std::string &Text )//Text not by const reference so that the function can be used with a 
{                               //character array as argument
    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

int main() {
    std::ifstream myfile("pairs_dev.txt");
    std::string line;
    std::ofstream output("pairs_dev_gold_standard.txt");
    if (myfile.is_open() && output.is_open()) {
        int curr_index = 1;
        while (getline (myfile,line)) {
            if (line.size() > 1) {
                std::vector<std::string> temp = split(line, ' ');
                output << curr_index << ' ' << curr_index - StringToNumber<int>(temp[6]) << std::endl;
                curr_index++;
            } else {
                curr_index = 1;
                output << std::endl;
            }
        }
        myfile.close();
        output.close();
    } else {
        std::cout << "Either the input or the output file fails to open" << std::endl;
    }

    return 0;
}
