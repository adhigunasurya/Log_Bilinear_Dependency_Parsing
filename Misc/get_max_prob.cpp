#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <random>
#include <stdlib.h>
#include <cmath>
#include <utility>
#include <ctime>

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

template <typename T>
T StringToNumber (const std::string &Text )//Text not by const reference so that the function can be used with a 
{                               //character array as argument
	std::stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

template <typename T>
std::string NumberToString ( T Number )
{
	std::stringstream ss;
	ss << Number;
	return ss.str();
}


int main() {
	std::ifstream prob("prob_result.txt");
	std::map<std::string, double> postag_frequency_map;
	std::map<std::string, std::string> postag_postag_map;
	std::ifstream fine_pos_tags("fine_pos_tags.txt");
	std::string line;
	//initialise
	while (getline(fine_pos_tags, line)) {
		postag_frequency_map[line] = -9999.0;
	}
	fine_pos_tags.close();
	//get the frequencies
	while (getline(prob, line)) {
		std::vector<std::string> temp = split(line, ' ');
		if (StringToNumber<double> (temp[4]) > postag_frequency_map.at(temp[0])) {
			postag_frequency_map[temp[0]] = StringToNumber<double> (temp[4]);
			postag_postag_map[temp[0]] = temp[2];
		}
	}
	std::ofstream check_result ("check_result.txt");
	for (std::map<std::string, std::string>::iterator it = postag_postag_map.begin(); it != postag_postag_map.end(); it++) {
		check_result << it->first << " | " << it->second << std::endl;
	}
	prob.close();
	return 0;
}