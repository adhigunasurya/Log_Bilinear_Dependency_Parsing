#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <iostream>
#include <stdlib.h>

#ifndef DUMMY_HEAD
#define DUMMY_HEAD -1
#endif

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
	std::string line;
	std::vector<std::string> temp;
	std::vector<std::string> pos_tags; //represent the sentence as a vector of fine-grained pos tags
	std::vector<int> head;
	std::ifstream myfile ("english_wsj_train_edited.conll");
	std::ofstream output_file("output_file.txt");
	int result = -1;
	if (!output_file.is_open()) {
		std::cout << "Fail to open output file" << std::endl;
	}
	bool first = true;
	//get the length of sentences
	if (myfile.is_open()) {
	    while (getline (myfile,line))
	    {
	      if (line.size() > 1) { //if the line is not empty
	      		temp = split(line, '\t'); //split the string with tab (\t) as its delimiter
	      		const char * c = temp[0].c_str();
	      		result = atoi(c);
	      		if (result == 1) { //this indicates the beginning of a new sentence
	      			//wrap up the previous vector of pos tags by adding "END"
	      			if (first == true) {
	      				first = false;
	      			} else {
	      				pos_tags.push_back("END");
	      				head.push_back(DUMMY_HEAD);
	      				//write the vector into an external file
		      			for (std::vector<std::string>::iterator it = pos_tags.begin(); it != pos_tags.end(); it++) {
		      				output_file << *it << ' ';
		      			}
		      			output_file << '\n';
		      			for (std::vector<int>::iterator it = head.begin(); it != head.end(); it++) {
		      				output_file << *it << ' ';
		      			}
		      			output_file << '\n';
	      			}
	      			//clear the vector to prepare it for the new sentence
	      			pos_tags.clear();
	      			head.clear();
	      			pos_tags.push_back("START");
	      			head.push_back(DUMMY_HEAD);
	      		} 
	      		pos_tags.push_back(temp[4]);
	      		const char * ch = temp[6].c_str();
	      		head.push_back(atoi(ch));
	      }
	    }
	    pos_tags.push_back("END");
	    head.push_back(DUMMY_HEAD);
	    for (std::vector<std::string>::iterator it = pos_tags.begin(); it != pos_tags.end(); it++) {
		    output_file << *it << ' ';
		}
		output_file << '\n';
		for (std::vector<int>::iterator it = head.begin(); it != head.end(); it++) {
			output_file << *it << ' ';
		}
		output_file << '\n';
	    myfile.close();
	    output_file.close();
	 } else {
  		std::cout << "Unable to open file"; 
  	 }
	return 0;
}