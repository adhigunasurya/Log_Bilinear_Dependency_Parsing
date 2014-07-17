#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <iostream>
#include <stdlib.h>

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
	std::vector<std::string> arr_of_coarse_pos_tags;
	std::vector<std::string> arr_of_fine_pos_tags;
	std::ifstream myfile ("english_wsj_train.conll");
	std::string coarse_pos_tag;
	std::string fine_pos_tag;
	//get the length of sentences
	if (myfile.is_open()) {
	    while (getline (myfile,line))
	    {
	      if (line.size() > 1) { //if the line is not empty
	      		temp = split(line, '\t'); //split the string with tab (\t) as its delimiter
	      		coarse_pos_tag = temp[3];
	      		fine_pos_tag = temp[4];
	      		if (!contains<std::string> (coarse_pos_tag, arr_of_coarse_pos_tags)) {
	      			arr_of_coarse_pos_tags.push_back(coarse_pos_tag);
	      		}
	      		if (!contains<std::string> (fine_pos_tag, arr_of_fine_pos_tags)) {
	      			arr_of_fine_pos_tags.push_back(fine_pos_tag);
	      		}
	      		temp.clear();
	      }
	    }
	    myfile.close();
	    myfile.open("english_wsj_dev.conll");
	  	 if (myfile.is_open()) {
	  	 	while (getline (myfile,line))
		    {
		      if (line.size() > 1) { //if the line is not empty
		      		temp = split(line, '\t'); //split the string with tab (\t) as its delimiter
		      		coarse_pos_tag = temp[3];
		      		fine_pos_tag = temp[4];
		      		if (!contains<std::string> (coarse_pos_tag, arr_of_coarse_pos_tags)) {
		      			arr_of_coarse_pos_tags.push_back(coarse_pos_tag);
		      		}
		      		if (!contains<std::string> (fine_pos_tag, arr_of_fine_pos_tags)) {
		      			arr_of_fine_pos_tags.push_back(fine_pos_tag);
		      		}
		      		temp.clear();
		      }
		    }
		    myfile.close();
	  	 } else {
	  		std::cout << "Unable to open file"; 
	  	 }
	    std::ofstream myfile_2("coarse_pos_tags.txt");
	    std::ofstream myfile_3("fine_pos_tags.txt");
	    if (myfile_2.is_open() && myfile_3.is_open()) {
	    	for (std::vector<std::string>::iterator it = arr_of_coarse_pos_tags.begin(); it != arr_of_coarse_pos_tags.end(); it++) {
	    		myfile_2 << *it << '\n';
	    	}
	    	for (std::vector<std::string>::iterator it = arr_of_fine_pos_tags.begin(); it != arr_of_fine_pos_tags.end(); it++) {
	    		myfile_3 << *it << '\n';
	    	}
	    myfile_2.close();
	    myfile_3.close();
	    } else {
	    	std::cout << "Unable to create output files" << std::endl;
	    }
	 } else {
  		std::cout << "Unable to open file"; 
  	 }
	return 0;
}