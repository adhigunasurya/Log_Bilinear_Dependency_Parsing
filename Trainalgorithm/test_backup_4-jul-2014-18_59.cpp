#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <random>
#include <stdlib.h>
#include "graph.cpp"
#include <math.h>

#ifndef FEATURE_SIZE
#define FEATURE_SIZE 10
#endif

#ifndef NO_OF_POSTAGS
#define NO_OF_POSTAGS 50
#endif

#ifndef MAX_DISTANCE
#define MAX_DISTANCE 49
#endif

#ifndef NO_OF_DISTANCE
#define NO_OF_DISTANCE 24
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 100
#endif

#ifndef LEARNING_RATE
#define LEARNING_RATE 0.06
#endif

//no need for learning decay, use adagrad
#ifndef LEARNING_DECAY
#define LEARNING_DECAY 0.00001
#endif

#ifndef NO_OF_EPOCH
#define NO_OF_EPOCH 10
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

std::map<std::string, int> postag_mapping(std::string file_name) {
	std::string line;
	std::map<std::string, int> result;
	std::ifstream fine_pos_tags (file_name);
	int j = -1;
	if (!fine_pos_tags.is_open()) {
		std::cout << "Error: file fails to open" << std::endl;
	} else {
		while (getline(fine_pos_tags, line)) {
			j++;
			result[line] = j;
		}
		j++;
		result["START"] = j++;
		result["END"] = j++;
		result["ROOT"]  = j++;
		result["ROOT+"] = j++;
		result["ROOT-"] = j++;
	}
	fine_pos_tags.close();
	return result;
}

//ignore punctuation accuracy
double get_validation_accuracy(std::vector<Eigen::ArrayXXd> Ci, Eigen::ArrayXXd R, std::string validation_input_file, std::string output, std::map<std::string, int> postag_mapping, std::map<std::string, int> big_postag_mapping, Eigen::ArrayXd bv) {
	std::ifstream input_file(validation_input_file);
	std::ofstream output_file(output);
	if (!input_file.is_open() || !output_file.is_open()) {
		std::cout << "either input or output file fails to open" << std::endl;
	} else {
		std::string line;
		//get the input sentences from the validation_input_file
		int j = -1;
		while (getline(input_file, line)) {
			std::vector<std::string> sentence;
			sentence.push_back("START");
			std::vector<std::string> temp = split(line, ' ');
			int length = temp.size();
			for (std::vector<std::string>::iterator it = temp.begin(); it != temp.end(); it++) {
				sentence.push_back(*it);
			}
			sentence.push_back("END");
			//construct a graph with the given number of tokens
			Graph G(length); 
			//get the probabilities between each token
			for (int i = 0; i < sentence.size(); i++) { //i indicates all the possible modifiers
				for (int j = 0; j < sentence.size(); j++) { //j indicates all the possible heads
					if (sentence[i] != "START" && sentence[i] != "END" && sentence[j] != "START" && sentence[j] != "END" && i != j) {
						//calculate the probability
						int dist = abs(i-j);
						std::ostringstream s;
						s << dist;
						std::string distance = s.str();
						if (big_postag_mapping.find(sentence[i]+distance) != big_postag_mapping.end()) {
							double bias_term = bv(big_postag_mapping[sentence[i] + distance]);
						} else {
							std::cout << "postag not found" << std::endl;
						}
						//connect the graph
					}
				}
			}
			
			//get the MST
		}
		//calculate accuracy (ignoring punctuations)
		input_file.close();
	}
	return 0.0;
}

//returns a mapping between the pos tags, concatenated with all possible distance, with integers
//the difference with function postag_mapping is that this time the mapping is concatenated with all the possible distance (from 1-49), 
//while in the function postag_mapping, the pos tags are "vanilla" and not concatenated with the distance
std::map<std::string, int> get_postag_bias(Eigen::ArrayXd& bias, std::string distance_count_file, std::string postags_file) {
	std::string line;
	std::ifstream distance_count(distance_count_file);
	std::ifstream postags(postags_file);
	std::vector<std::string> temp;
	std::map<std::string, int> result_map;
	if (!distance_count.is_open() || !postags.is_open()) {
		std::cout << "Error: file containing distance counts fail to open" << std::endl;
	} else {
		int total = 0;
		//first, create the mapping between list of pos tags (appended with distance) and their counts, initialised to 0
		int j = 0;
		while (getline(postags, line)) {
					for (int i = -MAX_DISTANCE; i <= MAX_DISTANCE; i++) {
						if (i != 0) {
							if (i == -40 || i == -30 || i == -25 || i == -20 || i == -15 || i == -10 || i == -5 || i == -4 || i == -3 || i == -2 || i == -1 || i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6 || i == 11 || i == 16 || i == 21 || i == 26 || i == 31 || i == 41) {
								j++;
							}
							std::stringstream sstm;
							sstm << line << i;
							result_map[sstm.str()] = j;
						}
					}
		}
		/*
		std::ofstream test_output("newer_mapping.txt");
		for (std::map<std::string, int>::iterator it = result_map.begin(); it != result_map.end(); it++) {
			test_output << it->first << " " << "is mapped to " << it->second << std::endl;
		} */
		//calculate the total occurences of postags
		while (getline(distance_count, line)) {
			temp = split(line, ' ');
			total += std::stoi(temp[4]);
			temp.clear();
		}
		std::cout << total << std::endl;
		//std::cout << "the total counts for all distance is: " << total << std::endl;
		distance_count.close();
		//re-open the file and assign the probabilities (count of each pos tag divided by total) to the bias array
		distance_count.open(distance_count_file);
		bias = Eigen::ArrayXd::Zero(result_map.size());
		while (getline(distance_count, line)) {
			temp = split(line, ' ');
			if (result_map.find(temp[2]) == result_map.end()) {
				std::cout << "ERROR" << std::endl;
			} else {
				int idx = result_map[temp[2]];
				int pos_tag_count = std::stoi(temp[4]);
				bias(idx) = bias(idx) + ((double) pos_tag_count);
			}
			temp.clear();
		}
		bias /= ((double) total);
		distance_count.close();
		postags.close();
	}
	return result_map;
}

int main()
{
	//initialisation
	Eigen::ArrayXXd R = Eigen::ArrayXXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE); //matrix R, dimension: Nw(number of pos tags conjoined with distance) x Nf
	std::default_random_engine generator;
	std::normal_distribution<double> distr(0.0, 0.1);
	//initialise the matrix R with small, random, non-zero numbers
	for (int i = 0; i < R.rows(); i++) {
		for (int j = 0; j < R.cols(); j++) {
			double d = distr(generator);
			//avoid having 0's as the weight, try to get all the coefficients initialised with small numbers
			while (d == 0.0) {
				d = distr(generator);
			}
			R(i, j) = d;
		}
	}
		//initialise the matrix Q, of size feature_size x (no_of_postags * no_of_distance) 
	Eigen::ArrayXXd Q = Eigen::ArrayXXd::Zero(NO_OF_POSTAGS * NO_OF_DISTANCE, FEATURE_SIZE);
	for (int i = 0; i < Q.rows(); i++) {
		for (int j = 0; j < Q.cols(); j++) {
			double d = distr(generator);
			//avoid having 0's as the weight, try to get all the coefficients initialised with small numbers
			while (d == 0.0) {
				d = distr(generator);
			}
			Q(i, j) = d;
		}
	}
	//create an array of Ci (where each Ci is a matrix of Nf x Nf that governs the interaction between the ith feature and vn)
	std::vector<Eigen::ArrayXXd> Ci;
	for (int i = 0; i < 5; i++) {
		Eigen::ArrayXXd Ci_temp = Eigen::ArrayXXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
		std::default_random_engine generator;
		std::normal_distribution<double> distr(0.0, 0.1);
		for (int j = 0; j < Ci_temp.rows(); j++) {
			for (int k = 0; k < Ci_temp.cols(); k++) {
				double d = distr(generator);
				//avoid having 0's as the weight, try to get all the coefficients initialised with small numbers
				while (d == 0.0) {
					d = distr(generator);
				}
				Ci_temp(j, k) = d;
			}
		}
		Ci.push_back(Ci_temp);
	}
	//initialise vi. We use a context of 5 (the postag of the head itself, the postag before the modifier, the postag after the modifier, the postag before the head, and the postag after the head)
	//initialise vn. Both vi and vn are column vectors with all 0's and a single 1, which indicates the POS-tag we desire to get
	//FOR NOW, we ignore the term br^transpose * R^transpose*vn, but find out what it is
	//initialise the bias term bv
	Eigen::ArrayXd bv;
	//print the array of biases into an external file
	std::ofstream final_postag_bias("final_postag_bias.txt"); 
	if (!final_postag_bias.is_open()) {
		std::cout << "fail to create final postag bias file" << std::endl;
	} else {
		std::map<std::string, int> big_pos_tag_map = get_postag_bias(bv, "distance_count.txt", "fine_pos_tags.txt");
		for (std::map<std::string, int>::iterator it = big_pos_tag_map.begin(); it != big_pos_tag_map.end(); it++) {
			final_postag_bias << "Pos tag " << it->first << " has bias value " << bv(it->second) << std::endl;
		}
	} 
	//std::cout << "one example: " << bv(1) << std::endl;
	//main bit
	std::map<std::string, int> mapping = postag_mapping("fine_pos_tags.txt");
	std::ifstream pairs ("real_pairs_small.txt");
	//initialise the gradient matrix
	std::vector<Eigen::MatrixXd> dCi;
	for (int i = 0; i < 5; i++) {
		Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
		dCi.push_back(temp);
	}
	Eigen::MatrixXd dR = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE);
	int no_of_iter = 0;
	std::ofstream parameters ("parameters_each_epoch_small_trial_6.txt");
	if (!parameters.is_open()) {
		std::cout << "Error: output file not opening" << std::endl;
	}
	Eigen::MatrixXd tempC = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
	Eigen::MatrixXd tempR_D = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE);
	Eigen::MatrixXd tempR_M = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE);
	Eigen::MatrixXd temp_matrix = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE);
	Eigen::MatrixXd temp_matrix_2 = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE);
	double normalising_constant = 0.0;
	//calculate the normalising constant for the first time
	for (int i = 0; i < NO_OF_EPOCH; i++) {
		if (!pairs.is_open()) {
			std::cout << "Fail to open input file" << std::endl;
		} else {
			std::vector<std::string> temp;
			std::string line;
			while (getline(pairs, line)) {
				if (line.size() > 1) {
					no_of_iter++;
					temp = split(line, ' ');
					//initialise an array with the elements as the indices of: postag of the modifier, postag of the head, postag of before modifier, postag of after modifier, postag of before head, postag of after head
					int arr_temp[6] = {mapping[temp[0]], mapping[temp[1]], mapping[temp[2]], mapping[temp[3]], mapping[temp[4]], mapping[temp[5]]};
					std::cout << no_of_iter << std::endl;
					//calculate all the gradient changes for each of the Ci, from dCi[0] to dCi[4]
					//this operation is already efficient
					for (int j = 0; j < 5; j++) {
						for (int k = 0; k < R.rows(); k++) {
							tempC += (R.row(arr_temp[j+1]).transpose().matrix() * R.row(k).matrix());
						}
						//normalise
						tempC /= NO_OF_POSTAGS; //wrong
						dCi[j] += ((R.row(arr_temp[j+1]).transpose().matrix() * R.row(arr_temp[0]).matrix()) - tempC);
						//reset the matrix tempC to prepare to calculate the next dCi[j]
						tempC.setZero();
					}
					//calculate dR 
					//calculate tempD
					for (int j = 0; j < 5; j++) {
						//source of inefficiency: temp_matrix & temp_matrix_2 is a matrix of all zeros except one row!
						//no need to multiply everything!
						temp_matrix.row(arr_temp[0]) = R.row(arr_temp[j+1]);
						temp_matrix_2.row(arr_temp[j+1]) = R.row(arr_temp[0]);
						tempR_D += ((temp_matrix * Ci[j].matrix()) + (temp_matrix_2 *  Ci[j].matrix().transpose()));
						temp_matrix.setZero();
						temp_matrix_2.setZero();
					}
					//calculate tempM
					for (int j = 0; j < NO_OF_POSTAGS; j++) {
						for (int k = 0; k < 5; k++) {
							temp_matrix.row(j) = R.row(arr_temp[k+1]);
							temp_matrix_2.row(arr_temp[k+1]) = R.row(j);
							tempR_M += ((temp_matrix * Ci[k].matrix()) + (temp_matrix_2 *  Ci[k].matrix().transpose()));
							temp_matrix.setZero();
							temp_matrix_2.setZero();
						}
					}
					//normalise
					tempR_M /= NO_OF_POSTAGS; //wrong
					//calculate dR
					dR += (tempR_D - tempR_M);
					//reset tempD and tempM
					tempR_D.setZero();
					tempR_M.setZero();

					//update the gradient if no_of_iter mod mini batch size == 0, reset the gradient matrix to 0
					if (no_of_iter % BATCH_SIZE == 0) {
						//update normalising constant (unimplemented)
						double curr_learning_rate = LEARNING_RATE / (1 + no_of_iter * LEARNING_DECAY);
						for (int j = 0; j < 5; j++) {
							Ci[j] = Ci[j].matrix() + curr_learning_rate * dCi[j];
							//reset the matrix dCi[j]
							dCi[j].setZero();
						}
						//Ci[4] = Ci[4].matrix() + curr_learning_rate * dC5;
						R = R.matrix() + curr_learning_rate * dR;

						//reset the matrices dR
						dR.setZero();
					}
				}
				temp.clear();
			}
			pairs.close();
		}
		parameters << "Epoch number: " << i+1 << std::endl;
		for (int j = 0; j < 5; j++) {
			parameters << "Matrix C" << j+1 << std::endl; 
			parameters << Ci[j] << std::endl;
			parameters << std::endl;
			parameters << "End of matrix C" << j+1 << std::endl;
		}
		parameters << "Matrix R" << std::endl;
		parameters << R << std::endl;
		parameters << "End of matrix R" << std::endl;
		pairs.open("real_pairs_small.txt");
	}
	pairs.close();
	parameters.close();
	return 0;
}