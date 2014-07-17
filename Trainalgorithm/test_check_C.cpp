#include <Eigen/Dense>
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

#ifndef FEATURE_SIZE
#define FEATURE_SIZE 10
#endif

#ifndef NO_OF_POSTAGS
#define NO_OF_POSTAGS 47
#endif

#ifndef NO_OF_POSTAGS_WITHOUT_DUMMY
#define NO_OF_POSTAGS_WITHOUT_DUMMY 45
#endif

#ifndef MAX_DISTANCE
#define MAX_DISTANCE 49
#endif

#ifndef NO_OF_DISTANCE
#define NO_OF_DISTANCE 28
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 50
#endif

#ifndef LEARNING_RATE
#define LEARNING_RATE 0.0018 //LEARNING RATE USED FOR PREVIOUS EXPERIMENTS: 0.06
#endif


//no need for learning decay, use adagrad
#ifndef LEARNING_DECAY
#define LEARNING_DECAY 0.0000033
#endif 

#ifndef NO_OF_EPOCH
#define NO_OF_EPOCH 15
#endif

#ifndef EPSILON
#define EPSILON 0.0001
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

double get_maximum_bias(std::map<int, double> bias) {
	double curr_max = -999.0;
	for (std::map<int, double>::iterator it = bias.begin(); it != bias.end(); it++) {
		if (it->second > curr_max) {
			curr_max = it->second;
		}
	}
	return curr_max;
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
	}
	fine_pos_tags.close();
	return result;
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
		int j = -1;
		while (getline(postags, line)) {
					for (int i = -MAX_DISTANCE; i <= MAX_DISTANCE; i++) {
						if (i != 0) {
							if (i == -49 || i == -45 || i == -40 || i == -35 || i == -30 || i == -25 || i == -20 || i == -15 || i == -10 || i == -5 || i == -4 || i == -3 || i == -2 || i == -1 || i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6 || i == 11 || i == 16 || i == 21 || i == 26 || i == 31 || i==36 || i == 41 || i == 46) {
								j++;
							}
							std::stringstream sstm;
							sstm << line << i;
							result_map[sstm.str()] = j;
						}
					}
		}
		if (j == (NO_OF_POSTAGS_WITHOUT_DUMMY * NO_OF_DISTANCE - 1)) {
			std::cout << "The size of pos tag mapping is correct" << std::endl;
		} else {
			std::cout << "The size of pos tag mapping is INCORRECT" << std::endl;
			std::cout << "j is " << j << std::endl;
		}
		/*
		std::ofstream test_output("newer_mapping.txt");
		for (std::map<std::string, int>::iterator it = result_map.begin(); it != result_map.end(); it++) {
			test_output << it->first << " " << "is mapped to " << it->second << std::endl;
		} */
		//calculate the total occurences of postags
		while (getline(distance_count, line)) {
			temp = split(line, ' ');
			total += std::stoi(temp[3]);
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
	Eigen::ArrayXXd R = Eigen::ArrayXXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS); //matrix R, dimension: Nf x Nw
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
	Eigen::ArrayXXd Q = Eigen::ArrayXXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
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
	} 
	/*std::map<std::string, int> big_pos_tag_map = get_postag_bias(bv, "distance_count.txt", "fine_pos_tags.txt");
	std::cout << "size of big pos tag map = " << big_pos_tag_map.size() << std::endl;
	for (std::map<std::string, int>::iterator it = big_pos_tag_map.begin(); it != big_pos_tag_map.end(); it++) {
		final_postag_bias << "Pos tag " << it->first << " has bias value " << bv(it->second) << std::endl;
	} */

	/* create a pair between each "target pos tag" with a bias value (tested and correct)*/
	//Eigen::ArrayXd word_bias_pair = Eigen::ArrayXd::Zero(NO_OF_POSTAGS_WITHOUT_DUMMY * NO_OF_DISTANCE);
	int i_3 = -1;
	std::ifstream fine_pos_tags("fine_pos_tags.txt");
	std::string temp_string;
	/*
	while (getline(fine_pos_tags, temp_string)) {
		for (int i = -MAX_DISTANCE; i <= MAX_DISTANCE; i++) {
			if (i != 0) {
				if (i == -49 || i == -45 || i == -40 || i == -35 || i == -30 || i == -25 || i == -20 || i == -15 || i == -10 || i == -5 || i == -4 || i == -3 || i == -2 || i == -1 || i == 1 || i == 2 || i == 3 || i == 4 || i == 5 || i == 6 || i == 11 || i == 16 || i == 21 || i == 26 || i == 31 || i == 36|| i == 41 || i == 46) {
					i_3++;
					std::stringstream sstm;
					sstm << temp_string << i;
					if (big_pos_tag_map.find(sstm.str()) == big_pos_tag_map.end()) {
						std::cout << "ERROR" << std::endl;
					} else {
						int pos_tag_idx = big_pos_tag_map[sstm.str()];
						//word_bias_pair(i_3) = bv(pos_tag_idx);
					}
				}
			}
		}
	} */
	/* end */
	//main training part
	std::map<std::string, int> mapping = postag_mapping("fine_pos_tags.txt");
	fine_pos_tags.close();
	//print out both the small and big mapping
	//std::ofstream big_pos_tag_mapping("big_pos_tag_map.txt");
	std::ofstream small_pos_tag_mapping("small_pos_tag_map.txt");
	/*
	for (std::map<std::string, int>::iterator it = big_pos_tag_map.begin(); it != big_pos_tag_map.end(); it++) {
		big_pos_tag_mapping << it->first << " " << it->second << std::endl;
	} */
	for (std::map<std::string, int>::iterator it = mapping.begin(); it != mapping.end(); it++) {
		small_pos_tag_mapping << it->first << " " << it->second << std::endl;
	}
	small_pos_tag_mapping.close();

	//create the pairs
	std::ifstream pairs ("real_pairs_medium.txt");
	//initialise the gradient matrix
	std::vector<Eigen::MatrixXd> dCi;
	for (int i = 0; i < 5; i++) {
		Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
		dCi.push_back(temp);
	}
	Eigen::MatrixXd dR = Eigen::MatrixXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
	Eigen::MatrixXd dQ = Eigen::MatrixXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
	std::ofstream parameters ("parameters_each_epoch_small_trial_1_check_C.txt");
	if (!parameters.is_open()) {
		std::cout << "Error: output file not opening" << std::endl;
	}
	
	//initialise all the temporary variables outside the loop to save time on memory allocation on the heap (expensive)
	//all these temporary variables are reset to zero using .setZero() method at the end of each loop
	Eigen::MatrixXd tempC = Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE);
	Eigen::MatrixXd temp_R_Model = Eigen::MatrixXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
	Eigen::MatrixXd temp_R_Data = Eigen::MatrixXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
	Eigen::MatrixXd temp_Q = Eigen::MatrixXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
	Eigen::MatrixXd temp_Q_Model = Eigen::MatrixXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
	Eigen::MatrixXd temp_Q_Data = Eigen::MatrixXd::Zero(FEATURE_SIZE, NO_OF_POSTAGS);
	Eigen::VectorXd temp_Qw = Eigen::VectorXd::Zero(FEATURE_SIZE); 
	//Eigen::VectorXd dbias = Eigen::VectorXd::Zero(NO_OF_POSTAGS_WITHOUT_DUMMY * NO_OF_DISTANCE);
	//Eigen::VectorXd temp_dbias = Eigen::VectorXd::Zero(NO_OF_POSTAGS_WITHOUT_DUMMY * NO_OF_DISTANCE);
	Eigen::VectorXd tempCr = Eigen::VectorXd::Zero(FEATURE_SIZE);
	Eigen::ArrayXd temp_e = Eigen::ArrayXd::Zero(NO_OF_POSTAGS);
	/*
	Eigen::MatrixXd temp_matrix = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE);
	Eigen::MatrixXd temp_matrix_2 = Eigen::MatrixXd::Zero(NO_OF_POSTAGS, FEATURE_SIZE); */
	//calculate the normalising constant for the first time
	std::ofstream gradient_check_C("gradient_check_C_big.txt");
	std::vector<Eigen::MatrixXd> curr_dC; //DELETE THIS LATER
	for (int i = 0; i < 5; i++) {
		curr_dC.push_back(Eigen::MatrixXd::Zero(FEATURE_SIZE, FEATURE_SIZE));
	}
	int no_of_iter = 0;
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
					if (no_of_iter % 10000 == 0) {
						std::cout << no_of_iter << std::endl;
					}
					//calculate the normalising constant
					double normalising_constant = 0.0;
					for (int iter2 = 0; iter2 < 5; iter2++) {
							tempCr += (Ci[iter2].matrix() * R.col(arr_temp[iter2+1]).matrix()); //prediction
					}
					for (int iter = 0; iter < Q.cols(); iter++) {
						double temp = exp((tempCr.transpose() * Q.col(iter).matrix()));//dot product
						temp_e(iter) = temp; 
						normalising_constant += temp; 
					} 					
					//calculate all the gradient changes for each of the Ci, from dCi[0] to dCi[4]
					//this operation is already efficient
					for (int j = 0; j < 5; j++) {
						for (int k = 0; k < Q.cols(); k++) {
							tempC += (temp_e(k) * (Q.col(k).matrix() * R.col(arr_temp[j+1]).transpose().matrix()));
						}
						//normalise
						tempC /= normalising_constant;
						if (no_of_iter  % 1000 == 0) {
							curr_dC[j] = (Q.col(arr_temp[0]).matrix() * R.col(arr_temp[j+1]).transpose().matrix() - tempC);
						}
						dCi[j] += (Q.col(arr_temp[0]).matrix() * R.col(arr_temp[j+1]).transpose().matrix() - tempC);
						//reset the matrix tempC to prepare to calculate the next dCi[j]
						tempC.setZero();
					}
					
					//calculate dQ VERIFY
					for (int j = 0; j < 5; j++) {
						temp_Qw += (Ci[j].matrix() * R.col(arr_temp[j+1]).matrix());
					}
					for (int j = 0; j < Q.cols(); j++) {
						temp_Q_Model.col(j) += (temp_e(j) * temp_Qw);
					}
					temp_Q_Model /= normalising_constant;
					temp_Q_Data.col(arr_temp[0]) = temp_Qw;
					dQ += (temp_Q_Data - temp_Q_Model);
					temp_Qw.setZero();
					temp_Q_Model.setZero();
					temp_Q_Data.setZero();

					//calculate dR
					//calculate dR with respect to data
					for (int j = 0; j < 5; j++) {
						temp_R_Data.col(arr_temp[j+1]) += Ci[j].transpose().matrix() * Q.col(arr_temp[0]).matrix();
					}
					//calculate dR with respect to model
					for (int k = 0; k < Q.cols(); k++) {
						for (int j = 0; j < 5; j++) {
							temp_R_Model.col(arr_temp[j+1]) += (temp_e(k) * (Ci[j].transpose().matrix() * Q.col(k).matrix()));
						}
					}
					temp_R_Model /= normalising_constant;
					dR += (temp_R_Data - temp_R_Model);
					temp_R_Data.setZero();
					temp_R_Model.setZero();

					//calculate dbias
					//temp_dbias(arr_temp[0]) = 1.0;
					//dbias += (temp_dbias - temp_e.matrix() / normalising_constant);
					//temp_dbias(arr_temp[0]) = 0.0;

					tempCr.setZero();
					temp_e.setZero();

					if (no_of_iter % 1000 == 0) {
							Eigen::VectorXd temp = Eigen::VectorXd::Zero(FEATURE_SIZE);
							gradient_check_C <<  "current data point = " << no_of_iter <<std::endl;
							for (int iter = 0; iter <Ci.size(); iter++) {
								for (int j = 0; j < Ci[iter].rows(); j++) {
									for (int k = 0; k < Ci[iter].cols(); k++) {
										Eigen::ArrayXXd C_copy = Ci[iter];
										double result_1 = 0.0;
										double result_2 = 0.0;
										double temp_double = 0.0;
										gradient_check_C << curr_dC[iter](j, k) << " ";
										C_copy(j, k) += EPSILON;
										for (int iter2 = 0; iter2 < 5; iter2++) {
											if (iter2 == iter) {
												temp += (C_copy.matrix() * R.col(arr_temp[iter2+1]).matrix());
											} else {
												temp += (Ci[iter2].matrix() * R.col(arr_temp[iter2+1]).matrix());
											}
										}
										result_1 += temp.transpose() * Q.col(arr_temp[0]).matrix();
										for (int l = 0; l < Q.cols(); l++) {
											temp_double += exp(temp.transpose() * Q.col(l).matrix());
										}
										result_1 -= log(temp_double);
										temp_double = 0.0;
										//calculate result_2
										C_copy(j, k) -= (2 * EPSILON);
										temp.setZero();
										for (int iter2 = 0; iter2 < 5; iter2++) {
											if (iter2 == iter) {
												temp += (C_copy.matrix() * R.col(arr_temp[iter2+1]).matrix());
											} else {
												temp += (Ci[iter2].matrix() * R.col(arr_temp[iter2 + 1]).matrix());
											}	
										}
										result_2 += temp.transpose() * Q.col(arr_temp[0]).matrix();
										for (int l = 0; l < Q.cols(); l++) {
											temp_double += exp(temp.transpose() * Q.col(l).matrix());
										}
										result_2 -= log(temp_double);
										temp.setZero();
										gradient_check_C << (result_1 - result_2) / (2 * EPSILON) << std::endl;

									}
							}
							curr_dC[iter].setZero();
						}
				}

					//update the gradient if no_of_iter mod mini batch size == 0, reset the gradient matrix to 0

					if (no_of_iter % BATCH_SIZE == 0) {
						double curr_learning_rate = LEARNING_RATE / (1 + no_of_iter * LEARNING_DECAY);
						for (int j = 0; j < 5; j++) {
							Ci[j] = Ci[j].matrix() + curr_learning_rate * dCi[j];
							//reset the matrix dCi[j]
							dCi[j].setZero();
						}
						//word_bias_pair = word_bias_pair + curr_learning_rate * dbias.array();
						R = R.matrix() + curr_learning_rate * dR;
						Q = Q.matrix() + curr_learning_rate * dQ;
						dR.setZero();
						dQ.setZero();
						//dbias.setZero();
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
		parameters << "Matrix Q" << std::endl;
		parameters << Q << std::endl;
		parameters << "End of matrix Q" << std::endl;
		//parameters << "Bias vector" << std::endl;
		//parameters << word_bias_pair << std::endl;
		//parameters << "End of bias vector" << std::endl;
		pairs.open("real_pairs_medium.txt");
	} 
	pairs.close();
	parameters.close();
	return 0;
}