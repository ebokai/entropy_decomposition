#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <bit>

using namespace std;

const int n = 20;
const int n_states = pow(2,n);
const string fname = "../data/bfpt_data_binary_n20_CSN_OPN.dat"; 

vector<uint64_t> read_data_file(int *N);
vector<double> get_pdata(vector<uint64_t> data, int N);
vector<double> get_entropy(vector<double> pdata, vector<double> log_pdata);
vector<double> get_log_pdata(vector<double> pdata);

// ==============================================

int main() {
	int N = 0;

	auto data = read_data_file(&N);
	auto pdata = get_pdata(data, N);
	auto log_pdata = get_log_pdata(pdata);
	auto entropies = get_entropy(pdata, log_pdata);

	return 0;
}

// ==============================================

vector<uint64_t> read_data_file(int *N) {
	string line, subline;

	int i = 0;
	uint64_t nline; 


	// get number of datapoints 
	ifstream datafile(fname);

	while (getline (datafile, line)) {
		(*N)++;
	}

	const int size = (*N);

	cout << *N << endl;

	// specify vector of length N and assign data
	vector<uint64_t> data(size);

	// reset file cursor to beginning
	datafile.clear();
	datafile.seekg(0, datafile.beg);

	while (getline (datafile, line)) {

		subline = line.substr(0,n);
		nline = bitset<n>(subline).to_ulong();
		data[i] = nline;
		i++;

	}

	return data;

}

// ==============================================

vector<double> get_pdata(vector<uint64_t> data, int N) {

	uint64_t state = 0;

	vector<double> pdata(n_states);

	for (int i = 0; i < N; i++) {

		state = data[i];
		pdata[state] += 1/ static_cast <float> (N);

	}

	return pdata;
}

vector<double> get_log_pdata(vector<double> pdata) {

	vector<double> log_pdata(n_states);
	double pd; 

	for (int i = 0; i < n_states; i++){

		pd = pdata[i];


		if (pd > 0) {

			log_pdata[i] = log(pd);
		} else {

			log_pdata[i] = 0;
		}
	}

	return log_pdata;
}

// ==============================================

vector<double> get_entropy(vector<double> pdata, vector<double> log_pdata) {

	vector<double> entropies(n_states);

	int op_bits;
	int eval_op;

	double phi = 0;
	double g = 0;
	
	double pd;
	double lpd;

	for (int op = 0; op < n_states; op++) {

		phi = 0;
		g = 0;

		if (op % 512 == 0){
			cout << op << endl;
		}

		for (int state = 0; state < n_states; state++) {

			pd = pdata[state];

			if (pd > 0) {
				lpd = log_pdata[state];
				op_bits = bitset<n>(op & state).count();
				// op_bits = popcount(op & state);
				op_bits %= 2;
				eval_op = 1 - 2 * op_bits;
				phi += pd * eval_op;
				g += lpd * eval_op / static_cast <float> (n_states);
			} 
		}

		entropies[op] = -g * phi;
		// cout << op << " " << bitset<n>(op) << " " << entropies[op] << endl;
	}

	return entropies;
}