#include <xerus.h>
#include <iostream>
#include <string>
#include <fstream>
#include <boost/timer/timer.hpp>

using namespace xerus;

const size_t MAX_NUM_PER_SITE = 6;

Tensor create_M() {
	Tensor M = -1*Tensor::identity({MAX_NUM_PER_SITE, MAX_NUM_PER_SITE});
	for (size_t i = 0; i < MAX_NUM_PER_SITE-1; ++i) {
		M[{i+1, i}] = 1.0;
	}
	return M;
}

Tensor create_L() {
	Tensor L({MAX_NUM_PER_SITE, MAX_NUM_PER_SITE});
	L.modify_diagonal_entries([](value_t& _x, const size_t _i) { 
		_x = double(_i)/double(_i+5); 
	});
	return L;
}

Tensor create_S() {
	Tensor S({MAX_NUM_PER_SITE, MAX_NUM_PER_SITE});
	
	// Set diagonal
	for (size_t i = 0; i < MAX_NUM_PER_SITE; ++i) {
		S[{i, i}] = -double(i);
	}
	
	// Set offdiagonal
	for (size_t i = 0; i < MAX_NUM_PER_SITE-1; ++i) {
		S[{i, i+1}] = double(i+1);
	}
	return 0.07*S;
}



TTOperator create_operator(const size_t _degree) { 
	const Index i, j, k, l;
	
	// Create matrices
	const Tensor M = create_M();
	const Tensor S = create_S();
	const Tensor L = create_L();
	const Tensor Sstar = 0.7*M+S;
	const Tensor I = Tensor::identity({MAX_NUM_PER_SITE, MAX_NUM_PER_SITE});
	
	// Create empty TTOperator
	TTOperator A(2*_degree);
	
	Tensor comp;
	
	// Create first component
	comp(i, j, k, l) = 
		Sstar(j, k)*Tensor::dirac({1, 3}, 0)(i, l) 
		+   L(j, k)*Tensor::dirac({1, 3}, 1)(i, l) 
		+   I(j, k)*Tensor::dirac({1, 3}, 2)(i, l);
    
	A.set_component(0, comp); 
	
	// Create middle components
	comp(i, j, k, l) =
		  I(j, k) * Tensor::dirac({3, 3}, {0, 0})(i, l)
		+ M(j, k) * Tensor::dirac({3, 3}, {1, 0})(i, l)
		+ S(j, k) * Tensor::dirac({3, 3}, {2, 0})(i, l)
		+ L(j, k) * Tensor::dirac({3, 3}, {2, 1})(i, l)
		+ I(j, k) * Tensor::dirac({3, 3}, {2, 2})(i, l);
	
	for(size_t c = 1; c+1 < _degree; ++c) {
		A.set_component(c, comp);
	}
	
	// Create last component
	comp(i, j, k, l) = 
		  I(j, k)*Tensor::dirac({3, 1}, 0)(i, l) 
		+ M(j, k)*Tensor::dirac({3, 1}, 1)(i, l) 
		+ S(j, k)*Tensor::dirac({3, 1}, 2)(i, l);
    
	A.set_component(_degree-1, comp);
	
	return A;
}

/// @brief calculates the one-norm of a TTTensor, assuming that all entries in it are positive
double one_norm(const TTTensor &_x) {
	Index j;
	return double(_x(j&0) * TTTensor::ones(_x.dimensions)(j&0));
}

std::vector<Tensor> implicit_euler(const Tensor& _A, Tensor _x, 
			const double _stepSize, const size_t _n) 
{ 
	const Tensor op = Tensor::identity(_A.dimensions)-_stepSize*_A;
	
	Index j,k;
// 	auto ourALS = ALS_SPD;
// 	ourALS.convergenceEpsilon = 0;
// 	ourALS.numHalfSweeps = 2;
	
	std::cout << "PASSED!" << std::endl;
	
	std::vector<Tensor> results;
	Tensor nextX = _x;
	results.push_back(_x);
    
	for(size_t i = 0; i < _n; ++i) {
		//ourALS(op, nextX, _x);
		nextX(k&0) = _x(j&0)/op(j/2,k/2);
        
		// Normalize
		double norm = one_norm(nextX);
		nextX /= norm;
		
		XERUS_LOG(iter, "Done itr " << i 
				<< " residual: " << frob_norm(op(j/2,k/2)*nextX(k&0) - _x(j&0)) 
				<< " norm: " << norm);
		
		_x = nextX;
		results.push_back(_x);
	}
	
	return results;
}

double get_mean_concentration(const Tensor& _res, const size_t _i) { 
	const Index k,l;
	TensorNetwork result(_res);
		const Tensor weights({MAX_NUM_PER_SITE}, [](const size_t _k){ 
		return double(_k); 
	});
	const Tensor ones = Tensor::ones({MAX_NUM_PER_SITE});
	
	for (size_t j = 0; j < _res.degree(); ++j) {
		if (j == _i) {
			result(l&0) = result(k, l&1) * weights(k);
		} else {
			result(l&0) = result(k, l&1) * ones(k);
		}
	}
	// at this point the degree of 'result' is 0, so there is only one entry
	return result[{}]; 
}

void print_mean_concentrations_to_file(const std::vector<Tensor> &_result) {
	std::fstream out("mean_noTT.dat", std::fstream::out);
	for (const auto& res : _result) {
		for (size_t k = 0; k < res.degree(); ++k) {
			out << get_mean_concentration(res, k) << ' ';
		}
		out << std::endl;
	}
}

int main() {
	const size_t numProteins = 5;
	const size_t numSteps = 200;
	const double stepSize = 0.1;
	const size_t rankX = 3;
	
	auto start = Tensor::dirac(
			std::vector<size_t>(numProteins, MAX_NUM_PER_SITE), 
			0
	);
	start.use_dense_representation();
	start += 1e-14 * Tensor::random(
			start.dimensions
//			std::vector<size_t>(start.degree()-1, rankX-1)
	);
	const auto A = create_operator(numProteins);

	Tensor A2(A);
	using xerus::misc::operator<<;
	std::cout << "degree: " << A2.degree() << " dim: " << A2.dimensions << std::endl;
// 	
	boost::timer::cpu_timer timer;
	const auto results = implicit_euler(A2, start, stepSize, numSteps);
	boost::timer::cpu_times times = timer.elapsed();
	
	std::cout << "WALL TIME: " << std::setprecision(4) << times.wall / 1e9 << "s" << std::endl;
// 	
	print_mean_concentrations_to_file(results);
}
