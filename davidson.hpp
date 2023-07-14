#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <chrono>
#include <iomanip>

//#define PROF true

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

template<class T>
void mgs(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &V) {
    const Eigen::Vector<T, Eigen::Dynamic> v = V.col(V.cols()-1);
    Eigen::Vector<T, Eigen::Dynamic> U = v;
    for (int i=0; i<V.cols()-1; i++) {
	U -= v.dot(V.col(i)) / V.col(i).dot(V.col(i)) * V.col(i);
    }
    V.col(V.cols()-1) = U/U.norm();
}

template<class T>
std::vector<Eigen::Vector<T, Eigen::Dynamic>> davidsonSolve(
		    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A,
		    std::vector<Eigen::Vector<T, Eigen::Dynamic>> b,
		    const double tol) {
    // init dimension of matrix
    const int n = A.rows();
    // init maximum number of iterations
    constexpr int iterlimit = 300;
    // number of roots
    const int nroots = b.size();
    int nconv; // number of converged roots


    // check for zero norm
    bool converged[nroots];
    for (int i=0; i<nroots; i++) converged[i] = false;
    for (int i=0; i<nroots; i++) {
	const double l2 = b[i].norm();
	if (l2 < 1e-18) converged[i] = true;
    }

    
    // DAVIDSON
    std::cout << "\n\n :: starting Davidson iterations...\n" << std::flush;
    std::cout << "\tnumber of rhs vectors: " << nroots << "\n";
    std::cout << "\tsize of vector space:  " << n << "\n";
    std::cout << "\tsize of rhs space:     " << b[0].rows() << "\n";
    std::cout << "\tconvergence criterion: " << tol << "\n\n" << std::flush;

    auto tstart = high_resolution_clock::now();

    // set up D Matrix
    T Darr[n];
    for (int i=0; i<n; i++) { Darr[i] = 1.0/A(i, i); }

    std::vector<Eigen::Vector<T, Eigen::Dynamic>> alld(nroots),
						allalpha(nroots),
						res(nroots);
    for (auto &x: res) x = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);

    // choose start vector
    Eigen::Vector<T, Eigen::Dynamic> v1 =
	Eigen::Vector<T, Eigen::Dynamic>::Unit(n, 0)
	+ 0.001*Eigen::Vector<T, Eigen::Dynamic>::Unit(n, 1);
    v1 = v1/v1.norm();

    // define basis
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> V;
    for (int i=0; i<nroots; i++) V.push_back(v1);

    // subspace matrix
    std::vector<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>
						    H(nroots), subL(nroots);
    Eigen::Vector<T, Eigen::Dynamic> delta(n);
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Vexpand(n, nroots);

    auto tsetup = high_resolution_clock::now();
    duration<double, std::milli> ms = tsetup - tstart;
#if PROF
    std::cout << "\tsetup done after " << ms.count() << " ms\n" << std::flush;
#endif

    // davidson iterations
    int iter = 1;
    int nexpand = 1;
    double maxres;
    for (int j=1; j<=iterlimit; j++) {
        std::cout << " ============= Davidson iteration " << iter
					<< " =============\n" << std::flush;
	nconv = 0;
	maxres = -10000000000;
	iter++;

	// iterate over all roots
	for (int root=0; root<nroots; root++) {
	    if (converged[root]) {nconv++; continue;}
	    // Project to subspace
	    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp
	    	    = V[root].col(j-1).adjoint() * A;
	    subL[root].conservativeResize(j, n);
	    subL[root].row(j-1) = tmp;
	    Eigen::Vector<T, Eigen::Dynamic> expansion 
	    	    = subL[root] * V[root].col(j-1);
	    H[root].conservativeResize(j, j);
	    H[root].col(j-1) = expansion;
	    H[root].row(j-1) = expansion.conjugate();
	    Eigen::Vector<T, Eigen::Dynamic> d = V[root].adjoint() * b[root];

	    // solve small system
	    Eigen::LLT<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>
								dec(H[root]);
	    const auto alpha = dec.solve(d);

	    // project out
	    const Eigen::Vector<T, Eigen::Dynamic> x = V[root] * alpha;
	    const Eigen::Vector<T, Eigen::Dynamic> r = b[root] - A*x;

	    maxres = std::max(maxres, r.norm());
            if (r.norm() <= tol) {nconv++; res[root] = x; converged[root]=true;}
            

	    for (int i=0; i<n; i++) { delta(i) = r(i) * Darr[i]; }

	    V[root].conservativeResize(V[root].rows(), V[root].cols()+1);
            V[root].col(V[root].cols()-1) = delta;
	    mgs(V[root]);
	}
        if (nconv == nroots) {std::cout << "CONVERGED!\n"; goto done;}
        std::cout << "\tmax. residual norm: " << maxres << "\n";
	std::cout << "\tnumber of converged roots: " << nconv << "\n";

	std::cout << "\n\n" << std::flush;
    } // end of davidson iterations

    done:
    std::cout << "\tmax. residual norm: " << maxres << "\n";
    std::cout << "\tnumber of converged roots: " << nconv << "\n";

    auto tend = high_resolution_clock::now();
    ms = tend - tstart;
    std::cout << "\ndone after " << ms.count() << " ms total\n" << std::flush;

    return res;
}
