#ifndef RUN_GSA_H
#define RUN_GSA_H

/* *****************************************************************************
Copyright (c) 2016-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

/**
 *  @author  Sang-ri Yi
 *  @date    8/2021
 *  @section DESCRIPTION
 *  Runs global sensitivity analysis. See: Hu, Z. and Mahadevan, S. (2019). Probability models for data-driven global sensitivity analysis. Reliability Engineering & System Safety, 187, 40-57.
 */


#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <set>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric> // std::iota
#include <armadillo>
#include "jsonInput.h"

#include "writeErrors.h"
extern writeErrors theErrorFile; // Error log
using namespace arma;
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
using std::vector;

class runGSA
{
public:
	runGSA();
	runGSA(vector<vector<double>> xval,
			vector<vector<double>> gval,
			vector<vector<int>> combs_tmp,
			int Kos,
			int procno,
			int nprocs);
	~runGSA();
	void writeOutputs(jsonInput inp, double dur, int procno);
	void writeTabOutputs(jsonInput inp, int procno);


	//vector<double> Si;

	vector<vector<double>> xval;
	vector<vector<double>> gval;
	vector<vector<int>> combs_tmp;
	char Opt;
	//int Kos;
	vector<vector<double>> Simat;
	vector<vector<double>> Stmat;

private:
	double mvnPdf(mat x, mat mu, mat cov);
	double calMean(vector<double> x);
	double calVar(vector<double> x);
	const double PI = 3.1415926535897932384626433;
	vector<double> doGSA(vector<double> gval, int Kos, char Opt);
	int nrv;
	int ncombs;
	int nmc;
};

#endif //RUN_GSA_H
