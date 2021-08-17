
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

#include "runGSA.h"
#include <iterator>

runGSA::runGSA() {}

runGSA::runGSA(vector<vector<double>> xval,
	vector<vector<double>> gmat,
	vector<vector<int>> combs_tmp,
	int Kos,
	int procno,
	int nprocs)
{
	this->xval = xval;
	this->gval = gmat;
	this->combs_tmp = combs_tmp;
	nmc = xval.size();

	nrv = xval[0].size();
	ncombs = combs_tmp.size();
	int nqoi = gmat[0].size();
	int Kos_base = std::min(Kos, int(nmc / 5));
	Kos = Kos_base;
	//
	// MPI
	//
	int chunkSize = std::ceil(double(nqoi) / double(nprocs));
	//int lastChunk = inp.nmc - chunkSize * (nproc-1);
	double* SmAll = (double*)malloc(ncombs * chunkSize * nprocs * sizeof(double));
	double* SmTmp = (double*)malloc(ncombs * chunkSize * sizeof(double));
	double* StAll = (double*)malloc(ncombs * chunkSize * nprocs * sizeof(double));
	double* StTmp = (double*)malloc(ncombs * chunkSize * sizeof(double));
	// for each QoI
	for (int nq = 0; nq < chunkSize ; nq++) {
		int id = chunkSize * procno + nq;
		if (id >= nqoi) { // dummy
			for (int i = 0; i < ncombs; i++) {
				StTmp[nq * ncombs + i] = 0.;
				SmTmp[nq * ncombs + i] = 0.;
			}
			continue;
		}

		vector<double> gvec;
		double sqDiff = 0;
		gvec.reserve(nmc);
		for (int i = 0; i < nmc; i++) {
			gvec.push_back(gmat[i][id]);
		}

		// check if the variance is zero
		double mean = 0;
		for (int i = 0; i < nmc; i++)
			mean += gvec[i];

		mean = mean / double(nmc);
		for (int i = 0; i < nmc; i++)
			sqDiff += (gmat[i][id] - mean) * (gmat[i][id] - mean);

		//double var = sqDiff / nmc;
		if (sqDiff < 1.e-10) {
			//theErrorFile << "Error running FEM: the variance of output is zero. Output value is " << mean;
			//theErrorFile.close();
			//exit(1);
			//vector<double> zeros(ncombs, 0.0);
			//Simat.push_back(zeros);
			//Stmat.push_back(zeros);
			//continue;
			for (int i = 0; i < ncombs; i++) {
				StTmp[nq * ncombs + i] = 0.;
				SmTmp[nq * ncombs + i] = 0.;
			}
			continue;
		};

		vector<double> Sij, Stj;


		// repeat GSA with different Kos untill success
		double failIdx = -100, i = 1;
		while ((failIdx == -100) || (Kos / i < 0.5))
		{
			Sij = doGSA(gvec, ceil(Kos / i), 'M');
			failIdx = Sij[0];
			i *= 2;
		}
		Kos = Kos_base;

		failIdx = -100, i = 1;
		while ((failIdx == -100) || (Kos / i < 0.5))
		{
			Stj = doGSA(gvec, ceil(Kos / i), 'T');
			failIdx = Stj[0];
			i *= 2;
		}
		Kos = Kos_base;


		for (int i = 0; i < ncombs; i++) {
			StTmp[nq * ncombs + i] = Stj[i];
			SmTmp[nq * ncombs + i] = Sij[i];
		}
		//Simat.push_back(Stj);
		//Stmat.push_back(Sij);
	}
	MPI_Allgather(StTmp, ncombs * chunkSize, MPI_DOUBLE, StAll, ncombs * chunkSize, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(SmTmp, ncombs * chunkSize, MPI_DOUBLE, SmAll, ncombs * chunkSize, MPI_DOUBLE, MPI_COMM_WORLD);
	for (int i = 0; i < nqoi; i++) {
		vector<double> StVectmp(ncombs,0), SmVectmp(ncombs, 0);
		for (int j = 0; j < ncombs; j++) {
			StVectmp[j] = StAll[i * ncombs + j];
			SmVectmp[j] = SmAll[i * ncombs + j];
		}
		Stmat.push_back(StVectmp);
		Simat.push_back(SmVectmp);
	}


}

vector<double> runGSA::doGSA(vector<double> gval,int Kos,char Opt)
{
	vector<vector<int>> combs;

	if (Opt == 'T')
	{
		vector<int> allSet(nrv);
		std::iota(allSet.begin(), allSet.end(), 0);

		for (auto comb : combs_tmp)
		{
			vector<int> cnc;
			//std::set_difference(allSet.begin(), allSet.end(), comb.begin(), comb.end(), cnc.begin());
			std::set_difference(allSet.begin(), allSet.end(), comb.begin(), comb.end(), std::inserter(cnc, cnc.begin()));
			combs.push_back(cnc);
		}
	}
	else
	{
		combs = combs_tmp;
	}

	double V = calVar(gval);
	double Vi;
	vector<double> Si;
	Si.reserve(ncombs);

	for (int nc = 0; nc < ncombs; nc++)
	{

		const int endm = combs[nc].size(); // (nx+ng)-1
		const int endx = endm - 1;			// (nx)-1
		if (endm == 0)
		{
			if (Opt == 'T')
			{
				Si.push_back(1.); // total
			}
			else
			{
				Si.push_back(0.);   // main
			}
			printf("GSA i=%i, Si=%.2f, K=%i \n", nc + 1, Si[nc], Kos);
			continue;
		}
		else if (endm == nrv)
		{
			if (Opt == 'T')
			{
				Si.push_back(0.); // total
			}
			else
			{
				Si.push_back(1.);   // main
			}
			printf("GSA i=%i, Si=%.2f, K=%i \n", nc + 1, Si[nc], Kos);
			continue;
		}
		mat data(endm + 1, nmc);

		for (int ne = 0; ne < endm; ne++)
		{
			int idx = combs[nc][ne];

			if (idx > nrv - 1) {
				std::string errMsg = "Error running UQ engine: combination index exceeds the bound";
				theErrorFile.write(errMsg);
			}

			for (int ns = 0; ns < nmc; ns++)
			{
				data(ne, ns) = xval[ns][idx];
			}
		}


		for (int ns = 0; ns < nmc; ns++)
		{
			data(endm, ns) = gval[ns];
		}

		gmm_full model;
		bool status = model.learn(data, Kos, maha_dist, random_subset, 30, 100, V *1.e-3, false);

		if (status == false)
		{
			std::string errMsg = "Error running UQ engine: GSA learning failed";
			theErrorFile.write(errMsg);
		}

		if (Kos == 0)
		{
			std::string errMsg ="Error running UQ engine: GSA learning failed. Try with more number of samples.";
			theErrorFile.write(errMsg);
		}

		mat mu = model.means;   //nrv x Ko
		cube cov = model.fcovs; //nrv x nrv x Ko
		rowvec pi = model.hefts;   //1 x Ko 
		rowvec mug = mu.row(endm);    //1 x Ko

		vector<double> mui;
		mui.reserve(nmc);
		// Used to calculate conditional mean and covariance
		cube SiginvSig(1, endm, Kos);
		mat muk(endm, Kos);
		for (int k = 0; k < Kos; k++)
		{
			mat Sig12 = cov.subcube(0, endm, k, endx, endm, k);
			mat Sig11 = cov.subcube(0, 0, k, endx, endx, k);
			muk.col(k) = mu.submat(0, k, endx, k);
			SiginvSig.slice(k) = solve(Sig11, Sig12).t();
		}

		//model.means.print("means:");
		//model.fcovs.print("fcovs:");

		for (int i = 0; i < nmc; i++)
		{
			rowvec pik_tmp(Kos, fill::zeros);
			colvec muki(Kos);
			mat xi = data.submat(0, i, endx, i);

			for (int k = 0; k < Kos; k++)
			{
				mat tmp = SiginvSig.slice(k);
				mat muval = muk.col(k);
				muki.subvec(k, k) = mug(k) + SiginvSig.slice(k) * (xi - muval);
				pik_tmp(k) = pi(k) * mvnPdf(xi, muval, cov.subcube(0, 0, k, endx, endx, k));

			}

			rowvec piki = pik_tmp / sum(pik_tmp);
			mat tmp = piki * muki;
			mui.push_back(tmp(0, 0));
		}

		Vi = calVar(mui);
		//Si.push_back(Vi / V);

		if (Opt == 'T')
		{
			Si.push_back(1 - Vi / V); // total
		}
		else
		{
			Si.push_back(Vi / V);   // main
		}

		printf("GSA i=%i, Si=%.2f, K=%i \n", nc + 1, Si[nc], Kos);

		if (isinf(Si[nc]) || isnan(Si[nc]))
		{
			return { -100 };
		}
	}

	return Si;
}


double runGSA::mvnPdf(mat x, mat mu, mat cov) 
{
	
	double n = size(x)(1);
	double sqrt2pi = std::sqrt(2 * PI);
	mat xmu = x - mu;
	mat quadform = xmu.t() * inv(cov) * xmu;
	double norm = std::pow(sqrt2pi, -n)*std::pow(abs(det(cov)), -0.5);
	//std::cout << norm << std::endl;

	return norm * std::exp(-0.5 * quadform(0,0));
}

double runGSA::calMean(vector<double> x) {
	double sum = std::accumulate(std::begin(x), std::end(x), 0.0);
	return sum / x.size();

}

runGSA::~runGSA() {};

double runGSA::calVar(vector<double> x) {
	double m = calMean(x);
	double accum = 0.0;
	std::for_each(std::begin(x), std::end(x), [&](const double d) {
		accum += (d - m) * (d - m);
		});
	//std::cout << (accum / (x.size())) << std::endl;
	return (accum / (x.size()));
}

void runGSA::writeTabOutputs(jsonInput inp, int procno)
{
	if (procno==0) {
		// dakotaTab.out
		std::string writingloc1 = inp.workDir + "/dakotaTab.out";
		std::ofstream Taboutfile(writingloc1);

		if (!Taboutfile.is_open()) {

			std::string errMsg = "Error running UQ engine: Unable to write dakotaTab.out";
			theErrorFile.write(errMsg);
		}

		Taboutfile.setf(std::ios::fixed, std::ios::floatfield); // set fixed floating format
		Taboutfile.precision(10); // for fixed format

		Taboutfile << "idx         ";
		for (int j = 0; j < inp.nrv + inp.nco + inp.nre; j++) {
			Taboutfile << inp.rvNames[j] << "           ";
		}
		for (int j = 0; j < inp.nqoi; j++) {
			Taboutfile << inp.qoiNames[j] << "            ";
		}
		Taboutfile << std::endl;


		for (int i = 0; i < inp.nmc; i++) {
			Taboutfile << std::to_string(i + 1) << "    ";
			for (int j = 0; j < inp.nrv + inp.nco + inp.nre; j++) {
				Taboutfile << std::to_string(xval[i][j]) << "    ";
			}
			for (int j = 0; j < inp.nqoi; j++) {
				Taboutfile << std::to_string(gval[i][j]) << "    ";
			}
			Taboutfile << std::endl;
		}
	}
}

void runGSA::writeOutputs(jsonInput inp, double dur, int procno)
{
	if (procno == 0) {
		// dakota.out
		string writingloc = inp.workDir + "/dakota.out";
		std::ofstream outfile(writingloc);

		if (!outfile.is_open()) {

			std::string errMsg = "Error running UQ engine: Unable to write dakota.out";
			theErrorFile.write(errMsg);

		}

		outfile.setf(std::ios::fixed, std::ios::floatfield); // set fixed floating format
		outfile.precision(4); // for fixed format

		outfile << "* number of input combinations" << std::endl;
		outfile << inp.ngr << std::endl;

		outfile << "* input names" << std::endl;
		for (int i = 0; i < inp.ngr; i++) {
			for (int j = 0; j < inp.groups[i].size() - 1; j++) {
				outfile << inp.rvNames[inp.groups[i][j]] << ",";
			}
			outfile << inp.rvNames[inp.groups[i][inp.groups[i].size() - 1]] << std::endl;
		}

		outfile << "* number of outputs" << std::endl;
		outfile << inp.nqoi << std::endl;

		outfile << "* output names" << std::endl;
		for (int i = 0; i < inp.nqoi; i++) {
			outfile << inp.qoiNames[i] << std::endl;
		}

		outfile << "* ";
		for (int j = 0; j < inp.ngr; j++) {
			//outfile << "Sm(" << std::to_string(j + 1) << ")  ";
			outfile << "Sm(";
			for (int k = 0; k < inp.groups[j].size() - 1; k++) {
				outfile << inp.rvNames[inp.groups[j][k]] << ",";
			}
			outfile << inp.rvNames[inp.groups[j][inp.groups[j].size() - 1]] << ") ";
		}
		for (int j = 0; j < inp.ngr; j++) {
			//outfile << "St(" << std::to_string(j + 1) << ")  ";
			outfile << "St(";
			for (int k = 0; k < inp.groups[j].size() - 1; k++) {
				outfile << inp.rvNames[inp.groups[j][k]] << ",";
			}
			outfile << inp.rvNames[inp.groups[j][inp.groups[j].size() - 1]] << ") ";
		}
		outfile << std::endl;

		for (int i = 0; i < inp.nqoi; i++) {

			for (int j = 0; j < inp.ngr; j++) {
				outfile << Simat[i][j] << " ";
			}
			for (int j = 0; j < inp.ngr; j++) {
				outfile << Stmat[i][j] << " ";
			}
			outfile << std::endl;
		}

		outfile << "* number of samples" << std::endl;
		outfile << inp.nmc << std::endl;

		outfile << "* elapsed time:" << std::endl;
		outfile << dur << " s" << std::endl;

		outfile.close();

	}

}