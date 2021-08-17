
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
 *  Calcualtes the moments of QoIs and writes the results at dakotaTab.out
 */

#include "runForward.h"
#include <iterator>

runForward::runForward() {}
runForward::runForward(vector<vector<double>> xval,	vector<vector<double>> gmat, int procno)
{
	if (procno == 0) {

		this->xval = xval;
		this->gval = gmat;
		nmc = xval.size();
		nrv = xval[0].size();

		std::cout << "RV     Mean    StdDev  Skewness  Kurtosis" << std::endl;
		for (int nr = 0; nr < nrv; nr++) {
			vector<double> xvec;
			for (int ns = 0; ns < nmc; ns++) {
				xvec.push_back(xval[ns][nr]);
			}

			double mean_val = calMean(xvec);
			double stdDev_val = calStd(xvec, mean_val);
			double skewness_val = calSkewness(xvec, mean_val, stdDev_val);
			double kurtosis_val = calKurtosis(xvec, mean_val, stdDev_val);

			std::cout << "RV " << nr << ": ";
			std::cout << mean_val << " " << stdDev_val << " " << skewness_val << " " << kurtosis_val << std::endl;

			mean.push_back(mean_val);
			stdDev.push_back(stdDev_val);
			skewness.push_back(skewness_val);
			kurtosis.push_back(kurtosis_val);
		}	
	}
}

runForward::~runForward() {};

double runForward::calMean(vector<double> x) {
	double sum = std::accumulate(std::begin(x), std::end(x), 0.0);
	return sum / x.size();
}

double runForward::calStd(vector<double> x, double m) {
	double accum = 0.0;
	std::for_each(std::begin(x), std::end(x), [&](const double d) {
		accum += (d - m) * (d - m);
		});
	return std::sqrt(accum / (x.size()));
}

double runForward::calSkewness(vector<double> x, double m, double s) {
	double accum = 0.0;
	std::for_each(std::begin(x), std::end(x), [&](const double d) {
		accum += (d - m) * (d - m) * (d - m);
		});
	return (accum / (x.size()) / (s*s*s));
}

double runForward::calKurtosis(vector<double> x, double m, double s) {
	double accum = 0.0;
	std::for_each(std::begin(x), std::end(x), [&](const double d) {
		accum += (d - m) * (d - m) * (d - m) * (d - m);
		});
	return (accum / (x.size()) / (s * s * s * s));
	 
}
void runForward::writeOutputs(jsonInput inp, int procno)
{
	if (procno == 0) {

		// dakota.out
		string writingloc = inp.workDir + "/dakota.out";
		std::ofstream outfile(writingloc);


		if (!outfile.is_open()) {

			std::string errMsg = "Error running UQ engine: Unable to write dakota.out";
			theErrorFile.write(errMsg);
		}

		json outJson;

		outJson["rvNames"] = inp.rvNames;
		outJson["mean"] = mean;
		outJson["standardDeviation"] = stdDev;
		outJson["skewness"] = skewness;
		outJson["kurtosis"] = kurtosis;

		outfile << outJson.dump(4) << std::endl;
	}
}


void runForward::writeTabOutputs(jsonInput inp, int procno)
{
	if (procno==0) {
		// dakotaTab.out
		std::string writingloc1 = inp.workDir + "/dakotaTab.out";
		std::ofstream Taboutfile(writingloc1);

		if (!Taboutfile.is_open()) {

			std::string errMsg = "Error running UQ engine: Unable to write dakota.out";
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