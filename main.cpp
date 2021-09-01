
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
 *  Backend application to run SimCenterUQ engine in the quoFEM application developed by SimCenter. 
 */

#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <chrono>

#include "jsonInput.h"
#include "ERADist.h"
#include "ERANataf.h"
#include "runGSA.h"
#include "runForward.h"
#include <regex>

#include "writeErrors.h"
writeErrors theErrorFile; // Error log

#include <mpi.h>

int main(int argc, char** argv)
{
	MPI_Comm comm; int nprocs, procno;
	MPI_Init(&argc, &argv);
	comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &procno);
	MPI_Comm_size(comm, &nprocs);

    std::regex re(R"(\{([^}]+)\})");
    std::cerr << "testing\n";

	if ((argc != 4) && (procno == 0)) {
		std::string errMsg = "Number of the additional commend line arguments is " + std::to_string(argc - 1) + ", but 3 is required. The arguments should always include the working directory / os type / run type";
		std::cerr << errMsg << std::endl;
		theErrorFile.abort();
	}

	std::string workDir = argv[1];
	std::string osType = argv[2];
	std::string runType = argv[3];

	if (procno == 0) std::cerr << "* WORKDIR: " << workDir << "\n";
	if (procno == 0) std::cerr << "* OS: " << osType << "\n";
	if (procno == 0) std::cerr << "* RUN: " << runType << "\n\n";

	//std::string workDir = "C:/Users/yisan/Documents/quoFEM/LocalWorkDir/tmp.SimCenter";
	//std::string osType = "Windows";
	//std::string runType = "runningLocal";

	std::string  errorDir = workDir + "/dakota.err";
	theErrorFile.getFileName(errorDir, procno);

	//auto elapseStart = std::chrono::high_resolution_clock::now();
	auto elapseStart = std::chrono::high_resolution_clock::now();
	double elapsedTime;

	//
	//  (1) read JSON file
	//

	jsonInput inp(workDir, procno);

	//
	//	(2) Construct Nataf Object
	//

	ERANataf T(inp, procno);

	//
	//	(3-1) Random number generator Gaussian(mean=0,var=1) - (batch samples)
	//

	T.sample(inp, procno);

	//
	//	(4-1) FE Analysis - (parallel)
	//	
	
	T.simulateAppBatch(osType, runType, inp, procno, nprocs);

	//std::cout<<"Just testing this location 1\n";


	//
	//	(5) Run UQ analysis
	//

	if (!inp.uqType.compare("Sensitivity Analysis")) {

		//
		//	(5-1) Global sensitivity analysis - (parallel)
		//

		runForward ForwardResults(T.X, T.G, procno);
		ForwardResults.writeTabOutputs(inp, procno); 	//	Write dakotaTab.out

		runGSA GsaResults(T.X, T.G, inp.groups, 25, procno, nprocs); // int Kos = 25;

		//elapsedTime = (std::chrono::high_resolution_clock::now() - elapseStart).count();
		elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - elapseStart).count()/1.e3;
		GsaResults.writeOutputs(inp, elapsedTime, procno); //	Write dakota.out
	}
	else if (!inp.uqType.compare("Forward Propagation")) {
		//
		//	(5-2) Forward analysis
		//
		runForward ForwardResults(T.X, T.G, procno);
		ForwardResults.writeTabOutputs(inp, procno); 	//	Write dakotaTab.out
		// = (std::chrono::high_resolution_clock::now() - elapseStart).count();
		elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - elapseStart).count() /1.e3;
		ForwardResults.writeOutputs(inp, procno);		//	Write dakota.out <- curretly not being used anywhere
	}
	theErrorFile.close();
	if (procno == 0) std::cout << "Elapsed time: " << elapsedTime << " s\n";
	MPI_Finalize();
	return 0;
}

