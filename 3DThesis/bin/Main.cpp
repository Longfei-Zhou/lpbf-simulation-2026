/****************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC                                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of 3dThesis. 3dThesis is distributed under a           *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <vector>
#include <cmath>
#include <chrono>
#include <iostream>

#include "DataStructs.h"
#include "Init.h"
#include "Util.h"
#include "Run.h"
#include "Out.h"
#include "Grid.h"
#include "ThesisConfig.h"

#ifdef Thesis_ENABLE_MPI
#include "MpiStructs.h"
#include <mpi.h>
#endif

using std::vector;
using std::string;

using std::chrono::high_resolution_clock;
using std::chrono::duration;

inline void run(int argc, char * argv[])
{
	auto start_in = high_resolution_clock::now();

	if (argc <= 1) { throw std::runtime_error( "Input file argument required: ./3DThesis ParamInput.txt" ); }
        string inputFile = argv[1];

	Simdat sim;

#ifdef Thesis_ENABLE_MPI
	ThesisMPI mpi(MPI_COMM_WORLD);
	mpi.setPrint(sim);
#endif

	Init::GetFileNames(sim.files, inputFile, sim.print);	
	
	Init::ReadSimParams(sim);

#ifdef Thesis_ENABLE_MPI
	mpi.makeLocalBounds(sim);
#endif

	Grid grid(sim); 

	sim.util.approxEndTime = sim.util.allScansEndTime;

	auto stop_in = high_resolution_clock::now();
	if (sim.print)
		std::cout << "Initialization time (s): " << (duration<double, std::milli>(stop_in - start_in).count())/1000.0 << "\n\n";

	auto start_sim = high_resolution_clock::now();

	Run::Simulate(grid, sim);

	auto stop_sim = high_resolution_clock::now();
	if (sim.print) {
                std::cout << "Version: " << Out::version() << "\n";
                std::cout << "Commit hash: " << Out::commitHash() <<	"\n";
                std::cout << "Execution time (s): " << (duration<double, std::milli>(stop_sim - start_sim).count())/1000.0 << "\n\n";
	}
	auto start_out = high_resolution_clock::now();

std::string rank_name = "";
#ifdef Thesis_ENABLE_MPI
if (sim.mpi)
	rank_name = "." + mpi.name;
#endif

	if (sim.param.mode=="Solidification"){ 
		if (sim.param.tracking!="Stork"){
			grid.Output(sim, "Solidification.Final" + rank_name); 
		}
		else{
			grid.Output_RRDF_csv(sim, "RRDF" + rank_name);
		}
	}
	if (sim.output.T_hist) { 
		grid.Output_T_hist(sim, "T.hist" + rank_name); 
	}
	if (sim.output.RDF) { 
		grid.Output_RDF(sim, "RDF.Final" + rank_name); 
	}

	auto stop_out = high_resolution_clock::now();
	if (sim.print)
		std::cout << "Output time (s): " << (duration<double, std::milli>(stop_out - start_out).count()) / 1000.0 << "\n\n";
}

#ifdef Thesis_ENABLE_MPI
int main(int argc, char * argv[]) {	
    MPI_Init(&argc, &argv);
	run(argc, argv);
    MPI_Finalize();

	return 0;
}
#else
int main(int argc, char * argv[]) {	
	run(argc, argv);
	return 0;
}
#endif
