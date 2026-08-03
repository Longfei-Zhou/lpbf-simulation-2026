#ifndef THESIS_MPI_STRUCTS_H
#define THESIS_MPI_STRUCTS_H

#include <mpi.h>
#include <string>

#include "DataStructs.h"
#include "Util.h"

using std::string;

class ThesisMPI{
private:
    MPI_Comm comm;

    int rank;
    int nproc;

    int i_min, i_max;
    int j_min, j_max;

public:
    ThesisMPI(MPI_Comm inComm){
        comm = inComm;

        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &nproc);

        name = Util::ZeroPadNumber(rank, 1+nproc/10);

        MPI_Dims_create(nproc, 2, dims);

        int periods[2] = {0,0}; // non-periodic boundaries
        MPI_Comm cart_comm;
        MPI_Cart_create(comm, 2, dims, periods, 0, &cart_comm); 

        MPI_Cart_coords(cart_comm, rank, 2, coords);  
    }

    int dims[2] = {0,0};    // Dimensionality of decomposition
    int coords[2] = {0,0};   // Own coordinates

    string name;

    int size() { return nproc; }

    void setPrint(Simdat& sim) {
        sim.print = rank == 0;
        sim.mpi = size() > 1;
    }

    void makeLocalBounds(Simdat& sim){
        const int I = dims[0];
        const int J = dims[1];
        
        const int i = coords[0];
        const int j = coords[1];

        if (sim.param.mode=="Stork" || sim.settings.mpi_overlap){ 
            i_min = ((sim.domain.xnum-1)*i)/I;
            i_max = ((sim.domain.xnum-1)*(i+1))/I;

            j_min = ((sim.domain.ynum-1)*j)/J;
            j_max = ((sim.domain.ynum-1)*(j+1))/J;
        }
        else{
            i_min = ((sim.domain.xnum-1)*i)/I + (i!=0);
            i_max = ((sim.domain.xnum-1)*(i+1))/I;

            j_min = ((sim.domain.ynum-1)*j)/J + (j!=0);
            j_max = ((sim.domain.ynum-1)*(j+1))/J;
        }

        const double xmin = sim.domain.xmin + sim.domain.xres*i_min;
        const double xmax = sim.domain.xmin + sim.domain.xres*i_max;
        
        const double ymin = sim.domain.ymin + sim.domain.yres*j_min;
        const double ymax = sim.domain.ymin + sim.domain.yres*j_max;

        sim.domain.xmin = xmin; sim.domain.xmax = xmax;
        sim.domain.xnum = 1 + int(0.5 + (sim.domain.xmax - sim.domain.xmin) / sim.domain.xres);
        sim.domain.xmax = sim.domain.xmin + (sim.domain.xnum - 1) * sim.domain.xres;
        
        sim.domain.ymin = ymin; sim.domain.ymax = ymax;
        sim.domain.ynum = 1 + int(0.5 + (sim.domain.ymax - sim.domain.ymin) / sim.domain.yres);
		sim.domain.ymax = sim.domain.ymin + (sim.domain.ynum - 1) * sim.domain.yres;

        sim.domain.pnum = sim.domain.xnum*sim.domain.ynum*sim.domain.znum;
    }
};


#endif // THESIS_MPI_STRUCTS_H
