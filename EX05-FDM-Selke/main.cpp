#include <iostream>
#include <time.h>
#include "fdm.h"

int main(int argc, char *argv[])
{
  clock_t start, end;
  double cpuTime;
  std::ofstream aux_file;
  aux_file.open("cputime.txt");
  start = clock();
  //----------------------------------------------
  FDM* fdm = new FDM();
  fdm->SetActiveNodes();
  fdm->SetInactiveNodes();
  fdm->SetInitialConditions();
  fdm->SetBoundaryConditions();
  //----------------------------------------------
  int tn = 1000;
  for(int t=0;t<tn;t++)
  {
    fdm->RunTimeStep();
    fdm->OutputResults(t); //changed with save
    fdm->SaveTimeStep();
    if((t%100)==0)
      fdm->OutputResultsVTK(t);
  }
  //----------------------------------------------
  end = clock();
  cpuTime= (end-start)/ (double)(CLOCKS_PER_SEC);
  aux_file << "CPU time:" << cpuTime << std::endl;
  aux_file.close();
  return 0;
}
