#include <iostream>
#include <time.h>
#include "fem.h"

extern void Gauss(double*,double*,double*,int);

int main()
{
  clock_t start, end;
  double cpuTime;
  start = clock();
  //----------------------------------------------
  // File handling
  std::ofstream aux_file;
  aux_file.open("cputime.txt");
  std::ifstream msh_file;
  msh_file.open("column.msh");
  std::ofstream msh_file_test;
  msh_file_test.open("column_test.msh");
  std::ofstream matrix_file_test;
  matrix_file_test.open("element_matrices.txt");
  if(!msh_file.good())               // Check is file existing
  {
    std::cout << "! Error in STD::Read: file could not be opened" << std::endl;
    return 0;
  }
  msh_file.seekg(0L,std::ios::beg);       // Rewind file
  //----------------------------------------------
  FEM* fem = new FEM();
  fem->ReadMesh(msh_file);
  fem->OutputMesh(msh_file_test);
  fem->SetInitialConditions();
  fem->SetBoundaryConditions();
  //----------------------------------------------
  int tn = 1;
  for(int t=0;t<tn;t++)
  {
    fem->CalculateElementMatrices();
    fem->DumpElementMatrices(matrix_file_test);
    fem->AssembleEquationSystem();
    fem->IncorporateBoundaryConditions();
    fem->DumpEquationSystem();
    Gauss(fem->matrix,fem->vecb,fem->vecx,fem->nn);
    fem->SaveTimeStep();
    fem->OutputResults(t);
  }
  //----------------------------------------------
  end = clock();
  cpuTime= (end-start)/ (double)(CLOCKS_PER_SEC);
  aux_file << "CPU time:" << cpuTime << std::endl;
  aux_file.close();
  fem->out_file.close();
  return 0;
}

