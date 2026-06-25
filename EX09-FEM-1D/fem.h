#ifndef FEM_H
#define FEM_H

#include <vector>
#include <list>
#include <fstream>
#include <iostream>

class FEM
{
private:
public:
  //data structures
  //std::ifstream msh_file;
  std::vector<double>node_vector;
  std::vector<int*>element_vector;
  int *element_nodes;
  int nn,ne;
  double* element_matrix;
  std::vector<double*>element_matrix_vector;
  double h_initial;
  double h_top;
  double h_bottom;
  double x0,x1,L;
  double K[4];
  int dof;

  std::vector<double>u_new;
  std::vector<double>u;
  std::vector<double>u_bc;
  double u0;
  double S0,Kf,Q;
  std::ofstream out_file;
  std::vector<int>bc_nodes;
  std::vector<int>nodes_inactive;
  double* matrix;
  double* vecb;
  double* vecx;
public:
    FEM();
    void SetInitialConditions();
    void SetBoundaryConditions();
    void RunTimeStep();
    void SaveTimeStep();
    void OutputResults(int);
    void AssembleEquationSystem();
    void DumpEquationSystem();
    void IncorporateBoundaryConditions();
    void CalculateElementMatrices();
    void CalculateElementMatrix(int);
    void ReadMesh(std::ifstream&);
    void OutputMesh(std::ofstream&);
    void DumpElementMatrices(std::ofstream&);
};

#endif // FEM_H
