#include "fem.h"
#include <cmath>

FEM::FEM()
{
  h_initial = 5.;
  h_top = 12.;
  h_bottom = 0.;
  K[0] = 1e-5;
  K[1] = 2e-5;
  K[2] = 1e-5;
  K[3] = 1e-5;
  dof = 5;
  u.resize(dof);
  u_new.resize(dof);
  out_file.open("out.txt");
  matrix = new double[dof*dof];
  vecb = new double[dof];
  vecx = new double[dof];
}

void FEM::ReadMesh(std::ifstream& msh_file)
{
  std::string input_line;
  char buffer[256]; // MAX_LINE
  std::ios::pos_type position;
  double x;
  while(!msh_file.eof())
  {
    position = msh_file.tellg();
    msh_file.getline(buffer,256);
    input_line = buffer;
    if(input_line.size()<1) // empty line
      continue;
    // Dealing with keywords
    if(input_line.find("#STOP")!=std::string::npos) // keyword found
      return; // position;
    // Dealing with subkeywords
    if(input_line.find("$NODES")!=std::string::npos)
    {
      msh_file >> nn;
      for(int i=0;i<nn;i++)
      {
        msh_file >> x;
        node_vector.push_back(x);
      }
    }
    if(input_line.find("$ELEMENTS")!=std::string::npos)
    {
      msh_file >> ne;
      for(int i=0;i<ne;i++)
      {
        element_nodes = new int[2];
        msh_file >> element_nodes[0];
        msh_file >> element_nodes[1];
        element_vector.push_back(element_nodes);
      }
    }
  }
}

void FEM::OutputMesh(std::ofstream& msh_file_test)
{
  //-----------------------------------------------------------------------
  msh_file_test << "#FEM_MSH" << std::endl;
  msh_file_test << "$NODES" << std::endl;
  for(int n=0;n<(int)node_vector.size();n++)
  {
    msh_file_test << node_vector[n] << std::endl;
  }
  msh_file_test << "$ELEMENTS" << std::endl;
  for(int e=0;e<(int)element_vector.size();e++)
  {
    msh_file_test << element_vector[e][0] << " " << element_vector[e][1] << std::endl;
  }
  msh_file_test << "#STOP" << std::endl;
}

void FEM::SetInitialConditions()
{
  for(int n=0;n<(int)node_vector.size();n++)
  {
    u[n] = h_initial;
    u_new[n] = h_initial;
  }
}

void FEM::SetBoundaryConditions()
{
  bc_nodes.push_back(0); u[bc_nodes[0]] = h_top;  u_new[bc_nodes[0]] = h_top;
  bc_nodes.push_back(4); u[bc_nodes[1]] = h_bottom;  u_new[bc_nodes[1]] = h_bottom;
}

void FEM::CalculateElementMatrices()
{
  for(int e=0;e<(int)element_vector.size();e++)
  {
    CalculateElementMatrix(e);
  }
}

void FEM::CalculateElementMatrix(int e)
{
  element_matrix = new double[4];
  element_nodes = element_vector[e];
  x0 = node_vector[element_nodes[0]];
  x1 = node_vector[element_nodes[1]];
  L = x1-x0;
  element_matrix[0] = K[e]/L;
  element_matrix[1] = -K[e]/L;
  element_matrix[2] = -K[e]/L;
  element_matrix[3] = K[e]/L;
  element_matrix_vector.push_back(element_matrix);
}

void FEM::DumpElementMatrices(std::ofstream& file)
{
  int ii;
  for(int e=0;e<(int)element_matrix_vector.size();e++)
  {
    file << "--------------------------" << std::endl;
    element_matrix = element_matrix_vector[e];
    for(int j=0;j<2;j++)
    {
      for(int i=0;i<2;i++)
      { 
        ii= 2*j+i;
        file << element_matrix[ii] << " ";
      }
      file << std::endl;
    }
  }
}

void FEM::AssembleEquationSystem()
{
  // Initialize matrix
  for(int i=0;i<nn;i++)
  {
    vecx[i] = u[i];
    vecb[i] = 0.0;
    for(int j=0;j<nn;j++)
    {
      matrix[i*nn+j] = 0.0;
    }
  }
  // Matrix entries
  int i,j,ij;
  for(int e=0;e<ne;e++)
  {
    element_nodes = element_vector[e];
    element_matrix = element_matrix_vector[e];
    i = element_nodes[0]; //0
    j = element_nodes[1]; //1
    ij = i*nn+i; //0
    //out_file << ij << " ";
    matrix[ij] += element_matrix[0];
    ij = i*nn+i+1; //1
    //out_file << ij << " ";
    matrix[ij] += element_matrix[1];
    ij = j*nn+j-1;//5
    //out_file << ij << " ";
    matrix[ij] += element_matrix[2];
    ij = j*nn+j;//6
    //out_file << ij << std::endl;
    matrix[ij] += element_matrix[3];
  }
}

void FEM::SaveTimeStep()
{
  for(int n=0;n<(int)node_vector.size();n++)
  {
    u_new[n] = vecx[n];
    u[n] = u_new[n];
  }
}

void FEM::OutputResults(int t)
{
  if((t%1)==0)
  {
    for(int n=0;n<(int)node_vector.size();n++)
    {
      out_file << u_new[n] << std::endl;
    }
  }
}

void FEM::DumpEquationSystem()
{
  for(int i=0;i<nn;i++)
  {
    for(int j=0;j<nn;j++)
    {
      out_file << matrix[i*nn+j] << "\t";
    }
    out_file << "b:" << vecb[i] << " ";
    out_file << std::endl;
  }
}

void FEM::IncorporateBoundaryConditions()
{
  //------------------------------------
  size_t i_bc;
  int i_row, k;
  for(i_bc=0;i_bc<bc_nodes.size();i_bc++)
  {
     i_row = bc_nodes[i_bc];
     // Null off-diangonal entries of the related row and columns
     // Apply contribution to RHS by BC
     for(int j=0;j<nn;j++)
     {
        if(i_row == j)          
           continue; // do not touch diagonals
        matrix[i_row*(nn)+j] = 0.0; // NUll row
        k = j*(nn)+i_row;
        // Apply contribution to RHS by BC
        vecb[j] -= matrix[k]*u[i_row];
        matrix[k] = 0.0; // Null column
     }
     // Apply Dirichlet BC
     vecb[i_row] = u[i_row]*matrix[i_row*(nn)+i_row];
  } 
}

