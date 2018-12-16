# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa047.hpp"

#include <armadillo>

using namespace std;

//global variables
arma::Row<double> x;	//vector-argument of objective function
int n = 2;				//number of vector x components

arma::Row<double> step_size;

//objective function
double objective_fn ( double z[]);


//****************************************************************************80

void minimize ( )
{
  timestamp ( );

  int i;
  int icount;
  int ifault;
  int kcount;
  int konvge;
  //int n;
  int numres;
  double reqmin;
  double *start;
  double *step;
  double *xmin;
  double ynewlo;

  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  //assigning values of arma::x to double x[] 
  //				and arma::step_size to double step[]
  //				to be used in nelmin()
  x.print("Starting point x0:");
  for (i = 0; i < n; i++)
  {
	  start[i] = x(i);
	  step[i] = step_size(i);
  }

  reqmin = 1.0E-08;
  konvge = 10;
  kcount = 1000;
  ynewlo = objective_fn(start);
  cout << "  F(x0) = " << ynewlo << "\n";

  nelmin(objective_fn, n, start, xmin, &ynewlo, reqmin, step,
	  konvge, kcount, &icount, &numres, &ifault);

  // assigning optimal vector xmin[] to arma::x (back to normal...)
  for (i = 0; i < n; i++) x(i) = xmin[i];
  x.print("\nResult x*:");
  cout << "Return code IFAULT = " << ifault << "\n";
  cout << "Resulting  F(x*) = " << ynewlo << "\n";
  cout << "  Number of iterations = " << icount << "\n";
  cout << "  Number of restarts =   " << numres << "\n";
  delete[] start;
  delete[] step;
  delete[] xmin;
  cout << "ASA047_PRB:\n";
  cout << "  Normal end of execution.\n";
  timestamp ( );
}
//****************************************************************************80
double objective_fn ( double z[])
{
	int i;
	double fx;
	double fx1;
	double fx2;
  // assign double z[] which comes from nelmin() to arma x
  // which is global and is used as argument of the objective function
	for (i = 0; i < n; i++) x(i) = z[i];
	fx1 = x(1) - x(0) * x(0) ;
	fx2 = 1.0 - x(0);
	fx = 100.0 * fx1 * fx1 +  fx2 * fx2;
//	cout << endl << fx;
	return fx;
}
//****************************************************************************80

int main()
{
	x.set_size(n);
	x.fill(0);
	x(0) = -0.05;
	x(1) = 5.5;

	step_size.set_size(n);
	step_size.fill(1);			//for now

	minimize();
	return 0;
}