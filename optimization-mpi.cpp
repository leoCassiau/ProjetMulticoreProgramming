/*
  Branch and bound algorithm to find the minimum of continuous binary 
  functions using interval arithmetic.

  MPI parallel version

  Authors: Cassiau Léo, Ugo Mahey and Frederic Goualard <Frederic.Goualard@univ-nantes.fr>
*/

#include <iostream>
#include <iterator>
#include <string>
#include <stdexcept>
#include "interval.h"
#include "functions.h"
#include "minimizer.h"
#if _OPENMP
#   include <omp.h>
#endif
#include <mpi.h>
#include <string.h>

using namespace std;

// Rank of actual processor and number of processors. MPI variables
int rang, numprocs;

// Next processor
int cpt = 1;

// Name of the function to optimize
char choice_fun[50];

// Split a 2D box into four subboxes by splitting each dimension
// into two equal subparts
void split_box(const interval& x, const interval& y,
	       interval &xl, interval& xr, interval& yl, interval& yr)
{
  double xm = x.mid();
  double ym = y.mid();
  xl = interval(x.left(),xm);
  xr = interval(xm,x.right());
  yl = interval(y.left(),ym);
  yr = interval(ym,y.right());
}

// Branch-and-bound minimization algorithm
void minimize(itvfun f,  // Function to minimize
	      const interval& x, // Current bounds for 1st dimension
	      const interval& y, // Current bounds for 2nd dimension
	      double threshold,  // Threshold at which we should stop splitting
	      double& min_ub,  // Current minimum upper bound
	      minimizer_list& ml) // List of current minimizers
{
  interval fxy = f(x,y);
     
  if (fxy.left() > min_ub) { // Current box cannot contain minimum?
    return ;
  }

  if (fxy.right() < min_ub) { // Current box contains a new minimum?
    min_ub = fxy.right();
    // Discarding all saved boxes whose minimum lower bound is 
    // greater than the new minimum upper bound
    auto discard_begin = ml.lower_bound(minimizer{0,0,min_ub,0});
	
		#pragma  omp critical 
		ml.erase(discard_begin,ml.end());
  }

  // Checking whether the input box is small enough to stop searching.
  // We can consider the width of one dimension only since a box
  // is always split equally along both dimensions
  if (x.width() <= threshold) { 
    // We have potentially a new minimizer
	 	#pragma  omp critical 
    ml.insert(minimizer{x,y,fxy.left(),fxy.right()});
    
    return ;
  }

  // The box is still large enough => we split it into 4 sub-boxes
  // and recursively explore them
  interval xl, xr, yl, yr;
  split_box(x,y,xl,xr,yl,yr);

  #pragma omp parallel  
	#pragma omp sections
	{
		#pragma omp section
		minimize(f,xl,yl,threshold,min_ub,ml);
		#pragma omp section
		minimize(f,xl,yr,threshold,min_ub,ml);
		#pragma omp section
		minimize(f,xr,yl,threshold,min_ub,ml);
		#pragma omp section
		minimize(f,xr,yr,threshold,min_ub,ml);
  }
}

void MPI_Send_Interval(const interval & x, const interval & y, double threshold, double min_ub) {
	double envoie_inter_thres_min[6] = {x.left(),x.right(),y.left(),y.right(),threshold,min_ub};
	MPI_Send(&choice_fun,50,MPI_CHAR,cpt,0,MPI_COMM_WORLD);
	MPI_Send(&envoie_inter_thres_min,6,MPI_DOUBLE,cpt,0,MPI_COMM_WORLD);	
}

void MPI_Recv_Interval(itvfun & f, interval & x, interval & y, double& threshold, double & min_ub)  {
	opt_fun_t fun;
	char recv_fun[50];
	double recv_inter_thres_min[6];

	MPI_Recv(&recv_fun,50,MPI_CHAR,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Reception du message contenant le nom de la fonction
	MPI_Recv(&recv_inter_thres_min,6,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Reception des intervales, de la précision et du minimum actuel

	x = interval(recv_inter_thres_min[0],recv_inter_thres_min[1]);
	y = interval(recv_inter_thres_min[2],recv_inter_thres_min[3]);
	threshold = recv_inter_thres_min[4];
	min_ub = recv_inter_thres_min[5];

  fun = functions.at(recv_fun);
	f = fun.f;
}

// Branch-and-bound minimization algorithm
void minimize_mpi(itvfun f,  // Function to minimize
	      const interval& x, // Current bounds for 1st dimension
	      const interval& y, // Current bounds for 2nd dimension
	      double threshold,  // Threshold at which we should stop splitting
	      double& min_ub,  // Current minimum upper bound
	      minimizer_list& ml) // List of current minimizers
{
  interval fxy = f(x,y);
     
  if (fxy.left() > min_ub) { // Current box cannot contain minimum?
    return ;
  }

  if (fxy.right() < min_ub) { // Current box contains a new minimum?
    min_ub = fxy.right();
    // Discarding all saved boxes whose minimum lower bound is 
    // greater than the new minimum upper bound
    auto discard_begin = ml.lower_bound(minimizer{0,0,min_ub,0});
	
		#pragma  omp critical 
		ml.erase(discard_begin,ml.end());
  }

  // Checking whether the input box is small enough to stop searching.
  // We can consider the width of one dimension only since a box
  // is always split equally along both dimensions
  if (x.width() <= threshold) { 
    // We have potentially a new minimizer
	 	#pragma  omp critical 
    ml.insert(minimizer{x,y,fxy.left(),fxy.right()});
    
    return ;
  }

  // The box is still large enough => we split it into 4 sub-boxes
  // and recursively explore them
  interval xl, xr, yl, yr;
  split_box(x,y,xl,xr,yl,yr);

	if(cpt < numprocs) {
		MPI_Send_Interval( xl, yl, threshold, min_ub);
		++cpt;
	} else {
		minimize(f,xl,yl,threshold,min_ub,ml);
	}
	
	if(cpt < numprocs) {
		MPI_Send_Interval( xl, yr,threshold,min_ub);
		++cpt;
	} else {
		minimize(f,xl,yr,threshold,min_ub,ml);
	}

	if(cpt < numprocs) {
		MPI_Send_Interval(xr, yl,threshold,min_ub);
		++cpt;
	} else {
		minimize(f,xr,yl,threshold,min_ub,ml);
	}

	if(cpt < numprocs) {
		minimize_mpi(f,xr,yr,threshold,min_ub,ml);
	} else {
		minimize(f,xr,yr,threshold,min_ub,ml);
	}
}


int main(int argc, char *argv[])
{
  cout.precision(16);
  // By default, the currently known upper bound for the minimizer is +oo
  double min_ub = numeric_limits<double>::infinity();
  // List of potential minimizers. They may be removed from the list
  // if we later discover that their smallest minimum possible is 
  // greater than the new current upper bound
  minimizer_list minimums;

  // Threshold at which we should stop splitting a box
  double precision;

	// MPI Initialization
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rang);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	if(rang == 0) {
		// The information on the function chosen (pointer and initial box)
		opt_fun_t fun;
		
		bool good_choice;
		// Asking the user for the name of the function to optimize
		do {
		  good_choice = true;

		  cout << "Which function to optimize?\n";
		  cout << "Possible choices: ";
		  for (auto fname : functions) {
		    cout << fname.first << " ";
		  }
		  cout << endl;
		  cin >> choice_fun;
		 	//strcpy(choice_fun, "goldstein_price");

		  try {
		    fun = functions.at(choice_fun);
		  } catch (out_of_range) {
		    cerr << "Bad choice" << endl;
		    good_choice = false;
		  }
		} while(!good_choice);

		// Asking for the threshold below which a box is not split further
		cout << "Precision? ";
		cin >> precision;
		//precision = 0.007;

  	minimize_mpi(fun.f,fun.x,fun.y,precision,min_ub,minimums);
	} else {
		itvfun f;
		interval x, y;
		MPI_Recv_Interval(f, x, y, precision, min_ub);

		minimize(f, x, y, precision, min_ub,minimums);
	}

	// Combining min_ub
	double total_min_ub;
	MPI_Reduce( &min_ub,&total_min_ub , 1, MPI_DOUBLE, MPI_MIN, 0,  MPI_COMM_WORLD);
	
	if(rang == 0) {
		// Displaying all potential minimizers
		/*copy(minimums.begin(),minimums.end(),
		     ostream_iterator<minimizer>(cout,"\n"));
		cout << "Number of minimizers: " << minimums.size() << endl;*/
		cout << "Upper bound for minimum: " << total_min_ub << endl;
	}

	MPI_Finalize();
}
