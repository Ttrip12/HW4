#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

using namespace std;

vector<double> DDx(vector<double>, double *, int *,int *,int);
vector<double> DDxDDx(vector<double>, double *, int *,int *,int);
vector<double> fill_gc(vector<double>, int,int);
vector<double> error(vector<double>,vector<double>);
void create(string,vector<double>,vector<double>,int);

int main(){
	
	int n, order, gc, i,i_max;
	double x_low,x_high,cfl,dt,dx;
	
	// *** INPUTS ***
	n = 10;                         // Number of Cells
	order = 2;                      // Order of Accuracy
	gc = order/2;                   // Find Number of GC
	x_low = 0, x_high = 6.28318531; // Range
	cfl = 0.1;						// cfl
	dt = 1.0/n*cfl;                 // Timestep
	i_max = 1000;                   // Total Number of Iterations


	// ***** Define Mesh *****
	dx = (x_high - x_low)/n;

	vector<double> u(n + 2*gc), unew(n + 2*gc), x(n + 2*gc),dudx(n + 2*gc);
	// ***** Initial Conditions *****	
		
		for ( i = 0; i < n + 2*gc; i++) {
			
			x[i] = x_low - gc*dx + i*dx + dx/2; 
			u[i] = sin(x[i]) + 2;

		}
		
		// Create sub directory for csv
		std::filesystem::create_directory("Data");
		std::filesystem::current_path("Data");
		string filename = "Iter_";

	// Solve
		for (int iter = 0; iter < i_max; iter++){
			
			dudx = DDx(u,&dx,&gc,&n,order);

			for ( i = 0; i < n + 2*gc; i++){
				unew[i] = u[i] - dt*dudx[i];
				u[i] = unew[i];
			}
			u = fill_gc(u,gc,n);

			if (iter%10 == 0){
				string s = to_string(iter);
				filename = "Iter_";
				filename.append(s);
				filename.append(".csv");
				create(filename,x,u,n);
			}
			
			

		}


	return 0;

}

vector<double> DDx(vector<double> A, double* dx,int* gc,int* n,int order)
{
    vector<double> derivative(*n + 2 * *gc);

		if (order == 2) {
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+1] - A[i-1])/2.0/ *dx;
	
		} 
		} else if (order == 4) {
		
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (-A[i+2] + 8*A[i+1] - 8*A[i-1] + A[i-2])/12.0/ *dx;
		} 
			
		} else if (order == 6) {

		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+3]/60 - 3*A[i+2]/20 + 3*A[i+1]/4 - 3*A[i-1]/4 + 3*A[i-2]/20 - A[i-3]/60)/ *dx;
		} 
		}

	return derivative;
}

vector<double> DDxDDx(vector<double> A, double* dx,int* gc,int* n,int order)
{
    vector<double> derivative(*n + 2 * *gc);

		if (order == 2) {
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+1] -2*A[i] + A[i-1])/ *dx/ *dx;
	
		} 
		} else if (order == 4) {
		
		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (-A[i+2]/12.0 + 4.0*A[i+1]/3.0 - 5.0*A[i]/2.0 + 4*A[i-1]/3 - A[i-2]/12.0)/ *dx/ *dx;
	
		} 
			
		} else if (order == 6) {

		for (int i = 0; i < *n + 2 * *gc; i++) {

		derivative[i] = (A[i+3]/90 - 3*A[i+2]/20 + 3*A[i+1]/2 - 49*A[i]/18 + 3*A[i-1]/2 - 3*A[i-2]/20 + A[i-3]/90)/ *dx/ *dx;
	
		} 
		}

	return derivative;
}

vector<double> fill_gc(vector<double> A, int gc,int n){
	
	int sz = n + 2*gc;
	for(int i = 0; i < gc; i++){
		A[i] = A[n - 1 + i];
		A[sz - i - 1] = A[2*gc - i];
	}	 
 
	return A;
}

vector<double> error(vector<double> Numer,vector<double> Exact){
	
	vector<double> error(Numer.size());
	
	for( int i = 0; i < Numer.size(); i++){
		error[i] = abs(Numer[i] - Exact[i]);
	}
	return error;
}

void create(string file, vector<double> x, vector<double> y,int n) 
{ 

	
    // file pointer 
    fstream fout; 

    // opens an existing csv file or creates a new file. 
    fout.open(file, ios::out | ios::app); 
  
    
  
    int i;
  
    // Read the input 
    for (i = 0; i < n; i++) { 
  
        fout << x[i] << ", "
            << y[i] << ", "
            <<"\n"; 
  
    } 
}
