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
	
	int n, order, gc, i,i_max,j;
	double x_low,x_high,cfl,dt,dx;
	string type;
	// *** INPUTS ***
	n = 10;                         // Number of Cells
	order = 4;                      // Order of Accuracy
	gc = order/2;                   // Find Number of GC
	x_low = 0, x_high = 6.28318531; // Range
	cfl = 0.01;						// cfl
	dt = 1.0/n*cfl;                 // Timestep
	i_max = 80000;                   // Total Number of Iterations
	type = "nonlinear_convective";

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
		j = 0;
		for (int iter = 0; iter < i_max; iter++){
			
			dudx = DDx(u,&dx,&gc,&n,order);

			for ( i = 0; i < n + 2*gc; i++){
				if(type =="linear_convective"){
					unew[i] = u[i] - dt*dudx[i];
					u[i] = unew[i];
				} else if(type =="nonlinear_convective"){
					unew[i] = u[i] - u[i]*dt*dudx[i];
					u[i] = unew[i];
				}
				
			}
			u = fill_gc(u,gc,n);

			if (iter%80 == 0){         //Taking snapshots of runs
				
				string s = to_string(j);   // Organizing files names so matlab can read them in order
				if (j < 10) {
					s.insert(0,1,'0');
					s.insert(0,1,'0');
					s.insert(0,1,'0');
					s.insert(0,1,'0');
				}else if (j<100){
					s.insert(0,1,'0');
					s.insert(0,1,'0');
					s.insert(0,1,'0');
				}
				else if (j<1000){
					s.insert(0,1,'0');
					s.insert(0,1,'0');
				}
				else if (j<10000){
					s.insert(0,1,'0');
				}
				filename = "Iter_";
				filename.append(s);
				filename.append(".csv");
				create(filename,x,u,n);
				j = j + 1;
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
