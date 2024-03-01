#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <chrono>

using namespace std;
std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

vector<double> DDx(vector<double>, double *, int *,int *,int);
vector<double> DDxDDx(vector<double>, double *, int *,int *,int);
vector<double> fill_gc(vector<double>, int,int);
vector<double> error(vector<double>,vector<double>);
double find_dt(vector<double>,double,double,double);
void create(string,vector<double>,vector<double>,int);
int main(){
	
	int n, order, gc, i,i_max,j;
	double x_low,x_high,cfl,dt,dx,nu;
	string type;
	ifstream input;
	// *** INPUTS ***
	input.open("Input.txt", ios::in);
	string line;

		//while (getline(input,line)) {
		while (1){
			// getline(input,line);
			// if(line.find("#") == 0) {
			// 	continue;
			// }
			input >> n;
			input >> order;
			input >> x_low;
			input >> x_high;
			input >> cfl;
			input >> type;
			input >> nu;
			if (input.eof()){
				break;
			}
			
		}
	
	gc = order/2;		            // Find Number of GC
	dt = 1.0/n*cfl;                 // Timestep
	i_max = n*10;                   // Total Number of Iterations

	// ***** Define Mesh *****
	dx = (x_high - x_low)/n;

	vector<double> u(n + 2*gc), unew(n + 2*gc), x(n + 2*gc), dudx(n + 2*gc), d2u_d2x(n + 2*gc);
	// ***** Initial Conditions *****	
		
		for ( i = 0; i < n + 2*gc; i++) {
			
			x[i] = (x_low - gc*dx + i*dx + dx/2); 
			
			// u[i] = sin(x[i]*x[i]);
			
			// if (i<(n + 2*gc)/2){
			// 	u[i] = 0.5*x[i]*x[i]; 
			// } else u[i] = u[i-1];
			
			if (i<(n + 2*gc)/2){
				u[i] = x[i]+1; 
			} else u[i] = x_high - x[i]+1;

			// if (i<(n + 2*gc)/2){
			// 	u[i] = x[i]; 
			// } else u[i] = u[((n + 2*gc)/2)-1] - x[n + 2*gc - i];

		}
		
		// Create sub directory for csv
		std::filesystem::remove_all("Data");
		std::filesystem::create_directory("Data");
		std::filesystem::current_path("Data");
		string filename = "Iter_";

	// Solve
		j = 0;
		int num = i_max/1000;
		for (int iter = 0; iter < i_max; iter++){
			
			dudx = DDx(u,&dx,&gc,&n,order);
			d2u_d2x = DDxDDx(u,&dx,&gc,&n,order);

				if(type == "linear_convective"){
					
					for ( i = 0; i < n + 2*gc; i++){
						unew[i] = u[i] - dt*dudx[i];
						u[i] = unew[i];
					}
					
				} else if(type == "nonlinear_convective"){
					
					dt =  find_dt(u,dx,cfl,nu);
					for ( i = 0; i < n + 2*gc; i++){
						unew[i] = u[i] - u[i]*dt*dudx[i];
						u[i] = unew[i];
					}

				} else if (type == "burger"){
					dt =  find_dt(u,dx,cfl,nu);
					for ( i = 0; i < n + 2*gc; i++){
						unew[i] = u[i] - dt*u[i]*dudx[i] + dt*nu*d2u_d2x[i];
						u[i] = unew[i];
					}
				}

			u = fill_gc(u,gc,n);

			if (iter%num == 0){         //Taking snapshots of runs
				
				string s = to_string(j);   // Organizing files names so matlab can read them in order
				if (j < 10) {
					s.insert(0,1,'0');
					s.insert(0,1,'0');
					s.insert(0,1,'0');
				}else if (j<100){
					s.insert(0,1,'0');
					s.insert(0,1,'0');
				}
				else if (j<1000){
					s.insert(0,1,'0');
				}

				filename = "Iter_";
				filename.append(s);
				filename.append(".csv");
				create(filename,x,u,n+2*gc);
				j = j + 1;
			}
			
		}

		std::filesystem::current_path("..");
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		double time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		cout << "Time Elapsed = " << time_elapsed << "[ms]" << endl;
		string file = "out.csv";
		fstream fout;
		fout.open(file, ios::out | ios::app);
		fout << time_elapsed << ",";
	
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
		A[i] = A[n + i];
		A[sz - i - 1] = A[2*gc - i - 1];
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

double find_dt(vector<double> u,double dx,double cfl,double nu) {
	double dt;
	double u_max = abs(*max_element(u.begin(),u.end()));
	double dtvelo = dx*cfl/u_max;
	double dtvisc = (dx*dx*cfl)/nu;
	
	if (dtvisc > dtvelo){
		dt = dtvelo;
	}else{ dt = dtvisc;
	}

	return dt;
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
