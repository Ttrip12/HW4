#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

vector<double> DDx(vector<double>, double *, int *,int *,int);
vector<double> DDxDDx(vector<double>, double *, int *,int *,int);
vector<double> fill_gc(vector<double>, int,int);
vector<double> error(vector<double>,vector<double>);
void create(string,vector<double>,vector<double>,int);

int main(){
int n, order, gc;

// *** INPUTS ***
n = 1000; // Number of Cells
order = 2;


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
