#include "Array.cpp"
#include "File_string.cpp"
#include "Eigen.cpp"
#include "FCI.cpp"
using namespace std;

int main() 
{
	string inputfile_name = "input\\FCI.input";
	ifstream file;
	string file_name;

	file.open(inputfile_name.c_str());
	getline(file, file_name);
	file.close();
	file_name = "input\\" + file_name;

	
	double **h = NULL;
	double ****g = NULL;
	double h_nuc = 0;
	int nelec = 0;
	int nOrb = 0;
	 
 	bool flag=read_int(file_name, nelec, nOrb, h_nuc, h, g);

 	if(!flag){
		cout << "NO";
		return 0;
	}

	
	CI FCI(h_nuc, h, g, nOrb);
	cout << nelec << " " << nOrb << endl;

	vector<Slater_det> CI_Array;
	int **Orbital = int2_new(nOrb, 2);
	CI_new(nelec, nOrb, nelec, 2*nOrb, Orbital, CI_Array);

	int nCI=CI_Array.size();
	int nJt = 30;
	double **H = double2_new(nCI, nCI);
	for (long long i = 0; i < nCI; i++){
		H[i][i] = FCI.H_ij(CI_Array[i], CI_Array[i])+h_nuc;
		for (long long j = i+1; j < nCI; j++){
			H[i][j] = FCI.H_ij(CI_Array[i], CI_Array[j]);
			H[j][i] = FCI.H_ij(CI_Array[i], CI_Array[j]);
		}
	}
	CI_Array.clear();

	double **dbVectors=double2_new(nCI, nCI);
	double dbEigenvalues[nCI]={0};
	eigenh(H, dbVectors, dbEigenvalues, nCI, nJt);

	double2_delete(h, nOrb);
	double2_delete(H, nOrb);
	double4_delete(g, nOrb, nOrb, nOrb);

	
	double2_delete(dbVectors, nCI);

	
}