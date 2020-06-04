#include "File_string.cpp"
#include "Eigen.cpp"
#include "FCI.cpp"
#include <iomanip>
using namespace std;

int main() 
{
	//输入文件，文件内为几份文件FCIDUMP的名称，方便随时更换文件
	ofstream logfile;
	logfile.open("data\\FCI.log");
	logfile.setf(ios::showpoint); //设置为始终输出小数点后的数字，就是说 a = 3，它也输出 3.00000 这样
	logfile.precision(6);
	logfile.setf(ios::fixed); //设置为小数位始终有 6 位，没有这个的话就会像上面那个代码那样固定的不是小数点后面的数字了。

	ifstream infile;
	string file_name;
	infile.open("data\\FCI.input");
	if(infile.is_open())
        logfile << "Successfully open the input file" << endl;
	else{
        logfile << "fail to open the input file" << endl;
		return 0;
	}
	getline(infile, file_name);
	infile.close();


	file_name = "data\\" + file_name + "\\FCIDUMP";
	int nelec = 0;							//电子数
	int nOrb = 0;							//轨道数
	double MS = 0;							//自旋z分量
	double h_nuc = 0;						//核积分项
	vector<double> h;						//单电子积分
	vector<double> g;						//双电子积分
	//读取积分文件中的积分值
	bool flag = read_int(file_name, nelec, nOrb, MS, h_nuc, h, g);
	if(flag){
		logfile << "successfully read the intergal file" << endl;
	}
	else{
		logfile << "fail to read the intergal file" << endl;
		return 0;
	}


	logfile << "the number of electrons: " << nelec << endl
		 << "the number of active orbitals: " << nOrb << endl
		 << "MS: " << MS << endl;
	cout << "the number of electrons: " << nelec << endl
		 << "the number of active orbitals: " << nOrb << endl
		 << "MS: " << MS << endl;


	//创建FCI类，存储分子轨道积分
	CI FCI(h_nuc, h, g, nOrb);
	//创建组态，CI_Array用于存储组态，Orbital表示暂时表示组态的数组
	vector<Slater_det> CI_Array;
	vector<int> Orbital(2*nOrb);
	flag = CI_new(nelec, 2*nOrb, MS, Orbital, CI_Array);

 	if(flag){
		logfile << "successfully build CIs" << endl;
	}
	else{
		logfile << "fail to build  CIs" << endl;
		return 0;
	}
	int nCI=CI_Array.size();
	logfile << "the number of CIs: " << nCI << endl;
	cout << "the number of CIs: " << nCI << endl;

	for (int i = 0; i < nCI;i++)
		logfile << CI_Array[i];

	//构建Hamilton矩阵
	double temp = 0;
	vector<double> H(nOrb * nOrb);
	for (int i = 0; i < nCI; i++){
		H[i*nOrb+i] = FCI.H_ij(CI_Array[i], CI_Array[i])+h_nuc;
		for (int j = i+1; j < nCI; j++){
			temp = FCI.H_ij(CI_Array[i], CI_Array[j]);
			H[i*nOrb+j] = temp;
			H[j*nOrb+i] = temp;
		}
	}

	CI_Array.clear();
	logfile << "Hamilton matrix";
	logfile << H;

	//对角化Hamilton矩阵
	int nJt = 30;
	vector<double> dbVectors(nCI * nCI);
	vector<double> dbEigenvalues(nCI);
	flag=eigenh(H, dbVectors, dbEigenvalues, nJt);
	if(flag){
		logfile << "Successfully diagonalize the Hamilton matrix" << endl;
	}
	else{
		logfile << "fail to diagonalize the Hamilton matrix" << endl;
		return 0;
	}



	ofstream outfile;
	outfile.open("data\\FCI.output");
	outfile.setf(ios::showpoint); //设置为始终输出小数点后的数字，就是说 a = 3，它也输出 3.00000 这样
	outfile.precision(6);
	outfile.setf(ios::fixed); //设置为小数位始终有 6 位，没有这个的话就会像上面那个代码那样固定的不是小数点后面的数字了。

	if(outfile.is_open())
        logfile << "Successfully open the output file" << endl;
	else{
        logfile << "fail to open the output file" << endl;
		return 0;
	}
	
	outfile << "Eigenvalues:";
	for (int i = 0; i < nCI;i++)
		outfile << dbEigenvalues[i] << "\t";

	outfile << endl
			<< "Eigenvectors: ";
	outfile << dbVectors;
	outfile << endl;

	outfile.close();
	logfile.close();

}