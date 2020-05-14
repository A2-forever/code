#include "Array.cpp"
#include "File_string.cpp"
#include "Eigen.cpp"
#include "FCI.cpp"
#include <iomanip>
using namespace std;

int main() 
{
	//输入文件，文件内为几份文件FCIDUMP的名称，方便随时更换文件
	string inputfile = "data\\FCI.input";
	ifstream infile;
	ofstream logfile;
	string file_name;

	logfile.open("data\\FCI.log");
	logfile.setf(ios::showpoint); //设置为始终输出小数点后的数字，就是说 a = 3，它也输出 3.00000 这样
	logfile.precision(6);
	logfile.setf(ios::fixed); //设置为小数位始终有 6 位，没有这个的话就会像上面那个代码那样固定的不是小数点后面的数字了。

	infile.open(inputfile.c_str());
	
	if(infile.is_open())
        logfile << "Successfully open the input file" << endl;
	else{
        logfile << "fail to open the input file" << endl;
		return 0;
	}
	getline(infile, file_name);
	infile.close();
	file_name = "data\\" + file_name + "\\FCIDUMP";

	int nelec = 0;
	int nOrb = 0;
	//读取文件中的轨道数与电子数，由于无法解决的bug：空指针作为参数输入函数后，得到的还是空指针，故采用此法
    ifstream file;
    file.open(file_name.c_str());
    string line;

    //第一行格式如右，用于读取轨道个数与电子个数:"&FCI NORB=114,NELEC=42,MS2=0,"
    getline(file, line);
	file.close();
    string delim = "=";
    vector<string> v = split(line, delim);

    delim = ",";
    vector<string> v1 = split(v[1], delim);
    vector<string> v2 = split(v[2], delim);
    nelec = atof(v2[0].c_str());
    nOrb = atof(v1[0].c_str());

	double h_nuc = 0;
	double **h = double2_new(nOrb, nOrb);
	double ****g = double4_new(nOrb, nOrb, nOrb, nOrb);


	//读取积分文件中的积分值
 	bool flag=read_int(file_name, nelec, nOrb, h_nuc, h, g);

 	if(flag){
		logfile << "successfully read the intergal file" << endl;
	}
	else{
		logfile << "fail to read the intergal file" << endl;
		return 0;
	}
	logfile << "the number of electron: " << nelec << endl
		 << "the number of active orbital: " << nOrb << endl;
	cout << "the number of electron: " << nelec << endl
		 << "the number of active orbital: " << nOrb << endl;

	logfile << "one electron intergal";
	output(logfile, h, nOrb);
	logfile << "two electrons intergal";
	output(logfile, g, nOrb);



	//创建FCI类，存储分子轨道积分
	CI FCI(h_nuc, h, g, nOrb);
	//创建组态，CI_Array用于存储组态，Orbital表示暂时表示组态的数组
	vector<Slater_det> CI_Array;
	int **Orbital = int2_new(nOrb, 2);
	flag = CI_new(nelec, nOrb, nelec, 2*nOrb, Orbital, CI_Array);

 	if(flag){
		logfile << "successfully build CIs" << endl;
	}
	else{
		logfile << "fail to build  CIs" << endl;
		return 0;
	}
	int nCI=CI_Array.size();
	logfile << "the number of CI: " << nCI << endl;
	cout << "the number of CI: " << nCI << endl;

	for (int i = 0; i < nCI;i++)
		logfile << CI_Array[i];

	//构建Hamilton矩阵
	double temp = 0;
	double **H = double2_new(nCI, nCI);
	for (int i = 0; i < nCI; i++){
		H[i][i] = FCI.H_ij(CI_Array[i], CI_Array[i])+h_nuc;
		for (int j = i+1; j < nCI; j++){
			temp = FCI.H_ij(CI_Array[i], CI_Array[j]);
			H[i][j] = temp;
			H[j][i] = temp;
		}
	}

	CI_Array.clear();
	logfile << "Hamilton matrix";
	output(logfile, H, nCI);
	
	
	//对角化Hamilton矩阵
	int nJt = 3;
	double **dbVectors=ones(nCI); 
	double dbEigenvalues[nCI]={0};
	flag=eigenh(H, dbVectors, dbEigenvalues, nCI, nJt);
	if(flag){
		logfile << "Successfully diagonalize the Hamilton matrix" << endl;
	}
	else{
		logfile << "fail to diagonalize the Hamilton matrix" << endl;
		return 0;
	}

	double2_delete(h, nOrb);
	double2_delete(H, nCI);
	double4_delete(g, nOrb, nOrb, nOrb);


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
	output(outfile, dbEigenvalues, nCI);

	outfile << endl
			<< "Eigenvectors: ";
	output(outfile, dbVectors, nCI);
	outfile << endl;

	outfile.close();
	logfile.close();
	double2_delete(dbVectors, nCI);

}