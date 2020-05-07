#include "Array.cpp"
#include "File_string.cpp"
#include "Eigen.cpp"
#include "FCI.cpp"
using namespace std;

int main() 
{
	//输入文件，文件内为几份文件FCIDUMP的名称，方便随时更换文件
	string inputfile = "data\\FCI.input";
	ifstream infile;
	ofstream logfile;
	string file_name;

	logfile.open("data\\FCI.log");
	infile.open(inputfile.c_str());
	
	if(infile.is_open())
        logfile << "成功打开输入文件" << endl;
	else{
        logfile << "未成功打开输入文件" << endl;
		return 0;
	}
	getline(infile, file_name);
	infile.close();
	file_name = "data\\" + file_name+"\\FCIDUMP";

	
	double h_nuc = 0;
	int nelec = 0;
	int nOrb = 0;

    ifstream file;
    file.open(file_name.c_str());
    string line;

    //第一行格式如右，用于读取轨道个数与电子个数:"&FCI NORB=114,NELEC=42,MS2=0,"
    getline(file, line);
    string delim = "=";
    vector<string> v = split(line, delim);

    delim = ",";
    vector<string> v1 = split(v[1], delim);
    vector<string> v2 = split(v[2], delim);
    nelec = atof(v2[0].c_str());
    nOrb = atof(v1[0].c_str());

	double **h = double2_new(nOrb, nOrb);
	double ****g = double4_new(nOrb, nOrb, nOrb, nOrb);

	
	
	//读取积分文件中的积分值
 	bool flag=read_int(file_name, nelec, nOrb, h_nuc, h, g);

 	if(flag){
		logfile << "积分文件读取成功" << endl;
	}
	else{
		logfile << "积分文件读取失败" << endl;
		return 0;
	}
	logfile << "电子数: " << nelec << endl
		 << "活性轨道数: " << nOrb << endl;
	cout << "电子数: " << nelec << endl
		 << "活性轨道数: " << nOrb << endl;


	logfile << "单电子积分";
	output(logfile, h, nOrb);
	logfile << "双电子积分";
	output(logfile, g, nOrb);


	//创建FCI类，存储分子轨道积分
	CI FCI(h_nuc, h, g, nOrb);
	//创建组态，CI_Array用于存储组态，Orbital表示暂时表示组态的数组
	vector<Slater_det> CI_Array;
	int **Orbital = int2_new(nOrb, 2);
	flag = CI_new(nelec, nOrb, nelec, 2*nOrb, Orbital, CI_Array);

 	if(flag){
		logfile << "组态构建成功" << endl;
	}
	else{
		logfile << "组态构建失败" << endl;
		return 0;
	}
	int nCI=CI_Array.size();
	logfile << "组态个数：" << nCI << endl;
	cout << "组态个数：" << nCI << endl;

	for (int i = 0; i < nCI;i++)
		logfile << CI_Array[i];


	//构建Hamilton矩阵
	int nJt = 30;
	double **H = double2_new(nCI, nCI);
	for (int i = 0; i < nCI; i++){
		H[i][i] = FCI.H_ij(CI_Array[i], CI_Array[i])+h_nuc;
		for (int j = i+1; j < nCI; j++){
			H[i][j] = FCI.H_ij(CI_Array[i], CI_Array[j]);
			H[j][i] = FCI.H_ij(CI_Array[i], CI_Array[j]);
		}
	}
	CI_Array.clear();
	
	

	//对角化Hamilton矩阵
	double **dbVectors=double2_new(nCI, nCI);
	double dbEigenvalues[nCI]={0};
	flag=eigenh(H, dbVectors, dbEigenvalues, nCI, nJt);
	if(flag){
		logfile << "Hamilton矩阵对角化成功" << endl;
	}
	else{
		logfile << "Hamilton矩阵对角化失败" << endl;
		return 0;
	}

	double2_delete(h, nOrb);
	double2_delete(H, nOrb);
	double4_delete(g, nOrb, nOrb, nOrb);


	ofstream outfile;
	outfile.open("file\\FCI.output");

	if(outfile.is_open())
        logfile << "成功打开输出文件" << endl;
	else{
        logfile << "未成功打开输出文件" << endl;
		return 0;
	}
	
	outfile << "本征能量：" << endl;
	for (int i = 0; i < nCI; i++){
		outfile << dbVectors[i] << "\t";
	}
	
	outfile << "特征向量：" << endl;
	for (int i = 0; i < nCI; i++){
		for (int j = 0; j < nCI; j++){
			outfile << dbVectors[i][j] << "\t";
		}
		outfile << endl;
	}

	outfile.close();
	logfile.close();
	double2_delete(dbVectors, nCI);

	
}