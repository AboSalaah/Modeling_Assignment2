#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <map>
#include <set>
#include <bitset>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <cmath>
#include<cstring>
#include <cstdio>
#include <stack>
#include<iomanip>
#include<queue>
#include <unordered_set>
#include<functional>
#include<iterator>
#include<new>
#include<cstdlib>
#include<math.h>
#include<fstream>
#include"Matrix.h"
using namespace std;
struct in {
	int sourceNode, distNode;
	double val, initialVal;
	in(int src,int dis,double v,double i):sourceNode(src),distNode(dis),val(v),initialVal(i){}
};
vector<in>Resistors;
vector<in>VoltageSources;
vector<in>CurrentSources;
vector<in>Capacitors;
vector<in>Inductors;
vector<vector<double>>output_data_V;
vector<vector<double>>output_data_I;
int nodecnt;
const double TIMESTEP = 0.1;
const int NUMBEROFITERATIONS = 20;
void readFile() {
	ifstream input("input1.txt");
	string compType; int src, dis; double v, init;
	while (input >> compType >> src >> dis >> v >> init) {
		if(compType[0]=='R')
		Resistors.push_back(in(src, dis, v, init));
		else if(compType.size()==1&&compType[0]=='I')
			Inductors.push_back(in(src, dis, v, init));
		else if(compType.size()>1&&compType[0]=='I')
			CurrentSources.push_back(in(src, dis, v, init));
		else if(compType[0]=='C')
			Capacitors.push_back(in(src, dis, v, init));
		else if(compType[0]=='V')
			VoltageSources.push_back(in(src, dis, v, init));
		nodecnt = max(nodecnt, max(src, dis));
	}
}
Matrix buildMatrixG() {
	Matrix G(nodecnt + 1, nodecnt + 1);
	for (int i = 0; i < Resistors.size(); ++i) {
		G.add(Resistors[i].sourceNode,Resistors[i].sourceNode , (1.0 / Resistors[i].val));
		G.add(Resistors[i].sourceNode,Resistors[i].distNode, -1*(1.0 / Resistors[i].val));
		G.add(Resistors[i].distNode,Resistors[i].sourceNode, -1*(1.0 / Resistors[i].val));
		G.add(Resistors[i].distNode,Resistors[i].distNode ,(1.0 / Resistors[i].val));
	}
	for (int i = 0; i < Capacitors.size(); ++i) {
		G.add(Capacitors[i].sourceNode,Capacitors[i].sourceNode,  (1.0 / (TIMESTEP/Capacitors[i].val)));
		G.add(Capacitors[i].sourceNode,Capacitors[i].distNode, -1*(1.0 / (TIMESTEP/Capacitors[i].val)));
		G.add(Capacitors[i].distNode,Capacitors[i].sourceNode, -1* (1.0 / (TIMESTEP / Capacitors[i].val)));
		G.add(Capacitors[i].distNode,Capacitors[i].distNode ,(1.0 /(TIMESTEP / Capacitors[i].val)));
	}
	//G.printMatrix();
	return G;
}
Matrix buildMatrixB() {
	Matrix B(nodecnt + 1, (int)VoltageSources.size() + (int)Inductors.size() + 1);
	for (int i = 0; i < VoltageSources.size(); ++i) {
		B.set(VoltageSources[i].sourceNode,i + 1, 1);
		B.set(VoltageSources[i].distNode,i + 1 , -1);
	}
	for (int i = 0; i < Inductors.size(); ++i) {
		B.set(Inductors[i].sourceNode,(int)VoltageSources.size() + i + 1 , 1);
		B.set(Inductors[i].distNode,(int)VoltageSources.size() + i + 1,  -1);
	}
	//B.printMatrix();
	return B;
}
Matrix buildMatrixC(Matrix B) {
	Matrix C((int)VoltageSources.size() + (int)Inductors.size() + 1, nodecnt + 1);
	for (int i = 0; i < C.getrows(); ++i) {
		for (int j = 0; j < C.getcols(); ++j) {
			C.set(i,j,B.get(j,i));
		}
	}
	//C.printMatrix();
	return C;
}
Matrix buildMatrixD() {
	Matrix D((int)VoltageSources.size() + (int)Inductors.size() + 1, (int)VoltageSources.size() + (int)Inductors.size() + 1);
	for (int i = 0; i < Inductors.size(); ++i) {
		D.set((int)VoltageSources.size() + i + 1,(int)VoltageSources.size() + i + 1 , (-1 * Inductors[i].val) / TIMESTEP);
	}
	//D.printMatrix();
	return D;
}
Matrix buildMatrixZ(bool f, Matrix X) {
	Matrix Z(nodecnt + (int)VoltageSources.size() + (int)Inductors.size() + 1, 2);
	
		for (int i = 0; i < CurrentSources.size(); ++i) {
			Z.add(CurrentSources[i].sourceNode, 1, CurrentSources[i].val);
			Z.add(CurrentSources[i].distNode, 1, -1 * CurrentSources[i].val);
		}
	
		for (int i = 0; i < Capacitors.size(); ++i) {
			if (f == 0) {
				Z.add(Capacitors[i].sourceNode, 1, Capacitors[i].initialVal*(Capacitors[i].val / TIMESTEP));
				Z.add(Capacitors[i].distNode, 1, -1 * (Capacitors[i].initialVal*(Capacitors[i].val / TIMESTEP)));
			}
			else {
				Z.add(Capacitors[i].sourceNode, 1, (X.get(Capacitors[i].sourceNode, 0) - X.get(Capacitors[i].distNode, 0))*(Capacitors[i].val / TIMESTEP));
				Z.add(Capacitors[i].distNode, 1, -1 * (X.get(Capacitors[i].sourceNode, 0) - X.get(Capacitors[i].distNode, 0))*(Capacitors[i].val / TIMESTEP));
			}
		}
		
			for (int i = 0; i < VoltageSources.size(); ++i) {
				Z.set(nodecnt + i + 1, 1, VoltageSources[i].val);
			}
	for (int i = 0; i < Inductors.size(); ++i) {
		if(f==0)Z.set(nodecnt + (int)VoltageSources.size() + i + 1, 1, (-1 * (Inductors[i].initialVal)*(Inductors[i].val / TIMESTEP)));
		else Z.set(nodecnt + (int)VoltageSources.size() + i + 1,1, (-1 * (X.get(nodecnt + (int)VoltageSources.size() + i+1,0))*(Inductors[i].val / TIMESTEP)));
	}
	//Z.printMatrix();
	return Z;
}
Matrix buildMatrixA(Matrix G, Matrix C,Matrix B, Matrix D) {
	Matrix A(nodecnt + (int)VoltageSources.size() + (int)Inductors.size()+1, (nodecnt + (int)VoltageSources.size() + (int)Inductors.size()+1));
	
	for (int i = 1; i <= nodecnt; ++i) {
		for (int j = 1; j <= nodecnt; ++j)A.set(i-1,j-1, G.get(i,j));
	}
	for (int i = nodecnt + 1; i < A.getrows(); ++i) {
		for (int j = 1; j <= nodecnt; ++j)A.set(i-1,j-1 ,C.get(i - nodecnt,j));
	}
	for (int i = 1; i <= nodecnt; ++i) {
		for (int j = nodecnt + 1; j < A.getrows(); ++j)A.set(i-1,j-1,  B.get(i,j - nodecnt));
	}
	for (int i = nodecnt + 1; i < A.getrows(); ++i) {
		for (int j = nodecnt + 1; j < A.getrows(); ++j)A.set(i-1,j-1, D.get(i - nodecnt,j - nodecnt));
	}
	//cout << endl;
	//A.printMatrix();
	return A;
}
Matrix getMatrixX(Matrix Inv, Matrix Z) {
	
	Matrix X(Inv.getrows(), 2);
	for (int i = 1; i < Inv.getrows() ; ++i) {
		double sum = 0;
		for (int j = 1; j < Inv.getrows(); ++j) {
			sum += Inv.get(i-1, j-1)*Z.get(j, 1);
		}
		X.set(i, 0, sum);
	}
	
	return X;
}
void ouputAnswer(Matrix A,Matrix B, Matrix C, Matrix D, Matrix G, Matrix Z) {
	ofstream output;
	output.open("output.txt");
	Matrix Inv(A.getrows(), A.getrows());
	Matrix temp1;
	temp1.MatrixInversion(A.getMatrix(), A.getrows()-1 , Inv.getMatrix());
	
	if ((int)Capacitors.size() == 0 && (int)Inductors.size() == 0) {
		
		//temp1.MatrixMultiplication(Inv.getMatrix(), A.getrows(), Z.getMatrix(), X.getMatrix());
		Matrix X = getMatrixX(Inv,Z);
		/*for (int i = 1; i < A.getrows() ; i++) {
			cout << X.get(i, 0) << endl;
		}
		*/
		for (int i = 1; i <=nodecnt; ++i) {
			output << "V" << i << endl;
			output << X.get(i, 0);
			output << endl << endl;
		}
		for (int i = 0; i < ((int)VoltageSources.size() + (int)Inductors.size()); ++i) {
			output << "Ivsrc" << i+1  << endl;
			output << X.get(nodecnt+i+1, 0);
			output << endl << endl;
		}
	}
	else {
		output_data_V.resize(nodecnt+1);
		output_data_I.resize((int)VoltageSources.size() + (int)Inductors.size() + 1);
		for (int i = 0; i < NUMBEROFITERATIONS; ++i) {
			Matrix X = getMatrixX(Inv, Z);
			for (int j = 1; j <=nodecnt; ++j) {
				output_data_V[j].push_back(X.get(j, 0));
			}
			for (int j = nodecnt+1; j < ((int)VoltageSources.size() + (int)Inductors.size() + nodecnt+1); ++j) {
				output_data_I[j - nodecnt].push_back(X.get(j, 0));
			}
			for (int j = 0; j < Z.getrows(); ++j)Z.set(j,1,0.0);
			Z = buildMatrixZ(1, X);
		}
		double cnt = 0.1;
		for (int i = 1; i <= nodecnt; ++i) {
			output << "V" << i<< endl;
			for (int j = 0; j < output_data_V[i].size(); ++j) {
				output << cnt<<" " << setprecision(15) <<output_data_V[i][j];
				output << endl;
				cnt += 0.1;
			}

			output << endl;
		}
		cnt = 0.1;
		for (int i = 0; i < ((int)VoltageSources.size() + (int)Inductors.size()); ++i) {
			output << "I_L" << i << endl;
			for (int j = 0; j < output_data_I[i + 1].size(); ++j) {
				output << cnt<<" "<< setprecision(15)<<output_data_I[i + 1][j];
				output << endl;
				cnt += 0.1;
			}
			output << endl;
		}

	}
	output.close();
}

int main() {
	readFile();
	Matrix G=buildMatrixG();
	Matrix B=buildMatrixB();
	Matrix C=buildMatrixC(B);
	Matrix D=buildMatrixD();
	Matrix temp;
	Matrix Z=buildMatrixZ(0,temp);
	Matrix A=buildMatrixA(G,C,B,D);
	
	ouputAnswer(A, B, C, D, G, Z);
}