#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<math.h>

using namespace std;

int spot;	//點個數
string *place;		//地名
double* closetPair(double* P, double* Q, int n, int start);		//找到最短路徑

int main() {
	cout << "總共有幾個地點:";
	cin >> spot;
	fstream inputStream;
	double **coordinate = new double*[2];	//二為座標陣列
	place = new string[spot];
	for (int i = 0; i < spot; ++i)
		coordinate[i] = new double[spot];
	inputStream.open("map.txt");
	for (int i = 0; i < spot; ++i) {	//讀地名及座標
		inputStream >> place[i];	
		inputStream >> coordinate[1][i] >> coordinate[0][i];
		cout << place[i] << '(' << fixed << setprecision(6) << coordinate[0][i] << ',' << coordinate[1][i] << ')' << endl;
	}
	for (int i = spot - 1; i >0; --i)		//照x座標排序，對應的y座標和地名順序與x座標排序相同
		for (int j = 0; j < i; ++j) {
			if (coordinate[0][j] > coordinate[0][j + 1]) {
				double temp = 0.0;
				temp = coordinate[0][j];
				coordinate[0][j] = coordinate[0][j + 1];
				coordinate[0][j + 1] = temp;
				temp = coordinate[1][j];
				coordinate[1][j] = coordinate[1][j + 1];
				coordinate[1][j + 1] = temp;
				string change = "";
				change = place[j];
				place[j] = place[j + 1];
				place[j + 1] = change;
			}
		}

	double *minD = new double[3];	//最短距離及兩點的索引值
	minD = closetPair(coordinate[0], coordinate[1], spot, 0);	//遞迴呼叫求出最短距離與兩點
	cout << "--------------------------------------第一題-----------------------------------------------------------------------" << endl;
	cout << "哪兩點最靠近:" << place[((int)minD[1])] << "--" << place[((int)minD[2])] << endl;
	cout << "最短距離:" << sqrt(minD[0]) << endl;
	
	system("pause");
	return 0;
}

double* closetPair(double* P, double* Q, int n, int start) {
	if (n <= 3) {		//divide到剩下三個點以下則用暴力法求最短距離及兩點並回傳
		double *dmin = new double[3];
		dmin[0] = pow((P[0] - P[1]), 2) + pow((Q[0] - Q[1]), 2);	//只有兩點則回傳其距離
		dmin[1] = start;
		dmin[2] = start + 1;
		if (n == 3) {		//若有三點則比較哪兩點為最短距離
			if ((pow((P[0] - P[2]), 2) + pow((Q[0] - Q[2]), 2))<dmin[0]) {
				dmin[0] = pow((P[0] - P[2]), 2) + pow((Q[0] - Q[2]), 2);
				dmin[1] = start;
				dmin[2] = start + 2;
			}

			if ((pow((P[1] - P[2]), 2) + pow((Q[1] - Q[2]), 2))<dmin[0]) {
				dmin[0] = pow((P[1] - P[2]), 2) + pow((Q[1] - Q[2]), 2);
				dmin[1] = start + 1;
				dmin[2] = start + 2;
			}
		}
		return dmin;
	}
	else {		
		int o = ceil((double)n / 2);
		double *Pl = new double[o];		//左半邊陣列(x)
		double *Pr = new double[o];		//右半邊陣列(x)
		double *Ql = new double[o];		//左半邊陣列(y)
		double *Qr = new double[o];		//右半邊陣列(y)
		double *S = new double[n];		//在2d距離內的y座標
		double *R = new double[n];		//在2d距離內的x座標
		double *C = new double[n];		//在2d距離內的地名
		for (int i = 0; i <ceil((double)n / 2); ++i) {		//將陣列分成兩半
			Pl[i] = P[i];
			Ql[i] = Q[i];
		}
		int a = 0;
		for (int i = ceil((double)n / 2); i < n; ++i) {		//將陣列分成兩半
			Pr[a] = P[i];
			Qr[a] = Q[i];
			++a;
		}

		double *dl = new double[3];
		dl = closetPair(Pl, Ql, ceil((double)n / 2), start);		//遞迴divide(左半部)
		double *dr = new double[3];
		dr = closetPair(Pr, Qr, (n / 2), (start + ceil((double)n / 2)));		//遞迴divide(右半部)
		double *d = new double[3];		
		if (dl[0] < dr[0]) {		//取dl和dr的最小值
			d[0] = dl[0];
			d[1] = dl[1];
			d[2] = dl[2];
		}
		else {
			d[0] = dr[0];
			d[1] = dr[1];
			d[2] = dr[2];
		}
		double m = P[o - 1];	//merge時的中線
		int j = 0;	//有j個元素在2d的範圍內
		for (int i = 0; i < n; ++i)		//保存x和中線在2d距離內的點
			if (abs(P[i] - m) < sqrt(d[0])) {
				S[j] = Q[i];
				R[j] = P[i];
				C[j] = start + i;
				++j;
			}

		for (int i = 0; i < j - 1; ++i) {		//判斷中線兩側是否有比d距離更近的點
			int k = i + 1;
			while (k <= j - 1) {
				if ((pow((S[k] - S[i]), 2) < d[0]))
					if ((pow((S[k] - S[i]), 2) + pow((R[k] - R[i]), 2)) < d[0]) {
						d[0] = pow((S[k] - S[i]), 2) + pow((R[k] - R[i]), 2);
						d[1] = C[i];
						d[2] = C[k];
					}
				k++;
			}
		}
		return d;
	}
}
