#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<math.h>

using namespace std;

int spot;	//�I�Ӽ�
string *place;		//�a�W
double* closetPair(double* P, double* Q, int n, int start);		//���̵u���|

int main() {
	cout << "�`�@���X�Ӧa�I:";
	cin >> spot;
	fstream inputStream;
	double **coordinate = new double*[2];	//�G���y�а}�C
	place = new string[spot];
	for (int i = 0; i < spot; ++i)
		coordinate[i] = new double[spot];
	inputStream.open("map.txt");
	for (int i = 0; i < spot; ++i) {	//Ū�a�W�ήy��
		inputStream >> place[i];	
		inputStream >> coordinate[1][i] >> coordinate[0][i];
		cout << place[i] << '(' << fixed << setprecision(6) << coordinate[0][i] << ',' << coordinate[1][i] << ')' << endl;
	}
	for (int i = spot - 1; i >0; --i)		//��x�y�бƧǡA������y�y�ЩM�a�W���ǻPx�y�бƧǬۦP
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

	double *minD = new double[3];	//�̵u�Z���Ψ��I�����ޭ�
	minD = closetPair(coordinate[0], coordinate[1], spot, 0);	//���j�I�s�D�X�̵u�Z���P���I
	cout << "--------------------------------------�Ĥ@�D-----------------------------------------------------------------------" << endl;
	cout << "�����I�̾a��:" << place[((int)minD[1])] << "--" << place[((int)minD[2])] << endl;
	cout << "�̵u�Z��:" << sqrt(minD[0]) << endl;
	
	system("pause");
	return 0;
}

double* closetPair(double* P, double* Q, int n, int start) {
	if (n <= 3) {		//divide��ѤU�T���I�H�U�h�μɤO�k�D�̵u�Z���Ψ��I�æ^��
		double *dmin = new double[3];
		dmin[0] = pow((P[0] - P[1]), 2) + pow((Q[0] - Q[1]), 2);	//�u�����I�h�^�Ǩ�Z��
		dmin[1] = start;
		dmin[2] = start + 1;
		if (n == 3) {		//�Y���T�I�h��������I���̵u�Z��
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
		double *Pl = new double[o];		//���b��}�C(x)
		double *Pr = new double[o];		//�k�b��}�C(x)
		double *Ql = new double[o];		//���b��}�C(y)
		double *Qr = new double[o];		//�k�b��}�C(y)
		double *S = new double[n];		//�b2d�Z������y�y��
		double *R = new double[n];		//�b2d�Z������x�y��
		double *C = new double[n];		//�b2d�Z�������a�W
		for (int i = 0; i <ceil((double)n / 2); ++i) {		//�N�}�C������b
			Pl[i] = P[i];
			Ql[i] = Q[i];
		}
		int a = 0;
		for (int i = ceil((double)n / 2); i < n; ++i) {		//�N�}�C������b
			Pr[a] = P[i];
			Qr[a] = Q[i];
			++a;
		}

		double *dl = new double[3];
		dl = closetPair(Pl, Ql, ceil((double)n / 2), start);		//���jdivide(���b��)
		double *dr = new double[3];
		dr = closetPair(Pr, Qr, (n / 2), (start + ceil((double)n / 2)));		//���jdivide(�k�b��)
		double *d = new double[3];		
		if (dl[0] < dr[0]) {		//��dl�Mdr���̤p��
			d[0] = dl[0];
			d[1] = dl[1];
			d[2] = dl[2];
		}
		else {
			d[0] = dr[0];
			d[1] = dr[1];
			d[2] = dr[2];
		}
		double m = P[o - 1];	//merge�ɪ����u
		int j = 0;	//��j�Ӥ����b2d���d��
		for (int i = 0; i < n; ++i)		//�O�sx�M���u�b2d�Z�������I
			if (abs(P[i] - m) < sqrt(d[0])) {
				S[j] = Q[i];
				R[j] = P[i];
				C[j] = start + i;
				++j;
			}

		for (int i = 0; i < j - 1; ++i) {		//�P�_���u�ⰼ�O�_����d�Z������I
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
