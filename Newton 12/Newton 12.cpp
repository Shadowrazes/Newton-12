#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include<iomanip>
#include <cmath>

using namespace std;

const double X = 2.25;

struct Eq {
    double x = 0;
    double y = 0;
};

void printVector(vector <Eq>& vect) {
    cout << "X: ";
    for (const auto& it : vect)
        cout << it.x << " ";
    cout << endl;
    cout << "Y: ";
    for (const auto& it : vect)
        cout << it.y << " ";
    cout << endl << endl;
}

void strokeWrap(string& strokeX, string& strokeY, vector <Eq>& equations) {
    while (strokeX.find(" ") != -1) {
        Eq buff = { stod(strokeX.substr(0, strokeX.find(" "))), stod(strokeY.substr(0, strokeY.find(" "))) };
        strokeX.erase(0, strokeX.find(" ") + 1);
        strokeY.erase(0, strokeY.find(" ") + 1);
        equations.push_back(buff);
    }
}

double fact(double x) {
    if (x == 0) return 1;
    return x * fact(x - 1);
}

double qCoef(double q, int ind) {
    if (ind == 1) return q;
    return q * qCoef(q - 1, --ind);
}

double qCoefRev(double q, int ind) {
    if (ind == 1) return q;
    return q * qCoefRev(q + 1, --ind);
}

void initDeltaMatrix(vector <vector <double>> &deltaY, double &h, int &tableLength) {
    for (int i = 0; i < tableLength - 1; i++) {
        vector <double> tmp;
        cout << "DeltaY_" << i + 1 << ": ";
        for (int j = 0; j < deltaY[i].size() - 1; j++) {
            tmp.push_back(deltaY[i][j + 1] - deltaY[i][j]);
            cout << deltaY[i][j + 1] - deltaY[i][j] << " ";
        }
        deltaY.push_back(tmp);
        cout << endl;
    }
}

void NewtonFirst(vector <vector <double>>& deltaY, double q) {
    cout << endl << "Newton's formula #1:";
    double func = deltaY[0][0];
    cout << endl << "f(" << X << ") = " << func;
    for (int i = 1; i < deltaY.size(); i++) {
        double tmp = deltaY[i][0] * qCoef(q, i) / fact((double)i);
        func += tmp;
        cout << " + ("  << tmp << ")";
    }
    cout << " = " << func << endl;
}

void NewtonSecond(vector <vector <double>>& deltaY, double qRev) {
    cout << endl << "Newton's formula #2:";
    double func = deltaY[0][deltaY.size() - 1];
    cout << endl << "f(" << X << ") = " << func;
    int j = 1;
    for (int i = deltaY.size() - 2; i >= 0; i--, j++) {
        double tmp = deltaY[j][i] * qCoefRev(qRev, j) / fact((double)j);
        func += tmp;
        cout << " + (" << tmp << ")";
    }
    cout << " = " << func << endl;
}

int main()
{
    fstream matrixData("Polynomial.txt");
    vector <Eq> equations;

    while (!matrixData.eof()) {
        string buffX;
        string buffY;
        getline(matrixData, buffX);
        getline(matrixData, buffY);
        strokeWrap(buffX, buffY, equations);
    }
    double h = 1;
    int tableLength = equations.size();
    vector <vector <double>> deltaY;

    vector <double> tmp;
    for (auto& it : equations)
        tmp.push_back(it.y);
    deltaY.push_back(tmp);
    printVector(equations);
    initDeltaMatrix(deltaY, h, tableLength);
    double q = (X - equations[0].x) / h;
    double qRev = (X - equations[equations.size() - 1].x) / h;
    NewtonFirst(deltaY, q);
    NewtonSecond(deltaY, qRev);
}