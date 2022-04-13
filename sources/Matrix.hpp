/**
 * Matrix Calculator Exercise
 * 
 * Author: Dor Baram
 * Since : 2022-04
 */

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <stdexcept>
using namespace std;

namespace zich
{
    class Matrix
    {
    public:
        int row;
        int col;
        map<int, map<int, int>> matrix;

        Matrix(vector<double> mat, int inrow, int incol);

        /*input validators*/
        friend void validator(const Matrix &m1, const Matrix &m2);      //check if the two matrix are same nXm size

        /*Input Output*/
        friend ostream &operator<<(ostream &out, const Matrix &m);
        friend istream &operator>>(istream &in, Matrix &m);

        /*Arithmetic*/
        friend Matrix operator+(const Matrix &m1, const Matrix &m2);    //addetive
        friend Matrix operator+=(Matrix &m1, const Matrix &m2);         //add equal
        friend Matrix operator+(const Matrix &m);                       //unari add
        friend Matrix operator-(const Matrix &m1, const Matrix &m2);    //subtract
        friend Matrix operator-=(Matrix &m1, const Matrix &m2);         //sub equal
        friend Matrix operator-(const Matrix &m);                       //unari sub

        /*Compare*/
        friend bool operator>(const Matrix &m1, const Matrix &m2);      //grater
        friend bool operator>=(const Matrix &m1, const Matrix &m2);     //grater equal
        friend bool operator<(const Matrix &m1, const Matrix &m2);      //lesser
        friend bool operator<=(const Matrix &m1, const Matrix &m2);     //lesser equal
        friend bool operator==(const Matrix &m1, const Matrix &m2);     //equal
        friend bool operator!=(const Matrix &m1, const Matrix &m2);     //not equal
    
        /*Increment Decrement*/
        friend Matrix operator++(Matrix &m,int);            //increments m++
        friend Matrix operator++(Matrix &m);                //increments ++m
        friend Matrix operator--(Matrix &m, int);           //decrements m--
        friend Matrix operator--(Matrix &m);                //decrements --m

        /*Mult*/
        friend Matrix operator*(double scalar, Matrix& m);      //scalar mult
        friend Matrix operator*(Matrix& m ,double scalar);      //scalar mult
        friend Matrix operator*=(Matrix& matrix, double scalar);//scalar mult equal
        friend Matrix operator*=(double scalar, Matrix& m);//scalar mult equal
        friend Matrix operator*(Matrix& m1, Matrix& m2);        //matrix mult
        friend Matrix operator*=(Matrix& m1, Matrix& m2);       //matrix mult equal
    };
}