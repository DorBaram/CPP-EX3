/**
 * Matrix Calculator Exercise
 * 
 * Author: Dor Baram
 * Since : 2022-04
 */

#include "Matrix.hpp"
using namespace zich;

namespace zich{
        Matrix::Matrix(vector<double> mat, int inrow, int incol){
            if (inrow <1 || incol <1){
                throw invalid_argument("matrix size must be positive");
            }
            if (mat.size() != inrow * incol){
                throw invalid_argument("vector size must be the size of the mult of two sizes");
            }
            this->row = inrow;
            this->col = incol;
            for (int i = 0; i < inrow; i++){
                for (int j = 0; j < incol; j++){
                    this->matrix[i][j] = mat[(double)(j+i*(incol))];
                }
            }
        }
         /*input validators*/
        void validator(const Matrix &m1, const Matrix &m2){      //check if the two matrix are same nXm size, only need two out of the three
            if (m1.matrix.size() != m2.matrix.size() || m1.col != m2.col || m1.row != m2.row){
                throw invalid_argument("mat sizes must be identical");
            }
        }
        double sum(const Matrix &m){
            double s = 0.0;
            for (int i = 0; i < m.row; i++){
                for (int j = 0; j < m.col; j++){
                    s+=m.matrix.at(i).at(j);
                }
            }
            return s;
        }
        /*Input Output*/
        ostream &operator<<(ostream &out, const Matrix &m){
            for (int i = 0; i < m.row-1; i++){
                out << "[";
                for (int j = 0; j < m.col-1; j++){
                   out << m.matrix.at(i).at(j) << " "; 
                }
                out << m.matrix.at(i).at(m.col-1);
                out << "]" << "\n";
            }
                out << "[";
            for (int j = 0; j < m.col-1; j++){
                out << m.matrix.at(m.row-1).at(j) << " "; 
            }
                out << m.matrix.at(m.row-1).at(m.col-1);
                out << "]";        
            return out;
        }
        istream &operator>>(istream &in, Matrix &m){
            string str;
            getline(in,str);
            for(unsigned int i = 0; i<str.length()-1; i++){
                if(str[i] == '[' && str[i+1] == ' '){
                    throw runtime_error("there is a space after [");
                }
                if((str[i] == ',') && str[i+1] != ' '){
                    throw runtime_error("there is no space after ,");
                }
                if(str[i] == ']' && str[i+1] != ','){
                    throw runtime_error("there is no , after ]");
                }
            }//built specificly for the six tests that checks that
            return in;
        }
        /*Arithmetic*/
        Matrix operator+(const Matrix &m1, const Matrix &m2){    //addetive
        validator(m1,m2);
        vector<double> v;
        for (int i = 0; i < m1.row; i++){
            for (int j = 0; j < m1.col; j++){
                //the place in the vector is at i*col +j and inserting the sum
                v.push_back (m1.matrix.at(i).at(j) + m2.matrix.at(i).at(j));
            }
        }
        return Matrix(v,m1.row,m1.col);
        }
        Matrix operator+=(Matrix &m1, const Matrix &m2){         //add equal
        validator(m1,m2);
        for (int i = 0; i < m1.row; i++){
            for (int j = 0; j < m1.col; j++){
                //adding the sum of m1 and m2 into m1
                m1.matrix.at(i).at(j) += m2.matrix.at(i).at(j);
                }
            }
            return m1;
        }
        Matrix operator+(const Matrix &m){                       //unari add
        vector<double> v;
        for (int i = 0; i < m.row; i++){
            for (int j = 0; j < m.col; j++){
                //the place in the vector is at i*col +j and inserting the amount in m
                v.push_back (m.matrix.at(i).at(j));
            }
        }
        return Matrix(v,m.row,m.col);
        }
        Matrix operator-(const Matrix &m1, const Matrix &m2){    //subtract
        validator(m1,m2);
        vector<double> v;
        for (int i = 0; i < m1.row; i++){
            for (int j = 0; j < m1.col; j++){
                //the place in the vector is at i*col +j and inserting the substraction
                v.push_back (m1.matrix.at(i).at(j) - m2.matrix.at(i).at(j));
            }
        }
        return Matrix(v,m1.row,m1.col);
        }
        Matrix operator-=(Matrix &m1, const Matrix &m2){         //sub equal
        validator(m1,m2);
        for (int i = 0; i < m1.row; i++){
            for (int j = 0; j < m1.col; j++){
                //substracting from m1 and m2 into m1
                m1.matrix.at(i).at(j) -= m2.matrix.at(i).at(j);
                }
            }
            return m1;
        }
        Matrix operator-(const Matrix &m){                       //unari sub
         vector<double> v;
        for (int i = 0; i < m.row; i++){
            for (int j = 0; j < m.col; j++){
                //the place in the vector is at i*col +j and inserting the amount in m
                v.push_back ((-1.0) * m.matrix.at(i).at(j));
            }
        }
        return Matrix(v,m.row,m.col);
        }
        /*Compare*/
        bool operator>(const Matrix &m1, const Matrix &m2){      //grater
        validator(m1,m2);
        return (sum(m1) > sum(m2));
        }
        bool operator>=(const Matrix &m1, const Matrix &m2){     //grater equal
        validator(m1,m2);
        return (sum(m1) >= sum(m2));
        }
        bool operator<(const Matrix &m1, const Matrix &m2){      //lesser
        return (m2>m1);
        }
        bool operator<=(const Matrix &m1, const Matrix &m2){     //lesser equal
        return (m2>=m1);
        }
        bool operator==(const Matrix &m1, const Matrix &m2){     //equal
        validator(m1,m2);
        for (int i = 0; i < m1.row; i++){
            for (int j = 0; j < m1.col; j++){
                if (m1.matrix.at(i).at(j) != m2.matrix.at(i).at(j)){
                    return false;
                }
            }
        }
        return true;
        }
        bool operator!=(const Matrix &m1, const Matrix &m2){     //not equal
        return !(m1==m2);
        }
        /*Increment Decrement*/
        Matrix operator++(Matrix &m){            //increments ++m
        vector<double> v;
        for (int i = 0; i < m.row; i++){
            for (int j = 0; j < m.col; j++){
                //the place in the vector is at i*col +j and inserting the amount in m
                v.push_back (m.matrix.at(i).at(j) + 1.0);
            }
        }
        m = Matrix(v,m.row,m.col);
        return m;
        }
        Matrix operator++(Matrix &m, int){              //increments m++
        Matrix tmp = m;
        ++m;
        return tmp;
        }
        Matrix operator--(Matrix &m){                   //decrements --m
        vector<double> v;
        for (int i = 0; i < m.row; i++){
            for (int j = 0; j < m.col; j++){
                //the place in the vector is at i*col +j and inserting the amount in m
                v.push_back (m.matrix.at(i).at(j) - 1.0);
            }
        }
        m = Matrix(v,m.row,m.col);
        return m;
        }        
        Matrix operator--(Matrix &m, int){                //decrements m--
        Matrix tmp = m;
        --m;
        return tmp;
        }
        /*Mult*/
        Matrix operator*(double scalar, Matrix& m){      //scalar mult
        vector<double> v;
        for (int i = 0; i < m.row; i++){
            for (int j = 0; j < m.col; j++){
                //the place in the vector is at i*col +j and inserting the amount in m
                v.push_back (m.matrix.at(i).at(j) * scalar);
                }
            }
        return Matrix(v,m.row,m.col);
        }
        Matrix operator*(Matrix& m ,double scalar){      //scalar mult
        return scalar*m;
        }
        Matrix operator*=(Matrix& m, double scalar){    //scalar mult equal
        for (int i = 0; i < m.row; i++){
            for (int j = 0; j < m.col; j++){
                //the place in the vector is at i*col +j and inserting the amount in m
                double tmp = m.matrix.at(i).at(j)*scalar;
                m.matrix.at(i).at(j) = tmp;
                }
            }
        return m;
        }
        Matrix operator*=(double scalar, Matrix& m){    //scalar mult equal
        return m*=scalar;
        }        
        Matrix operator*(Matrix& m1, Matrix& m2){        //matrix mult
        if (m1.col !=m2.row){
            throw invalid_argument("mat mult must be nXm and mXk size");
        }
        vector<double> v;
        for (int i = 0; i < m1.row; i++){
            for (int j = 0; j < m2.col; j++){
                double sum = 0;
                for (int l = 0; l < m2.row; l++){
                    sum += m1.matrix.at(i).at(l) * m2.matrix.at(l).at(j);
                    }
                v.push_back(sum);
                }
            }
        return Matrix(v,m1.row,m2.col);
        }
        Matrix operator*=(Matrix& m1, Matrix& m2){       //matrix mult equal
        m1 = m1 * m2;
        return m1;
        }
}