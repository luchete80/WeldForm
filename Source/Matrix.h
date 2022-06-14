#ifndef _MATRIX_H_
#define _MATRIX_H_

namespace SPH{

class Matrix{
public:
	Matrix (){} // if use a NON pointer matrix in solver we need this
	Matrix (int row, int col){
    row_count = row;
    col_count = col;
    //initialize 10 x 20 array:
    m_data = new double*[row];
    for (int i = 0; i < row; i++)
      m_data[i] = new double[col];
    for (int i=0;i<row;i++){
      for (int j=0;j<col;j++)
        m_data[i][j] = 0.;
    }
	}
	double **m_data;	
	
	double ** Get(){return m_data;}
  inline void Set(const int &i, const int &j, const double &val){m_data[i][j] = val;}
  inline void Add(const int &i, const int &j, const double &val){m_data[i][j] += val;}
	inline double operator()(int i, int j){return m_data[i][j];}

	
	inline Matrix Transpose();
	inline Matrix operator*(const Matrix Right) ;
	inline Matrix operator+(const Matrix Right) ;
  
private:
  
	int row_count,col_count;
};


inline Matrix Matrix::Transpose(){
	Matrix ret(col_count, row_count);
	for (int i=0;i<row_count;i++)
		for (int j=0;j<col_count;j++)		
			ret.m_data[j][i]=m_data[i][j];
	
	return ret;
}

inline Matrix Matrix::operator*(const Matrix Right)  {
	Matrix Ret (row_count, Right.col_count);
	
	if (col_count !=Right.row_count){
		//cout << "bad matrix operation *"<<endl;
		return Ret;
	}
	for (unsigned int i = 0 ; i < row_count ; i++) {
			for (unsigned int j = 0 ; j < Right.col_count ; j++) {
				for (unsigned int k = 0 ; k < col_count ; k++)
					Ret.m_data[i][j] += m_data[i][k]*Right.m_data[k][j];
			}
	}

	return Ret;
}

inline Matrix Matrix::operator+(const Matrix Right)  {
	Matrix Ret (row_count, col_count);

	for (unsigned int i = 0 ; i < row_count ; i++) {
			for (unsigned int j = 0 ; j < col_count ; j++) {
					Ret.m_data[i][j] += m_data[i][j]+Right.m_data[i][j];
			}
	}

	return Ret;
}

inline Matrix Identity(const int &c){
  Matrix ret(c,c);
  for (int i=0;i<c;i++)
    ret.m_data[i][i]=1.;
  return ret;
}

};

#endif