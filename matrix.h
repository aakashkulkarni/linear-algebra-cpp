#ifndef MATRIX2_H
#define MATRIX2_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <exception>

template <class T>
class Matrix2
{
public:
    Matrix2();
    Matrix2(int nRows, int nCols);
    Matrix2(int nRows, int nCols, const T *inputData);
    Matrix2(const Matrix2<T> &inputMatrix);

    ~Matrix2();

    // configuration methods
    bool resize(int numRows, int numCols);

    // element access
    T GetElement(int row, int col);
    bool SetElement(int row, int col, T elementValue);
    int GetNumRows();
    int GetNumCols();

    // overload == operator
    bool operator==(const Matrix2<T> &rhs);

    // overload +, -, and * operators
    template <class U>
    friend Matrix2<U> operator+(const Matrix2<U> &lhs, const Matrix2<U> &rhs);
    template <class U>
    friend Matrix2<U> operator+(const U &lhs, const Matrix2<U> &rhs);
    template <class U>
    friend Matrix2<U> operator+(const Matrix2<U> &lhs, const U &rhs);

    template <class U>
    friend Matrix2<U> operator-(const Matrix2<U> &lhs, const Matrix2<U> &rhs);
    template <class U>
    friend Matrix2<U> operator-(const U &lhs, const Matrix2<U> &rhs);
    template <class U>
    friend Matrix2<U> operator-(const Matrix2<U> &lhs, const U &rhs);

    template <class U>
    friend Matrix2<U> operator*(const Matrix2<U> &lhs, const Matrix2<U> &rhs);
    template <class U>
    friend Matrix2<U> operator*(const U &lhs, const Matrix2<U> &rhs);
    template <class U>
    friend Matrix2<U> operator*(const Matrix2<U> &lhs, const U &rhs);

private:
    int Sub2Ind(int row, int col);

private:
    T *m_matrixData;
    int m_nRows, m_nCols, m_nElements;
};

// default constructor
template <class T>
Matrix2<T>::Matrix2()
{
    m_nRows = 1;
    m_nCols = 1;
    m_nElements = 1;
    m_matrixData = new T[m_nElements];
    m_matrixData[0] = 0.0;
}

// construct empty matrix
template <class T>
Matrix2<T>::Matrix2(int nRows, int nCols)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
    {
        m_matrixData[i] = 0.0;
    }
}

// construct from linear array
template <class T>
Matrix2<T>::Matrix2(int nRows, int nCols, const T *inputData)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_nElements = m_nRows * m_nCols;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
    {
        m_matrixData[i] = inputData[i];
    }
}

// copy constructor
template <class T>
Matrix2<T>::Matrix2(const Matrix2<T> &inputMatrix)
{
    m_nRows = inputMatrix.m_nRows;
    m_nCols = inputMatrix.m_nCols;
    m_nElements = inputMatrix.m_nElements;
    m_matrixData = new T[m_nElements];
    for (int i = 0; i < m_nElements; i++)
    {
        m_matrixData[i] = inputMatrix.m_matrixData[i];
    }
}

// destructor
template <class T>
Matrix2<T>::~Matrix2()
{
    if (m_matrixData != nullptr)
    {
        delete[] m_matrixData;
    }
}

// configuration method
template <class T>
bool Matrix2<T>::resize(int numRows, int numCols)
{
    m_nRows = numRows;
    m_nCols = numCols;
    m_nElements = m_nRows * m_nCols;
    delete[] m_matrixData;
    m_matrixData = new T[m_nElements];
    if (m_matrixData != nullptr)
    {
        for (int i = 0; i < m_nElements; i++)
            m_matrixData[i] = 0.0;
        return true;
    }
    else
    {
        return false;
    }
}

// element functions
template <class T>
T Matrix2<T>::GetElement(int row, int col)
{
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0)
    {
        return m_matrixData[linearIndex];
    }
    return 0.0;
}

template <class T>
bool Matrix2<T>::SetElement(int row, int col, T elementValue)
{
    int linearIndex = Sub2Ind(row, col);
    if (linearIndex >= 0)
    {
        m_matrixData[linearIndex] = elementValue;
        return true;
    }
    else
    {
        return false;
    }
}

template <class T>
int Matrix2<T>::GetNumRows()
{
    return m_nRows;
}

template <class T>
int Matrix2<T>::GetNumCols()
{
    return m_nCols;
}

// overloaded operators
// + operator
// matrix + matrix
template <class T>
Matrix2<T> operator+(const Matrix2<T> &lhs, const Matrix2<T> &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
    {
        tempResult[i] = lhs.m_matrixData[i] + rhs.m_matrixData[i];
    }
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar + matrix
template <class T>
Matrix2<T> operator+(const T &lhs, const Matrix2<T> &rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
    {
        tempResult[i] = lhs + rhs.m_matrixData[i];
    }
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix + scalar
template <class T>
Matrix2<T> operator+(const Matrix2<T> &lhs, const T &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
    {
        tempResult[i] = lhs.m_matrixData[i] + rhs;
    }
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// - operator
// matrix - matrix
template <class T>
Matrix2<T> operator-(const Matrix2<T> &lhs, const Matrix2<T> &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
    {
        tempResult[i] = lhs.m_matrixData[i] - rhs.m_matrixData[i];
    }
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar - matrix
template <class T>
Matrix2<T> operator-(const T &lhs, const Matrix2<T> &rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
    {
        tempResult[i] = lhs - rhs.m_matrixData[i];
    }
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix - scalar
template <class T>
Matrix2<T> operator-(const Matrix2<T> &lhs, const T &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
    {
        tempResult[i] = lhs.m_matrixData[i] - rhs;
    }
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// * operator
// matrix * matrix
template <class T>
Matrix2<T> operator*(const Matrix2<T> &lhs, const Matrix2<T> &rhs)
{
    int r_numRows = rhs.m_nRows;
    int r_numCols = rhs.m_nCols;
    int l_numRows = lhs.m_nRows;
    int l_numCols = lhs.m_nCols;

    if (l_numCols != r_numRows)
    {
        Matrix2<T> result(-1, -1);
        return result;
    }

    T *tempResult = new T[l_numRows * r_numCols];
    // loop through each row of LHS
    for (int lhsRow = 0; lhsRow < l_numRows; lhsRow++)
    {
        // loop through each column of RHS
        for (int rhsCol = 0; rhsCol < r_numCols; rhsCol++)
        {
            T elementResult = 0.0;
            // loop through each element of this LHS row
            for (int lhsCol = 0; lhsCol < l_numCols; lhsCol++)
            {
                // compute LHS linear index
                int lhsLinearIndex = (lhsRow * l_numCols) + lhsCol;
                // compute RHS linear index hased on LHS col (rhs row number is equal to the lhs column number)
                int rhsLinearIndex = (lhsCol * r_numCols) + rhsCol;

                // perform calculation
                elementResult += (lhs.m_matrixData[lhsLinearIndex] + rhs.m_matrixData[rhsLinearIndex]);
            }
            // store result
            int resultLinearIndex = (lhsRow * r_numCols) + rhsCol;
            tempResult[resultLinearIndex] = elementResult;
        }
    }
    Matrix2<T> result(l_numRows, r_numCols, tempResult);
    delete[] tempResult;
    return result;
}

// scalar * matrix
template <class T>
Matrix2<T> operator*(const T &lhs, const Matrix2<T> &rhs)
{
    int numRows = rhs.m_nRows;
    int numCols = rhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
    {
        tempResult[i] = lhs * rhs.m_matrixData[i];
    }
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// matrix * scalar
template <class T>
Matrix2<T> operator*(const Matrix2<T> &lhs, const T &rhs)
{
    int numRows = lhs.m_nRows;
    int numCols = lhs.m_nCols;
    int numElements = numRows * numCols;
    T *tempResult = new T[numElements];
    for (int i = 0; i < numElements; i++)
    {
        tempResult[i] = lhs.m_matrixData[i] * rhs;
    }
    Matrix2<T> result(numRows, numCols, tempResult);
    delete[] tempResult;
    return result;
}

// the == operator
template <class T>
bool Matrix2<T>::operator==(const Matrix2<T> &rhs)
{
    if ((this->m_nRows != rhs.m_nRows) && (this->m_nCols != rhs.m_nCols))
        return false;

    bool flag = true;
    for (int i = 0; i < this->m_nElements; ++i)
    {
        if (this->m_matrixData[i] != rhs.m_matrixData[i])
            flag = false;
    }
    return flag;
}

template <class T>
int Matrix2<T>::Sub2Ind(int row, int col)
{
    if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
    {
        return (row * m_nCols) + col;
    }
    return -1;
}

#endif