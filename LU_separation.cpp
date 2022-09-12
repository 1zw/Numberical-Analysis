/*
    Course Name: Numerical Analysis 
    Time: 2022/09/08 15：28
    Target:Caculate any sitution towards Ax=b
        Limition: A(n*n) & x(n*1) & b(n*1)
    Requirements: 
        1. Use LU separation; 
        2. for unique solution, cout it;
        3. for many solution, select one of them and cout;
        4. for no solution, describe it;
    Author: @Zhang Wei
    Student ID:2200120323
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
using namespace std;
/*
    matrix类
    Crout分解：U为单位上三角矩阵
    Dolliute分解：L为下三角矩阵
*/
class LU_solution;
void show(vector<double>& result)
{
    // vector<double> result = f.solution();
    for(int i=0;i<result.size();i++)
    {
        cout<<result[i]<<endl;
    }
}
class matrix
{
private:
    vector<vector<double>> A;
    vector<vector<double>> L;
    vector<vector<double>> U;
protected:
    int size;
    bool crout = false,dolliute = false;
public:
    matrix(){}
    // 初始化使用
    void init(vector<vector<double>> input){
        A = input;
        size = A.size();
        L.resize(size);
        U.resize(size);
    }
    void do_crout()
    {   
        // air = do_sum + lir ,(i = r,r+1,...n;)
        // arj = do_sum + lrr*urj (j = r,r+1,...n;r!=n)           
        for(int r = 0;r<size;r++)
        { 
            // 先更新L第r列
            for(int i = 0;i<size;i++)
            {   
                if(r>i)
                    L.at(i).push_back(0);
                else
                {   
                    double temp = A[i][r] - do_sum(i,r);
                    if(i==r && temp == 0)
                    {   
                        cout<<"crout separation error!";
                        return;
                    }
                    L.at(i).push_back(temp);
                }
            }
            // 再更新U第r行
            vector<double> temp;
            for(int j = 0;j<size;j++)
            {
                if(j<r)temp.push_back(0);
                else if(j==r) temp.push_back(1);
                else temp.push_back((A[r][j]-do_sum(r,j))/L[r][r]);
            }
            U.at(r) = temp;
        }
        crout = true;
    }
    void do_dolliute()
    {
        // air = do_sum + uir ,(i = r,r+1,...n;)
        // arj = do_sum + urr*lrj (j = r,r+1,...n;r!=n)
        for(int r = 0;r<size;r++)
        { 
            // 先更新U第r行数
            vector<double> temp;
            for(int j = 0;j<size;j++)
            {   
                if(j<r)temp.push_back(0);
                else 
                {   
                    if(j == r && A[r][j]-do_sum(r,j) == 0)
                    {
                        dolliute = false;
                        cout<<"dolliute separation error!";
                        return;
                    }
                    temp.push_back(A[r][j]-do_sum(r,j));
                }
            }
            
            U.at(r) = temp;
            // 再更新L第r列      
            for(int i = 0;i<size;i++)
            {
                if(r>i)
                    L.at(i).push_back(0);
                else if(r == i) L.at(i).push_back(1);
                else
                    L.at(i).push_back((A[i][r] - do_sum(i,r))/U[r][r]);
            }
            
            
        }
        dolliute = true;
    }
    tuple<vector<vector<double>>,vector<vector<double>>> get_lu()
    {
        return make_tuple(L,U);
    }
    double do_sum(int i,int j)
    {
        int minus = min(i,j);
        double sum = 0;
        for(int k=0;k<minus;k++)
        {
            sum += L[i][k]*U[k][j];
        }
        return sum;
    }
    void show()
    {
        if(crout || dolliute)
        {
            cout<<"-------------"<<endl;
            cout<<"here is the L"<<endl;
            cout<<"-------------"<<endl;
            for (int i = 0; i < size; i++)
            {   
                // cout<<"\t";
                for(int j = 0; j<size;j++)
                    cout<<L[i][j]<<" ";
                cout<<endl;
            }
            cout<<"-------------"<<endl;
            cout<<"here is the U"<<endl;
            cout<<"-------------"<<endl;
            for (int i = 0; i < size; i++)
            {
                // cout<<"\t";
                for(int j = 0; j<size;j++)
                    cout<<U[i][j]<<" ";
                cout<<endl;
            }
        }
        else
        {
            cout<<"------------------------------"<<endl;
            cout<<"this matrix has no LU solution"<<endl;
            cout<<"------------------------------"<<endl;
        }
    }
    // 声明LU_solution 类为自己的友元
    friend LU_solution;
};
// solution类，完成AX=b的计算
class LU_solution
{
private:
    matrix A_;
    vector<double> b_,x_;
    vector<vector<double>> L_;
    vector<vector<double>> U_;
public:
    LU_solution(){};
    // 完成初始化
    void init(vector<vector<double>> a,vector<double> b)
    {
        A_.init(a);
        A_.do_crout();
        tuple<vector<vector<double>>,vector<vector<double>>> turtle = A_.get_lu();
        L_ = get<0>(turtle);U_ = get<1>(turtle);
        b_ = b;
    }
    vector<double> solution()
    {
        // Ly = b
        vector<double> y;
        for(int i = 0;i<A_.size;i++)
        {
            double sum = 0;
            for(int j=0;j<y.size();j++)
            {
                sum += L_[i][j]*y[j];
            }
            y.push_back((b_.at(i)-sum)/L_[i][i]);
        }
        // Ux = y
        x_.resize(A_.size);
        for(int i = A_.size-1;i>=0;i--)
        {
            double sum = 0;
            for(int j=1;j<x_.size()-i;j++)
            {
                sum += U_[i][i+j]*x_[i+j];
            }
            x_.at(i) = ((y.at(i)-sum)/U_[i][i]);
        }
        return x_;
    }
};

int main(int argc, char const *argv[])
{   
    vector<vector<double>> a={{1,1,1,1},
                                {1,2,3,4},
                                {1,10,100,1000},
                                {2,5,9,8}};
    vector<double> b = {1,1,1,1};
    LU_solution f;
    f.init(a,b);
    vector<double> result = f.solution();
    show(result);
    return 0;
}
/*
matrix类测试案例：
1,1,1,1
1,2,3,4
1,10,100,1000
2,5,9,8
-------------
Dolliute
-------------
here is the L
-------------
1 0 0 0 
1 1 0 0 
1 9 1 0 
2 3 0.0123457 1 
-------------
here is the U
-------------
1 1 1 1 
0 1 2 3 
0 0 81 972 
0 0 0 -15 
-------------
Crout
-------------
here is the L
-------------
1 0 0 0 
1 1 0 0 
1 9 81 0 
2 3 1 -15 
-------------
here is the U
-------------
1 1 1 1 
0 1 2 3 
0 0 1 12 
0 0 0 1 
----------------
LU_solution类测试
----------------
在b = {1,1,1,1}时候
result = 
0.333333
1.4
-0.8
0.0666667
*/
