/*
    Course Name: Numerical Analysis 
    Time: 2022/09/01 21:55
    Target:Caculate any sitution towards Ax=b
        Limition: A(n*n) & x(n*1) & b(n*1)
    Required: 
        1. Stable; 
        2. unique solution, cout it;
        3. many solution, select one of them and cout;
        4. no solution, describe it;
    Author: @Zhang Wei
    Student ID:2200120323
*/
/*
----------------说明：本代码实现了两个类-------------
矩阵类MatrixRC可以实现对矩阵的初等行变换和只有交换的列变换
     solution类可以通过主元高斯法则求解出x向量
-------------------------------------------------
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
namespace Matrix
{
template <class T>
class MatrixRC
{   
private:
    std::vector<std::vector<T>> matrix;
    int row_;
    int column_;
    bool is_inited = false;
    int rank_;
    std::vector<std::vector<double>> initial_matrix;
    std::vector<std::vector<int>> changed_list;
public:
    void init(std::vector<std::vector<T>>& v)
    {
        matrix = v;
        row_ = v.size();
        column_ = v.at(0).size();
        is_inited = true;
        std::tuple<int,std::vector<std::vector<double>>> turtle = Caculate_Rank();
        rank_ = std::get<0>(turtle);
        initial_matrix = std::get<1>(turtle);
    }
    void show()
    {   
        std::cout<<"-----------------------------------"<<std::endl;
        std::cout<<"Matrix\trow:"<<row_<<"\tcolumn:"<<column_<<std::endl;
        std::cout<<"\tRank:"<<rank_<<std::endl;
        std::cout<<"-----------------------------------"<<std::endl;

        for(int i=0;i<row_;i++)
        {
            for(int j=0;j<column_;j++)
                std::cout<<matrix[i][j]<<" ";
            std::cout<<"\n";
        }

        std::cout<<"-----------------------------------"<<std::endl;
        std::cout<<"\tAfter Initialized\t"<<std::endl;
        std::cout<<"-----------------------------------"<<std::endl;
        for(int i=0;i<row_;i++)
        {
            for(int j=0;j<column_;j++)
                std::cout<<initial_matrix[i][j]<<" ";
            std::cout<<"\n";
        }
        if(!changed_list.empty())
        {
            std::cout<<"-----------------------------------"<<std::endl;
            std::cout<<"\tChanged List\t"<<std::endl;
            std::cout<<"-----------------------------------"<<std::endl;
            for(int i=0;i<changed_list.size();i++)
            {   
                std::cout<<"("<<changed_list.at(i).at(0)<<","<<changed_list.at(i).at(1)<<")"<<std::endl;
            }
        }
    }
    // 本函数返回false的情况：
    //    temp[here][here]>=剩下所有元素    
    bool transform_matrix(std::vector<std::vector<double>>& temp,int here)
    {
        int temp_row = temp.size();
        int temp_column = temp.at(0).size();
        double max_number = temp.at(here).at(here);
        int changed_column_number,changed_row_number;
        bool has_changed = false;
        if(here < std::min(temp_row,temp_column))
        {   
            for(int i=here;i<temp_row;i++)
            for(int j=here;j<temp_column;j++)
            {
                if(fabs(max_number)<fabs(temp[i][j]))
                {
                    has_changed = true;
                    max_number = temp[i][j];
                    changed_column_number = j;
                    changed_row_number = i;
                }
            }
            if(has_changed)
            {
                if(changed_row_number!=here)
                {
                    std::vector<double> tp = temp.at(changed_row_number);
                    temp.at(changed_row_number) = temp.at(here);
                    temp.at(here) = tp;
                }
                if(changed_column_number!=here)
                {
                    std::vector<int> who_changed = {here,changed_column_number};
                    changed_list.push_back(who_changed);
                    for(int k=0;k<temp_row;k++)
                    {
                        double temp_number = temp[k][changed_column_number];
                        temp[k][changed_column_number] = temp[k][here];
                        temp[k][here] = temp_number;
                    }
                }
            }
        }
        return has_changed;
    }
    void line_transformation(std::vector<std::vector<double>>& temp,int here)
    {
        int temp_row = temp.size();
        int temp_column = temp.at(0).size();
        double shoot;
        if(here+1<temp_row)
        {
            for(int i=here+1;i<temp_row;i++)
            {   
                bool has_seen = false;
                if(!has_seen)
                {
                    shoot = temp[i][here]/temp[here][here];
                    has_seen = true;
                }
                for(int j=here;j<temp_column;j++)                    
                    temp[i][j] -= shoot*temp[here][j];  
            }
        }
    }
    std::tuple<int,std::vector<std::vector<double>>> Caculate_Rank()
    {   
        if(is_inited)
        {
            int full_rank = std::min(row_,column_);
            std::vector<std::vector<double>> temp = matrix;
            for(int i=0;i<temp.size();i++)
            {
                bool has_changed = transform_matrix(temp,i);
                if(has_changed)
                    line_transformation(temp,i);
                else{
                    if(temp[i][i])
                        line_transformation(temp,i);
                    else
                        return std::make_tuple(i,temp);
                }
            }
            return std::make_tuple(full_rank,temp);
        }
        else{
            std::cerr<<"your matrix haven't been inited!"<<std::endl;
            throw "matrix need been inited first!";
        }

    }
    int get_rank(){return rank_;}
    std::vector<std::vector<double>> get_initial_matrix(){return initial_matrix;}
    std::vector<std::vector<int>> get_changed_list(){return changed_list;}
};
class Axb_solution
{
private:
    std::vector<std::vector<double>> A_;
    std::vector<double> b_;
    std::vector<std::vector<double>> A1b_;
    MatrixRC<double> matrix_A;
    MatrixRC<double> matrix_A1b;
    std::vector<double> result_;
    bool is_solved = false;
    int count=0;
    int full_rank_;
public:
    Axb_solution(){};
    void init(std::vector<std::vector<double>>& A,std::vector<double>& b)
    {   
        A_ = A;b_ = b;
        full_rank_ = std::min(A.size(),A.at(0).size());
        for(int i=0;i<A_.size();i++)
        {
            std::vector<double> temp = A_.at(i);
            temp.push_back(b_.at(i));
            A1b_.push_back(temp);
        }
        matrix_A.init(A_);
        // matrix_A.show();
        matrix_A1b.init(A1b_);
        matrix_A1b.show();
    }
    bool have_solution()
    {
        if(matrix_A.get_rank()!=matrix_A1b.get_rank())
            return false;
        else
            return true;
    }
    std::vector<double> get_result(){return result_;}
    void restore_sequence(){
       std::vector<std::vector<int>> sequence = matrix_A1b.get_changed_list();
       if(!sequence.empty())
       {
        for(int i=sequence.size()-1;i>=0;i--)
        {
            int a = sequence.at(i).at(0);
            int b = sequence.at(i).at(1);
            double temp = result_.at(a);
            result_.at(a) = result_.at(b);
            result_.at(b) = temp;
        }
        is_solved = true;
       }
        else
        return;
    }
    bool solution()
    {   
        std::vector<std::vector<double>> initial_mat = matrix_A1b.get_initial_matrix();
        std::vector<double> reverse_temp;
        result_.clear();
        if(have_solution())
        {   
            /*full_rank*/
            if(matrix_A.get_rank() == full_rank_)
            {   
                std::cout<<"-----------------------------------"<<std::endl;
                std::cout<<"this linear equations have unique solution"<<std::endl;
                std::cout<<"-----------------------------------"<<std::endl;

                for(int i=matrix_A.get_rank()-1;i>=0;i--)
                {   
                    double sum = 0;   
                    for(int j=0;j<reverse_temp.size();j++)
                        sum -= (initial_mat.at(i).at(i+j+1))*reverse_temp.at(reverse_temp.size()-1-j);
                
                    double x = (initial_mat.at(i).at(full_rank_) + sum)/initial_mat.at(i).at(i);
                    reverse_temp.push_back(x);
                }
            }
            else{
                int unknown_number_count = full_rank_ - matrix_A1b.get_rank();
                std::cout<<"-----------------------------------"<<std::endl;
                std::cout<<"this linear equations have infinite solution"<<std::endl;
                std::cout<<"the "<<count+1<<" time to ask for answer"<<std::endl;
                std::cout<<"-----------------------------------"<<std::endl;
                for(int i = A1b_.size()-1;i>=0;i--)
                {   
                    double sum = 0;
                    if(i+1>matrix_A1b.get_rank())
                    {   
                        // sum -= count; 
                        reverse_temp.push_back(count);
                    }
                    else{
                        for(int j=0;j<reverse_temp.size();j++)
                            sum -= (initial_mat.at(i).at(i+j+1))*reverse_temp.at(reverse_temp.size()-1-j);
                        double x = (initial_mat.at(i).at(full_rank_) + sum)/initial_mat.at(i).at(i);
                        // sum -= initial_mat.at(i).at(i)*x;
                        reverse_temp.push_back(x);
                    }
                }
                count++;
            }
            for(int i=reverse_temp.size()-1;i>=0;i--)
                    result_.push_back(reverse_temp.at(i));
            restore_sequence();
            return true;
        }
        else{
            std::cout<<"-----------------------------------"<<std::endl;
            std::cerr<<"sorry,this matrix don't any have solution"<<std::endl;
            std::cout<<"-----------------------------------"<<std::endl;
            return false;
        }
    }
};
}
int main(int argc, char const *argv[])
{
    // Matrix::MatrixRC<double> m;
    // std::vector<std::vector<double>> v={{0,0,0},{0,1,1},{0,1,1}};
    // m.init(v);
    // m.show();
    Matrix::Axb_solution solver;
    std::vector<std::vector<double>> A = {{0,0,0},{0,1,1},{0,1,1}};
    std::vector<double> b = {0,1,0};
    solver.init(A,b);
    // int i=0;
    // while(i<3)
    // {
    if(solver.solution())
    {
        std::vector<double> result = solver.get_result();
        std::cout<<"result:"<<std::endl;
        for(int i=0;i<result.size();i++)
            std::cout<<"\t"<<result.at(i)<<std::endl;
    }
    // i++;
    // }
    return 0;
}
/*---------------MatrixRC测试用例-------------------*/
/*-------------------------------------------------
{{0,5,3},{2,3,1,},{0,9,1}};
{{1,5,3},{2,3,1,},{0,9,1}};
{{0,1,2},{3,4,5},{6,7,8}};
{{1,1,1},{2,2,2},{3,3,3},{4,4,4},{5,5,5}};
{{1,1,1,1,1},{2,2,2,2,2},{3,3,3,3,4}};
{{0,1,1},{0,1,1},{0,1,1}}
{{0,0,0},{0,1,1},{0,1,1}};
以{{0,0,0},{0,1,1},{0,1,1}}为例子，输入结果为
-----------------------------------
Matrix  row:3   column:3
        Rank:1
-----------------------------------
0 0 0 
0 1 1 
0 1 1 
-----------------------------------
        After Initialized
-----------------------------------
1 0 1 
0 0 0 
0 0 0 
-----------------------------------
        Changed List
-----------------------------------
(0,1)
符合预期
----------------------------------------------*/
/*--------------solution测试用例---------------*/
/*
1.
{{0.001,2.000,3},{-1.000,3.712,4.623},{-2,1.072,5.643}};
{1,2,3};
计算结果为：
-----------------------------------
this linear equations has unique solution
-----------------------------------
result:
        -0.490396
        -0.0510352
        0.36752
2.
{{0,0,0},{0,1,1},{0,1,1}};
{0,1,1};
-----------------------------------
this linear equations have infinite solution
the 1 time to ask for answer
-----------------------------------
result:
        0
        1
        0
-----------------------------------
this linear equations have infinite solution
the 2 time to ask for answer
-----------------------------------
result:
        1
        0
        1
-----------------------------------
this linear equations have infinite solution
the 3 time to ask for answer
-----------------------------------
result:
        2
        -1
        2
3. 
{{0,0,0},{0,1,1},{0,1,1}}
{0,1,0}
-----------------------------------
sorry,this matrix don't any have solution
-----------------------------------
*/


