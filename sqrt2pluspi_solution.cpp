/*
    Course Name: Numerical Analysis 
    Time: 2022/08/29 10.25
    Target:Calutate sqrt2 and pi
    Ideas: 
        1. dichotomy    
        2. Newton iteration
        3. secant iteration
        4. poplynomial 1)sinx = 0 (1/2*pi<x<3/2*pi)
                       2)x^2-2 =0 
    Author: @Zhang Wei
*/
#include <iostream>
#include <cmath>
using namespace std;
class dichotomy
{
private:
    bool pi_flag_,sqrt2_flag_;
    int iteration_times_;
public:
    dichotomy(bool pi_flag = false,bool sqrt2_flag = false,int iteration_times = 10 ):
    //选择计算对象,默认值为false
    pi_flag_(pi_flag),sqrt2_flag_(sqrt2_flag),iteration_times_(iteration_times){
        cout<<"---------------dichotomy begin------------"<<endl;
    };
    /* 异或检验 */
    bool xor_check(bool a, bool b)
    {
        if(((a==false) && (b==true)) || ((a==true) && (b==false)))
            return true;
        else
            return false;
    }
    /* 计算f(x) = x^2-2 的结果，正返回true，负返回false */
    bool sqrt2_polynomial(double figure){
        double this_result = pow(figure,2)-2;
        if(this_result > 0) return true;
        return false; 
    }
    bool pi_polynomial(double figure)
    {
        double this_result = sin(figure);
        if(this_result > 0) return true;
        return false;
    }
    void solution()
    {
        if(pi_flag_)
        {
            cout<<"pi use "<<pi_solution()<<endl;
        }
        if(sqrt2_flag_)
        {
            cout<<"sqrt2 use "<<sqrt2_solution()<<endl;
        }
    }
    double pi_solution()
    {
        double init_left = 1.58,init_right =  4.71; // 初始化二分左值和右值
        double mid_figure;
        double last_mid_figure = 1e9;
        for(int i = 0;i<iteration_times_;i++)
        {   
            mid_figure = 0.5*(init_left + init_right);
            bool in_left_interval  = xor_check(pi_polynomial(init_left),pi_polynomial(mid_figure));
            bool in_right_interval = xor_check(pi_polynomial(mid_figure),pi_polynomial(init_right));
            if(in_left_interval)
                init_right = mid_figure;
            if(in_right_interval)
                init_left = mid_figure;
            if(fabs(last_mid_figure-mid_figure) < 1e-6)
            {
                cout<<"iteration times:"<<i<<"  result is  ";
                break;
            }
            last_mid_figure = mid_figure;
        }
        return mid_figure;
    }
    double sqrt2_solution()
    {
        double init_left = 0.0,init_right = 2.0;
        double mid_figure;
        double last_mid_figure = 1e9;
        for(int i = 0;i<iteration_times_;i++)
        {   
            mid_figure = 0.5*(init_left + init_right);
            bool in_left_interval  = xor_check(sqrt2_polynomial(init_left),sqrt2_polynomial(mid_figure));
            bool in_right_interval = xor_check(sqrt2_polynomial(mid_figure),sqrt2_polynomial(init_right));
            if(in_left_interval)
                init_right = mid_figure;
            if(in_right_interval)
                init_left = mid_figure;
            if(fabs(last_mid_figure-mid_figure) < 1e-6)
            {
                cout<<"iteration times:"<<i<<"  result is  ";
                break;
            }
            last_mid_figure = mid_figure;
        }
        return mid_figure;
    }
};
class Newton
{
private:
    bool pi_flag_,sqrt2_flag_;
    int iteration_times_;
    double precision_ = 1e-6;
public:
    Newton(bool pi_flag = false,bool sqrt2_flag = false,int iteration_times = 100):
    pi_flag_(pi_flag),sqrt2_flag_(sqrt2_flag),iteration_times_(iteration_times)
    {
        cout<<"---------------Newton begin------------"<<endl;
    };
    double pi_polynomial(double result)
    {
        return result-sin(result)/cos(result);
    }
    double sqrt2_polynomial(double result)
    {
        return result-(pow(result,2)-2)/(2*result);
    }
    void solution()
    {
        if(pi_flag_)
        {
            cout<<"pi use "<<pi_solution()<<endl;
        }
        if(sqrt2_flag_)
        {
            cout<<"sqrt2 use "<<sqrt2_solution()<<endl;
        }
    }
    double pi_solution()
    {   
        double result = 2,last_result = 1e9;
        int i;
        for(i=0;i<iteration_times_;i++)
        {   
            result = pi_polynomial(result);
            if(fabs(last_result-result) < precision_)
                break;
        }
        cout<<"iteration times:"<<i<<"  result is  ";

        return result;
    }
    double sqrt2_solution()
    {
        double result = 2,last_result = 1e9;
        int i;
        for(i=0;i<iteration_times_;i++)
        {
            result = sqrt2_polynomial(result);
            if(fabs(last_result-result) < precision_)
            {
                break;
            }
        }
        cout<<"iteration times:"<<i<<"  result is  ";
        return result;
    }
};
class secant
{
private:
    bool pi_flag_,sqrt2_flag_;
    int iteration_times_;
    double precision_ = 1e-6;
    double step_length = 1e-3;
public:
    secant(bool pi_flag=false,bool sqrt2_flag=false,int iteration_times=100):
    pi_flag_(pi_flag),sqrt2_flag_(sqrt2_flag),iteration_times_(iteration_times)
    {
        cout<<"-----------secant begin------------"<<endl;
    }
    double pi_polynomial(double result)
    {
        return result-(sin(result)*step_length)/(sin(result)-sin(result-step_length));
    }
    double sqrt2_polynomial(double result)
    {
        return result-(pow(result,2)-2)*step_length/(pow(result,2)-2-(pow(result-step_length,2)-2));
    }
    double pi_solution()
    {
        double result = 2,last_result = 1e9;
        int i;
        for(i=0;i<iteration_times_;i++)
        {   
            result = pi_polynomial(result);
            if(fabs(last_result-result) < precision_)
                break;
        }
        cout<<"iteration times:"<<i<<"  result is  ";

        return result;
    }
    double sqrt2_solution()
    {
        double result = 2,last_result = 1e9;
        int i;
        for(i=0;i<iteration_times_;i++)
        {
            result = sqrt2_polynomial(result);
            if(fabs(last_result-result) < precision_)
            {
                break;
            }
        }
        cout<<"iteration times:"<<i<<"  result is  ";
        return result;
    }
    void solution()
    {
        if(pi_flag_)
        {
            cout<<"pi use "<<pi_solution()<<endl;
        }
        if(sqrt2_flag_)
        {
            cout<<"sqrt2 use "<<sqrt2_solution()<<endl;
        }
    }

};
int main(int argc, char const *argv[])
{
    dichotomy dichotomy_solutuion(true,true,100);
    dichotomy_solutuion.solution();
    Newton newton_solution(true,true,100);
    newton_solution.solution();
    secant secant_solution(true,true,100);
    secant_solution.solution();
    return 0;
}
/*
代码运行结果
dichotomy:(二分法)
---------------dichotomy begin------------
pi use iteration times:21  result is  3.14159
sqrt2 use iteration times:20  result is  1.41421
Newton iteration:(牛顿迭代)
---------------Newton begin------------
pi use iteration times:100  result is  3.14159
sqrt2 use iteration times:100  result is  1.41421
secant iteration：(割线法)
-----------secant begin------------
pi use iteration times:100  result is  3.14159
sqrt2 use iteration times:100  result is  1.41421
*/