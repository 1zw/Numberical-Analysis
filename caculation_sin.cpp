/*
    Course Name: Numerical Analysis 
    Time: 2022/08/22 15:55
    Target:Caculate sin1 and sin1000
    Ideas: 
        1. Taylor expansion
            sin taylor function:
                sigema(-1)^m-1 * x^(2m-1) /(2m-1)!
        2. function periodicity
    Author: @Zhang Wei
*/
#include <iostream>
#include <cmath>
using namespace std;
class sin_solution
{
private:
    const double pi = 3.1415926;
    const double precision = 1e-5;
    double figure;
    double error_line;
    double result;
    int iteration_time;
public:
    sin_solution(double number):
    iteration_time(1),
    figure(number),
    result(0.0)
    {};
    ~sin_solution(){};
   void limit_figure()
   {
        while(1)
        {
            if(figure>=2*pi)
                figure-=2*pi;
            else if (figure<0)
                figure+=2*pi;
            else
                break;
        }
   }
   int get_factorial(int n)
   {
    if(n == 0) return 1;
    return n*get_factorial(n-1);
   }
   double solution()
   {
        limit_figure();
        error_line = figure;
        cout<<"figure = "<<figure<<endl;
        while(fabs(error_line) >= precision)
        {   
            cout<<"---------------iteration---------------"<<endl;
            cout<<"time:"<<iteration_time<<endl;
            error_line = pow(-1,iteration_time-1)*pow(figure,2*iteration_time-1) \
                        /get_factorial(2*iteration_time-1);
            result += error_line;
            iteration_time++;
            cout<<"now_error_line = "<<error_line<<endl;
            cout<<"now_result = "<<result<<endl;

        }
        return result;
   }
};
int main(int argc, char const *argv[])
{
    sin_solution sin_1(1);
    sin_solution sin_1000(1000);
    double result_1 = sin_1.solution();
    double result_1000 = sin_1000.solution();
    cout<<"sin(1) = "<<result_1<<endl;
    cout<<"sin(1000) = "<<result_1000<<endl;

    return 0;
}
/*
sin(1)
figure = 1
---------------iteration---------------
time:1
now_error_line = 1
now_result = 1
---------------iteration---------------
time:2
now_error_line = -0.166667
now_result = 0.833333
---------------iteration---------------
time:3
now_error_line = 0.00833333
now_result = 0.841667
---------------iteration---------------
time:4
now_error_line = -0.000198413
now_result = 0.841468
---------------iteration---------------
time:5
now_error_line = 2.75573e-06
now_result = 0.841471
sin(1) = 0.841471

sin(1000)
figure = 0.973553
---------------iteration---------------
time:1
now_error_line = 0.973553
now_result = 0.973553
---------------iteration---------------
time:2
now_error_line = -0.15379
now_result = 0.819763
---------------iteration---------------
time:3
now_error_line = 0.00728815
now_result = 0.827051
---------------iteration---------------
time:4
now_error_line = -0.00016447
now_result = 0.826887
---------------iteration---------------
time:5
now_error_line = 2.16508e-06
now_result = 0.826889
sin(1000) = 0.826889
*/