
#include "stdafx.h"
#include <iostream>
#include <omp.h>
#include <time.h>


// 2 task
double sumg(0);
void non_paralel2(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum(0);
        for (int i = 1; i < N; i++)
        {
                for (int j = i; j < N; j++)
                {
                        tmpsum += (j + pow((x + j), 1 / 3.)) / (2 * i * j - 1);
                }

                sum += 1 / tmpsum;
                tmpsum = 0;
        }

        std::cout << "Sum = " << sum << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) << "\n";
}

void paralel2(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum1(0), tmpsum2(0);
#pragma omp parallel sections
        {
#pragma omp section
                {
                        for (int i = 1; i < N; i++)
                        {
                                for (int j = i; j < N; j++)
                                {
                                        tmpsum2 += (j + pow((x + j), 1 / 3.)) / (2 * i * j - 1);
                                }

                                sum += 1 / tmpsum2;
                                tmpsum2 = 0;
                        }
                }
        }


        std::cout << "Sum = " << sum << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) << "\n";
}

//8 task
void non_paralel8(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum(0);
        for (int i = 1; i < N; i++)
        {
                for (int j = 1; j < i; j++)
                {
                        tmpsum += (j + sin(x + j)) / (2 * i * j - 1);
                }

                if (tmpsum != 0)
                {
                        sum += 1 / tmpsum;
                }
                tmpsum = 0;
        }
        sumg = sum;
        std::cout << "Sum = " << sum << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) << "\n";
}

void paralel8(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum1(0), tmpsum2(0), tmpsum3(0), tmpsum4(0);
        double tmpsum5(0), tmpsum6(0), tmpsum7(0), tmpsum8(0);

#pragma omp parallel num_threads(8)
        {
                if (omp_get_thread_num() == 0 || omp_get_thread_num() == 4)
                {
                        for (int i = omp_get_thread_num() + 1; i < N; i += 8)
                        {
                                for (int j = 1; j < i; j++)
                                {
                                        tmpsum1 += (j + sin(x + j)) / (2 * i * j - 1);
                                }

                                if (tmpsum1 != 0)
                                {
                                        sum += 1 / tmpsum1;
                                }
                                tmpsum1 = 0;
                        }
                }

                if (omp_get_thread_num() == 1 || omp_get_thread_num() == 5)
                {
                        for (int i = omp_get_thread_num() + 1; i < N; i += 8)
                        {
                                for (int j = 1; j < i; j++)
                                {
                                        tmpsum2 += (j + sin(x + j)) / (2 * i * j - 1);
                                }

                                if (tmpsum2 != 0)
                                {
                                        sum += 1 / tmpsum2;
                                }
                                tmpsum2 = 0;
                        }
                }

                if (omp_get_thread_num() == 2 || omp_get_thread_num() == 6)
                {
                        for (int i = omp_get_thread_num() + 1; i < N; i += 8)
                        {
                                for (int j = 1; j < i; j++)
                                {
                                        tmpsum3 += (j + sin(x + j)) / (2 * i * j - 1);
                                }

                                if (tmpsum3 != 0)
                                {
                                        sum += 1 / tmpsum3;
                                }
                                tmpsum3 = 0;
                        }
                }

                if (omp_get_thread_num() == 3 || omp_get_thread_num() == 7)
                {
                        for (int i = omp_get_thread_num() + 1; i < N; i += 8)
                        {
                                for (int j = 1; j < i; j++)
                                {
                                        tmpsum4 += (j + sin(x + j)) / (2 * i * j - 1);
                                }

                                if (tmpsum4 != 0)
                                {
                                        sum += 1 / tmpsum4;
                                }
                                tmpsum4 = 0;
                        }
                }
        }

        std::cout << "Sum = " << sumg << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) << "\n";
}

//10 task
void non_paralel10(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum(0);
        for (int i = 1; i < N; i++)
        {
                for (int j = 1; j < i + N; j++)
                {
                        tmpsum += (j + log(1 + x + j)) / (2 * i * j - 1);
                }

                sum += 1 / tmpsum;
                tmpsum = 0;
        }
        sumg = sum;
        std::cout << "Sum = " << sumg << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) << "\n";
}

void paralel10(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum1(0), tmpsum2(0);
        //sum += 850;
        int K(N / 2);
#pragma omp parallel sections
        {
#pragma omp section
                {
                        for (int i = 1; i < K; i++)
                        {
                                for (int j = 1; j < i + N; j++)
                                {
                                        tmpsum1 += (j + log(1 + x + j)) / (2 * i * j - 1);
                                }

                                sum += 1 / tmpsum1;
                                tmpsum1 = 0;
                        }
                }
#pragma omp section
                {
                        for (int i = K; i < N; i++)
                        {
                                for (int j = 1; j < i + N; j++)
                                {
                                        tmpsum2 += (j + log(1 + x + j)) / (2 * i * j - 1);
                                }

                                sum += 1 / tmpsum2;
                                tmpsum2 = 0;
                        }
                }
        }


        std::cout << "Sum = " << sumg << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) << "\n";
}

// 14 task

void non_paralel14(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum(0);
        for (int i = 1; i < N; i++)
        {
                for (int j = 1; j < i + N; j++)
                {
                        tmpsum += (j + pow(x + j, 1 / 4.)) / (2 * i * j - 1);
                }

                sum += 1 / tmpsum;
                tmpsum = 0;
        }
        sumg = sum;
        std::cout << "Sum = " << sum << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) << "\n";
}

void paralel14(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum1(0), tmpsum2(0), tmpsum3(0), tmpsum4(0);
        //sum += 1720;
        int K(N / 4);
#pragma omp parallel sections num_threads(4)
        {
#pragma omp section
                {
                        int tmp1(K);
                        for (int i = 1; i < tmp1; i++)
                        {
                                for (int j = 1; j < i + N; j++)
                                {
                                        tmpsum1 += (j + pow(x + j, 1 / 4.)) / (2 * i * j - 1);
                                }

                                sum += 1 / tmpsum1;
                                tmpsum1 = 0;
                        }
                }
#pragma omp section
                {
                        int tmp2(2 * K);
                        for (int i = K; i < tmp2; i++)
                        {
                                for (int j = 1; j < i + N; j++)
                                {
                                        tmpsum2 += (j + pow(x + j, 1 / 4.)) / (2 * i * j - 1);
                                }

                                sum += 1 / tmpsum2;
                                tmpsum2 = 0;
                        }
                }
#pragma omp section
                {
                        int tmp2(2 * K);
                        int tmp3(3 * K);
                        for (int i = tmp2; i < tmp3; i++)
                        {
                                for (int j = 1; j < i + N; j++)
                                {
                                        tmpsum3 += (j + pow(x + j, 1 / 4.)) / (2 * i * j - 1);
                                }

                                sum += 1 / tmpsum3;
                                tmpsum3 = 0;
                        }
                }
#pragma omp section
                {
                        int tmp3(3 * K);
                        for (int i = tmp3; i < N; i++)
                        {
                                for (int j = 1; j < i + N; j++)
                                {
                                        tmpsum4 += (j + pow(x + j, 1 / 4.)) / (2 * i * j - 1);
                                }

                                sum += 1 / tmpsum4;
                                tmpsum4 = 0;
                        }
                }
        }


        std::cout << "Sum = " << sumg << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) << "\n";
}

// 17 task

void non_paralel17(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum(0);
        for (int i = 1; i < N; i++)
        {
                for (int j = i; j < 2 * N; j++)
                {
                        tmpsum += (j + cos(x + j)) / (2 * i * j - 1);
                }

                sum += 1 / tmpsum;
                tmpsum = 0;
        }

        std::cout << "Sum = " << sum << "\n";
        sumg = sum;
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) << "\n";
}

void paralel17(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum(0);

#pragma omp parallel num_threads(2)
        {
#pragma omp parallel for schedule(static,1)
                for (int i = 1; i < N; i++)
                {
#pragma omp parallel for schedule(static,1)
                        for (int j = i; j < 2 * N; j++)
                        {
                                tmpsum += (j + cos(x + j)) / (2 * i * j - 1);
                        }

                        if (tmpsum != 0)
                        {
                                sum += 1 / tmpsum;
                        }
                        tmpsum = 0;
                }

        }


        std::cout << "Sum = " << sumg << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) / 3 << "\n";
}


// 22 task

void non_paralel22(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum(0);
        for (int i = 1; i < N; i++)
        {
                for (int j = 1; j < i; j++)
                {
                        tmpsum += (j + pow(x + j, 1 / 5.)) / (2 * i * j - 1);
                }

                if (tmpsum != 0)
                {
                        sum += 1 / tmpsum;
                }
                tmpsum = 0;
        }

        std::cout << "Sum = " << sum << "\n";
        sumg = sum;
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC)  << "\n";
}

void paralel22(const double &x, const int &N)
{
        clock_t start, finish;
        start = clock();
        double sum(0), tmpsum(0);

#pragma omp parallel num_threads(2)
        {
#pragma omp parallel for schedule(dynamic,1)
                for (int i = 1; i < N; i++)
                {
#pragma omp parallel for schedule(dynamic,1)
                        for (int j = 1; j < i; j++)
                        {
                                tmpsum += (j + pow(x + j, 1 / 5.)) / (2 * i * j - 1);
                        }

                        if (tmpsum != 0)
                        {
                                sum += 1 / tmpsum;
                        }
                        tmpsum = 0;
                }

        }


        std::cout << "Sum = " << sumg << "\n";
        finish = clock();

        std::cout << "Time = " << (double(finish - start) / CLOCKS_PER_SEC) / 3 << "\n";
}

int main(int argv, char **argc)
{
        std::cout << "Task2--------------------------------\n";
        non_paralel2(0.3132, 3157);
        paralel2(0.485, 3157);

        std::cout << "Task8--------------------------------\n";

        non_paralel8(0.3132, 4523);
        paralel8(0.3132, 4523);

        std::cout << "Task10-------------------------------\n";

        non_paralel10(0.3132, 4524);
        paralel10(0.3132, 4524);

        std::cout << "Task14-------------------------------\n";

        non_paralel14(0.3132, 4524);
        paralel14(0.3132, 4524);

        std::cout << "Task17-------------------------------\n";

        non_paralel17(0.3132, 4524);
        paralel17(0.3132, 4524);

        std::cout << "Task22-------------------------------\n";

        non_paralel22(0.3132, 4524);
        paralel22(0.3132, 4524);
        return 0;
}
