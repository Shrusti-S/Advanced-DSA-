/*                         CS6013 - ADSA - ASSIGNMENT 1 - MULTIPLYING POLYNOMIALS USING FFT

Submitted By,
    Name : SHRUSTI
    Roll Number : CS22MTECH11017
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>

using namespace std;

int poly_1[16], poly_2[16], deg1, deg2;

// returns the product of the complex numbers a + ib and c + id
double *multiply_complex_num(double a, double b, double c, double d)
{
    double *mult = new double[2];
    double mult_real = a * c - b * d;
    double mult_img = a * d + b * c;

    mult[0] = mult_real;
    mult[1] = mult_img;

    return mult;
}

// returns the sum of the complex numbers a + ib and c + id
double *add_complex_num(double a, double b, double c, double d)
{
    double *sum = new double[2];
    double sum_real = a + c;
    double sum_img = b + d;

    sum[0] = sum_real;
    sum[1] = sum_img;

    return sum;
}

// The function takes N as an argument and computes the N-th roots of unity
double **find_complex_roots(int N)
{
    double **roots = new double *[N];
    for (int i = 0; i < N; i++)
    {
        roots[i] = new double[2];
    }

    double theta = (2 * ((double)22 / 7)) / N;

    for (int i = 0; i < N; i++)
    {
        roots[i][0] = cos(theta * i); // real
        roots[i][1] = sin(theta * i); // img
         // printf("(%.2f) + j(%.2f)\n", roots[i][0], roots[i][1]);
    }
    return roots;
}

// Recursive FFT : converts coefficient representation to point-value representation
double **eval(int poly[], int N)
{
    double **y = new double *[N];
    for (int i = 0; i < N; i++)
    {
        y[i] = new double[2];
    }

    int A_even[16], A_odd[16];

    // find roots
    double **roots = find_complex_roots(N);

    // base case
    if (N == 1)
    {
        y[0][0] = poly[0];
        y[0][1] = 0;
        return y;
    }

    // performing divide & conquer i.e; dividing polynomial into even and odd terms
    for (int i = 0; i < N / 2; i++)
    {
        A_even[i] = poly[i * 2];    // even indexed
        A_odd[i] = poly[i * 2 + 1]; // odd indexed
    }

    // recursive call for even indexes
    double **y_even = eval(A_even, N / 2);
    // recursive call for odd indexes
    double **y_odd = eval(A_odd, N / 2);

    for (int k = 0; k < N / 2; k++)
    {
        // y[i] = y_even[i] + wn^i * y_odd[i]
        double *mult = multiply_complex_num(roots[k][0], roots[k][1], y_odd[k][0], y_odd[k][1]);
        double *add = add_complex_num(y_even[k][0], y_even[k][1], mult[0], mult[1]);

        y[k][0] = add[0];
        y[k][1] = add[1];

        // y[i+N/2] = y_even[i] - wn^i * y_odd[i]
        add = add_complex_num(y_even[k][0], y_even[k][1], (-mult[0]), (-mult[1]));

        y[k + N / 2][0] = add[0];
        y[k + N / 2][1] = add[1];
    }
    return y;
}

// multiplying two polynomials
double **product_polynomial_evaluations(double **poly_1_evaluations, double **poly_2_evaluations, int N)
{
    double **product_poly_evaluations = new double *[N];
    for (int i = 0; i < N; i++)
    {
        product_poly_evaluations[i] = new double[2];
    }

    for (int i = 0; i < N; i++)
    {
        double *poly_mult = multiply_complex_num(poly_1_evaluations[i][0], poly_1_evaluations[i][1], poly_2_evaluations[i][0], poly_2_evaluations[i][1]);

        product_poly_evaluations[i][0] = poly_mult[0];
        product_poly_evaluations[i][1] = poly_mult[1];
    }
    return product_poly_evaluations;
}

// Inverse Recursive FFT : interpolates point-value representation of product polynomial to coefficient representation
double **evalIFFT(double **product_poly_evaluations, int N)
{
    double **y = new double *[N];
    for (int i = 0; i < N; i++)
        y[i] = new double[2];

    double **B_even = new double *[N / 2];
    double **B_odd = new double *[N / 2];
    for (int i = 0; i < N / 2; i++)
    {
        B_even[i] = new double[2];
        B_odd[i] = new double[2];
    }

    // find roots
    double **roots = find_complex_roots(N);

    //
    if (N == 1)
    {
        y[0][0] = product_poly_evaluations[0][0];
        y[0][1] = product_poly_evaluations[0][1];
        return y;
    }

    // using divide & conquer
    for (int i = 0; i < N / 2; i++)
    {
        B_even[i][0] = product_poly_evaluations[i * 2][0]; // even indexed
        B_even[i][1] = product_poly_evaluations[i * 2][1];

        B_odd[i][0] = product_poly_evaluations[i * 2 + 1][0]; // odd indexed
        B_odd[i][1] = product_poly_evaluations[i * 2 + 1][1];
    }

    // recursive call for even index
    double **y_even = evalIFFT((double **)B_even, N / 2);
    // recursive call for odd index
    double **y_odd = evalIFFT((double **)B_odd, N / 2);

    //
    for (int k = 0; k < N / 2; k++)
    {
        // real part
        // y[i] = y_even[i] + wn^i * y_odd[i]
        double *mult = multiply_complex_num(roots[k][0], -roots[k][1], y_odd[k][0], y_odd[k][1]);
        double *add = add_complex_num(y_even[k][0], y_even[k][1], mult[0], mult[1]);

        y[k][0] = (add[0] / 2);
        y[k][1] = (add[1] / 2);

        // img part
        // y[i+N/2] = y_even[i] - wn^i * y_odd[i]
        add = add_complex_num(y_even[k][0], y_even[k][1], (-mult[0]), (-mult[1]));

        y[k + N / 2][0] = (add[0] / 2);
        y[k + N / 2][1] = (add[1] / 2);
    }
    return y;
}

// finding N : computes the smallest power of 2 greater than or equal to deg1 + deg2 + 1
int find_N(int deg1, int deg2)
{
    return pow(2, ceil(log2(deg1 + deg2 + 1)));
}

// naive approach taking O(n^2) time
int *naive_polynomial_multiplication(int poly_1[], int poly_2[])
{
    int *naive_prod = new int[32];

    // initially, initialize naive_prod with zeros
    for (int i = 0; i < 32; i++)
    {
        naive_prod[i] = 0;
    }

    // multiplication of polynomials
    for (int i = 0; i <= deg1; i++)
    {
        for (int j = 0; j <= deg2; j++)
        {
            naive_prod[i + j] += poly_1[i] * poly_2[j];
        }
    }
    return naive_prod;
}

// function that prints the polynomial
void print_poly(int coeff[], int deg)
{

    for (int i = deg; i >= 0; i--)
    {
        if (coeff[i] > 0)
        {
            if (i != deg)
            {
                cout << " + ";
            }

            if (coeff[i] != 1 && i != 0)
            {
                cout << coeff[i] << "x*" << i;
            }
            else
            {
                if (i != 0)
                    cout << "x*" << i;
                else
                    cout << coeff[i];
            }
        }

        else if (coeff[i] < 0)
        {
            cout << " - ";

            if (i != 0)
            {
                if (coeff[i] == -1)
                {
                    cout << "x*" << i;
                }
                else
                {
                    cout << abs(coeff[i]) << "x*" << i;
                }
            }
            else
            {
                cout << abs(coeff[i]);
            }
        }
    }
    cout << endl;
}

// main function
int main()
{

    // taking two polynomials as input
    cout << " Enter degree of the first polynomial : " << endl;
    cin >> deg1;

    cout << "Enter the " << deg1 + 1 << " coefficients of the first polynomial in the increasing order of the degree of the monomials they belong to: " << endl;
    for (int i = 0; i <= deg1; i++)
        cin >> poly_1[i];

    cout << " Enter degree of the second polynomial : " << endl;
    cin >> deg2;

    cout << "Enter the " << deg2 + 1 << " coefficients of the first polynomial in the increasing order of the degree of the monomials they belong to: " << endl;
    for (int i = 0; i <= deg2; i++)
        cin >> poly_2[i];

    cout << "The first polynomial is : " << endl;
    print_poly(poly_1, deg1);

    cout << "The second polynomial is : " << endl;
    print_poly(poly_2, deg2);

    // calling naive approach function
    int *naive_prod = naive_polynomial_multiplication(poly_1, poly_2);
    cout << "The product of the two polynomials obtained via naive polynomial multiplication is: " << endl;
    print_poly(naive_prod, deg1 + deg2);

    // calling function which finds N
    int N = find_N(deg1, deg2);

    // finding complex roots of unity
    find_complex_roots(N);

    // evaluation of poly[] at wn^0, wn^1....wn^(n-1)
    double **poly_1_evaluations = eval(poly_1, N);

    double **poly_2_evaluations = eval(poly_2, N);

    // calling function to multiply polynomials
    double **product_poly_evaluations = product_polynomial_evaluations(poly_1_evaluations, poly_2_evaluations, N);

    // calling inverse FFT function
    double **invFFT = evalIFFT(product_poly_evaluations, N);
    
    cout << "The product of the two polynomials obtained via polynomial multiplication using FFT is: " << endl;

    // storing coefficients of the result after multiplying two polynomials
    int *fft_prod = new int[32];
    for (int i = 0; i < N; i++)
    {
        // real values only
        fft_prod[i] = (int)round(invFFT[i][0]);
    }
    
    print_poly(fft_prod, deg1 + deg2);

    return 0;
}