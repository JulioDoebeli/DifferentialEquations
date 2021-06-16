/*
 * main.c
 * 
 * Copyright 2021 jcld14 <jcld14@inf.ufpr.br>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <float.h>		//Contains defined constants specifying the implementation-specific properties of the floating-point library, such as the minimum difference between two different floating-point numbers (_EPSILON), the maximum number of digits of accuracy (_DIG) and the range of numbers which can be represented (_MIN, _MAX).
#include <math.h>		//For computing common mathematical functions -- see Further math or C++ Programming/Code/Standard C Library/Math for details.
#include <stdio.h>		//Provides the core input and output capabilities of the C language. This file includes the venerable printf function.
#include <stdlib.h>		//For performing a variety of operations, including conversion, pseudo-random numbers, memory allocation, process control, environment, signalling, searching, and sorting.
#include <string.h>		//For manipulating several kinds of strings.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
    #define M_PI_2 3.14159265358979323846/2
#endif

/*  Retorna tempo em milisegundos
    Forma de uso:
 
    real_t tempo1;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo2 = timestamp() - tempo1;
*/
double timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

typedef double real_t;
//Matriz tridiagonal
typedef struct {
	char name;
	//numero de pontos internos na malha
    int n; 
    int m;
    //intervalo
    real_t x0;
    real_t xn;
    //condições de contorno
    real_t y0;
    real_t yn;
    //derivadas 
    real_t (* p)(real_t);
    real_t (* q)(real_t, int);
    real_t (* r)(real_t, int);
}Edo;

typedef struct OrdDiffEq_s{
    char name;
	//numero de pontos internos na malha
    int n; 
    int m;
    //intervalo
    real_t x0;
    real_t xn;
    //condições de contorno
    real_t y0;
    real_t yn;
    //derivadas 
    real_t (* p)(real_t);
    real_t (* q)(real_t);
    real_t (* r)(real_t);
}OrdDiffEq_t;


typedef struct {
    char name;
	//numero de pontos internos na malha
    int n;
    int m;
    //intervalo
    real_t Lx;
    real_t Ly;
    //(0,y)
    real_t (* u1)(real_t, int);
    //(lx, y)
    real_t (* u2)(real_t, int);
    //(x, 0)
    real_t (* u3)(real_t, int);
    //(x, ly)
    real_t (* u4)(real_t, int);
    //
    real_t (* func)(real_t, real_t, int);
}Edo2;

typedef struct ParcDiffEq_s{
    char name;
	//numero de pontos internos na malha
    int n;
    int m;
    //intervalo
    real_t Lx;
    real_t Ly;
    //(0,y)
    real_t (* u1)(real_t);
    //(lx, y)
    real_t (* u2)(real_t);
    //(x, 0)
    real_t (* u3)(real_t);
    //(x, ly)
    real_t (* u4)(real_t);
    //
    real_t (* func)(real_t, real_t);
}ParcDiffEq_t;


real_t Ap(real_t x) { return 0; }
real_t Aq(real_t x) { return 0; }
real_t Ar(real_t x) { return 6*x - 0.5*x*x; }


real_t Bfunc(real_t x, real_t y) { return sin(x)*sin(x); }
real_t Bu1  (          real_t y) { return 20.0; }
real_t Bu2  (          real_t y) { return 45.0; }
real_t Bu3  (real_t x          ) { return 0; }
real_t Bu4  (real_t x          ) { return 100; }

real_t Cp(real_t x) { return 0; }
real_t Cq(real_t x) { return 1; }
real_t Cr(real_t x) { return 0; }


real_t Dfunc(real_t x, real_t y) { return -cos(x+y) -cos(x-y); }
real_t Du1  (          real_t y) { return cos(y); }
real_t Du2  (          real_t y) { return -cos(y); }
real_t Du3  (real_t x          ) { return cos(x); }
real_t Du4  (real_t x          ) { return 0; }



real_t func(real_t x, real_t y, int isFirst){//Eq a.
    //sin²(x)
    if(isFirst) return sin(x)*sin(x);
    else return -cos(x+y) -cos(x-y);
}

real_t u1(real_t var, int isFirst){
    if(isFirst) return 20.0;
    else return cos(var);
}

real_t u2(real_t var, int isFirst){
    if(isFirst) return 45.0;
    else return -cos(var);
}

real_t u3(real_t var, int isFirst){
    if(isFirst) return 0;
    else return cos(var);
}

real_t u4(real_t var, int isFirst){
    if(isFirst) return 100;
    else return 0;
}

void printn(real_t *v, int n){
    for(int i=0; i<n; i++)
        printf("%.7g ", v[i]);
    printf("\n");
}

void printnm(real_t **v, int n, int m){
    for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			printf("%.7g ", v[i][j]);
    printf("\n");
}

real_t gaussSeidel(OrdDiffEq_t *edoeq, real_t *Y) {
    real_t time1 = timestamp();
    char name = edoeq->name;
    //intervalo
    real_t x0 = edoeq->x0;
    real_t xn = edoeq->xn;
    real_t y0 = edoeq->y0;
    real_t yn = edoeq->yn;//condições de contorno
    real_t (* p)(real_t) = edoeq->p;
    real_t (* q)(real_t) = edoeq->q;
    real_t (* r)(real_t) = edoeq->r;
    
    int n = edoeq->n, k, i;
    real_t h, xi, bi, d, di, ds;
    real_t re[n], norma=0;
    real_t biV[n], diV[n], dV[n], dsV[n];
    
    //Largura do passo da malha
    h = (xn - x0)/(n+1);
    for(k=0; k<50; ++k){
        for(i=0; i<n; ++i){
			//valor xi da malha
            xi = x0 + (i+1)*h;
            //termo independente
            bi = h*h * r(xi);
            biV[i] = bi;
            di = 1 - h*p(xi)/2.0;//diagonal inferior
            d = -2 + h*h * q(xi);//diagonal principal
            ds = 1 + h*p(xi)/2.0;//diagonal superior
            if(i==0)          bi -= ds*Y[i+1] + y0 * (1 - h*p(x0+h)/2.0);
            else if(i == n-1) bi -= di*Y[i-1] + yn * (1 + h*p(xn-h)/2.0);
            else              bi -= ds*Y[i+1] + di*Y[i-1];

            Y[i] = bi/d; //calcula incognita
            diV[i] = di;
            dV[i] = d;
            dsV[i] = ds;
        }
    }
    for(int i=0; i<n; i++){//fazendo residuo
        if(i==0) re[i] = biV[i] - (dsV[i]*Y[i+1] + dV[i]*Y[i] + y0 * (1 - h*p(x0+h)/2.0));
        else if(i == n-1)  re[i] = biV[i] - (diV[i]*Y[i-1] + dV[i]*Y[i] + yn * (1 + h*p(xn-h)/2.0));
        else re[i] = biV[i] - (diV[i]*Y[i-1] + dV[i]*Y[i] + dsV[i]*Y[i+1]);
    }

    for(int i=0; i<n; i++)//fazendo norma
        norma += re[i]*re[i];
    norma = sqrt(norma);


    printf("***** item (%c): n = %d, H = %.7g\n",name, n, h);
    printf("SL:\n");
    
    printn(dsV, n);
    printn(dV, n);
    printn(diV, n);
    printn(biV, n);
    printf("Y: ");
    printn(Y, n);
    
    /*printf("Tridiagonal matrix:\n"); 
    printf("Upper Diagonal: %.7g\n", dsV[0]);
    printf("Main Diagonal:  %.7g\n", dV[0]);
    printf("Lower Diagonal: %.7g\n", diV[0]);
    printf("b: ");
	printf("{");
	printf("%.7g ", biV[0]);
	for (int i=1; i<n; i++) {
		printf(", %.7g", biV[i]);
	}
    printf("}\n");
    printf("{");
	printf("%.7g ", biV[0]);
	for (int i=1; i<n; i++) {
		printf(", %.7g", Y[i]);
	}
    printf("}\n");*/
    
    
    real_t time2 = timestamp();
    printf("Norma L2: %.7g, Tempo: %.7g ms\n\n\n", norma, time2-time1);
    
    return norma;
}




void gaussSeidel2(ParcDiffEq_t *edoeq, int isFirst) {
    real_t time1 = timestamp();
    char name = edoeq->name;
    int n = edoeq->n;
    int m = edoeq->m;
    real_t (* func)(real_t x, real_t y) = edoeq->func;
    real_t (*   u1)(          real_t y) = edoeq->u1;
    real_t (*   u2)(          real_t y) = edoeq->u2;
    real_t (*   u3)(real_t x          ) = edoeq->u3;
    real_t (*   u4)(real_t x          ) = edoeq->u4;
    real_t U[n][m];
    real_t r[n][m];
    real_t norma=0;
    for(int j=0; j<edoeq->m; j++)//zerando a matriz
        for(int i=0; i<edoeq->n; i++)
            U[i][j] = 0;
    
    real_t hx = 0;
	real_t hy = 0;
	real_t xi = 0;
	real_t bi = 0;
	real_t yj = 0;
	real_t d = 0;
	real_t di = 0;
	real_t ds = 0;
	real_t di2 = 0;
	real_t ds2 = 0;
	real_t biV[n][m];
	
    hx = edoeq->Lx/(n+1);
    hy = edoeq->Ly/(m+1);//largura do passo
    
    for(int k=0; k<50; ++k){
        
        for(int j=0; j<m; j++){
            
            for(int i=0; i<n; i++){
                xi = (i+1)*hx;//valor xi da malha
                yj = (j+1)*hy;//valor yj da malha

                bi = hx*hx*hy*hy*func(xi, yj);//termo indep.
                biV[i][j] = bi;
                di = hy*hy;
                di2 = hx*hx;
                d = -2.0*(hx*hx + hy*hy) - u3(0);
                ds = hy*hy;
                ds2 = hx*hx;


                if( i==0 ) {
					if((j == 0)) {
						bi -= di2*u3(xi) + di*u1(yj) + ds*U[i+1][j] + ds2*U[i][j+1];
					} else if((j != m-1)) {
						bi -= di2*U[i][j-1] + di*u1(yj) + ds*U[i+1][j] + ds2*U[i][j+1];
					} else if((j == m-1)) {
						bi -= di2*U[i][j-1] + di*u1(yj) + ds*U[i+1][j] + ds2*u4(xi);
					}					
				} else if (i==n-1) {
					if((j == 0)) {
						bi -= di2*u3(xi) + di*U[i-1][j] + ds*u2(yj) + ds2*U[i][j+1];
					} else if((j != m-1)) {
						bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*u2(yj) + ds2*U[i][j+1];
					} else if((j == m-1)) {
						bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*u2(yj) + ds2*u4(xi);
					}
					
				} else if (i!=n-1) {
					if((j == 0)) {
						bi -= di2*u3(xi) + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1];
					} else if((j != m-1)) {
						bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1];
					} else if((j == m-1)) {
						bi -= di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] +  ds2*u4(xi); 
					}
				}
				U[i][j] = bi/d; //calcula incognita
            }

        }
        

    
    }

    for(int j=0; j<m; j++){//fazendo residuo
        for(int i=0; i<n; i++){
            xi = (i+1)*hx;//valor xi da malha
            yj = (j+1)*hy;//valor yj da malha
            if((i == 0) && (j == 0)) {
				r[i][j] = biV[i][j] - (di2*u3(xi) + di*u1(yj) + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            } else if((i == 0) && (j != m-1)) {
				r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*u1(yj) + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            } else if((i == 0) && (j == m-1)) {
				r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*u1(yj) + ds*U[i+1][j] + ds2*u4(xi) + d*U[i][j]);
            }
            else if((i == n-1) && (j == 0))
				r[i][j] = biV[i][j] - (di2*u3(xi) + di*U[i-1][j] + ds*u2(yj) + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == n-1) && (j != m-1))
				r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*u2(yj) + ds2*U[i][j+1] + d*U[i][j]);
            else if((i == n-1) && (j == m-1))
				r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*u2(yj) + ds2*u4(xi) + d*U[i][j]);
            
            else if((i != n-1) && (j == 0))
				r[i][j] = biV[i][j] - (di2*u3(xi) + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i != n-1) && (j != m-1))
				r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] + ds2*U[i][j+1] + d*U[i][j]);
            else if((i != n-1) && (j == m-1))
				r[i][j] = biV[i][j] - (di2*U[i][j-1] + di*U[i-1][j] + ds*U[i+1][j] +  ds2*u4(xi) + d*U[i][j]);
        }
    }
    for(int j=0; j<m; j++) //fazendo norma
        for(int i=0; i<n; i++)
            norma+=r[i][j]*r[i][j];
    norma = sqrt(norma);
    
	printf("***** (%c): L = %.4f, W = %.4f, n = %d, m = %d, Hx = %.4f , Hy = %.4f\n", name, edoeq->Lx, edoeq->Ly, n, m, hx, hy);
	
	
	for(int i=0; i<m*n - 2; i++)
        printf("%f ", ds2);
    printf("\n");
    
    for(int i=0; i<m*n - 1; i++)
        printf("%f ", ds);
    printf("\n");
    
    for(int i=0; i<m*n; i++)
        printf("%f ", d);
    printf("\n");

    for(int i=0; i<m*n-1; i++)
        printf("%f ", di);
    printf("\n");

    for(int i=0; i<m*n-2; i++)
        printf("%f ", di2);
    printf("\n");

    for(int j=0; j<m; j++)
        for(int i=0; i<n; i++)
            printf("%f ", biV[i][j]);
    printf("\n");
    
    printf("\nT: ");
    for(int j=0; j<edoeq->m; j++)
        for(int i=0; i<edoeq->n; i++)
            printf("%.7g ", U[i][j]);
    printf("\n");
    
    real_t time2 = timestamp();
    printf("Norma L2: %.7g, Tempo: %.7g ms\n\n\n", norma, time2-time1);
    /*
    printf("Pentadiagonal Matrix (%dx%d)\n",n,n); 
	printf("Upper Diagonal 2: %f\n", ds2);
	printf("Upper Diagonal 1: %f\n", ds);
	printf("Main  Diagonal :  %f\n", d);
	printf("Lower Diagonal 1: %f\n", di);
	printf("Lower Diagonal 2: %f\n", di2);
    printf("b:");
    printf("{\n   ");
	for(int j=0; j<m; j++) {
		printf("{");
		printf("%.7g", biV[0][j]);
        for(int i=1; i<n; i++) {
            printf(", %.7g", biV[i][j]);
		}
		if (j+1!=m)
			printf("},\n   ");
		else
		    printf("}\n  ");
	}
    printf("}\n");
    printf("T:");
    printf("{\n   ");
	for(int j=0; j<m; j++) {
		printf("{");
		printf("%.7g", biV[0][j]);
        for(int i=1; i<n; i++) {
            printf(", %.7g", U[i][j]);
		}
		if (j+1!=m)
			printf("},\n   ");
		else
		    printf("}\n  ");
	}
    printf("}\n");
    printf("Norma L2: %.7g\n", norma);*/
}

   

int main(){
    //OrdDiffEq_t eq;
    OrdDiffEq_t a;
    OrdDiffEq_t c;
    a.name = 'a';	a.x0 = 0;	a.xn = 12;	a.n = 5;	a.m = 0;	a.y0 = 0;	a.yn = 0;	a.p = Ap;	a.q = Aq;	a.r = Ar;
    c.name = 'c';	c.x0 = 0;	c.xn = 1;	c.n = 5;	a.m = 0;	c.y0 = 0;	c.yn = 1;	c.p = Cp;	c.q = Cq;	c.r = Cr;
    
    ParcDiffEq_t b;
	ParcDiffEq_t d;
	b.name = 'b';	b.Lx = 6;		b.Ly = 8;		b.n = 5;	b.m = 3;	b.func = Bfunc;	b.u1 = Bu1;	b.u2 = Bu2;	b.u3 = Bu3;	b.u4 = Bu4;
	d.name = 'd';	d.Lx = M_PI;	d.Ly = M_PI_2;	d.n = 5;	d.m = 3;	d.func = Dfunc;	d.u1 = Du1;	d.u2 = Du2;	d.u3 = Du3;	d.u4 = Du4; 
    
    real_t *Y;
    Y = (real_t *)calloc((a.n), sizeof(real_t));
    
    /// Calcula Sistemas
    /// A
    // Calcula A, n=5
    gaussSeidel( &a, Y);
    
    // Calcula A, n=10
    a.n = 10;
    for(int i=0; i<a.n; i++) Y[i] = 0;
    gaussSeidel( &a, Y);
    
    /// B
    // Calcula B, n=5
    gaussSeidel2( &b, 1);
	
    // Calcula B, n=10
    b.n = 10;
    gaussSeidel2( &b, 1);
    
    /// C
    // Calcula C, n=5
    gaussSeidel( &c, Y);
    
    // Calcula C, n=10
    c.n = 10;
    for(int i=0; i<c.n; i++) Y[i] = 0;
    gaussSeidel( &c, Y);
    
    /// D
    // Calcula D, n=5
    gaussSeidel2( &d, 0);
    
    // Calcula D, n=10
    d.n = 10;
    gaussSeidel2( &d, 0);
    
}
