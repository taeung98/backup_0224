#include "header/legendre_poly.h"

//Legendre Transform_Trapezoid
/*
void AN_LT_T(double complex* func, double* x, double complex *A) {
	double (*p)(int, double);
	p = P;

	for(int i=0; i<Amax; i++){
		double dx = B/(L-1);
		double complex c=0;
		x[i] = (double)i;
		for(int j=0; j<L; j++){ 
			c += func[j]*p(i,2*(j*dx)/B-1);               // [0,B]
		//	c += -func[L-j]*p(i,-2*(-B+j*dx)/B-1);        //[-B,0]
		//	c += -func[L-j]*p(i,-1+j*dx/B)+func[j]*p(i,j*dx/B); //[-B,B]
	    }
		A[i] = c*(2*i+1)/B*dx;                      // [0,B],[-B,0]
	//  A[i] = c*(2*i+1)/(2*B)*dx;                  // [-B,B]
	}
}

void LT_T(double complex* func, double* x, double complex *A, double complex *Y) {
	double (*p)(int, double);
	p = P;

	for(int i=0;i<L; i++){
		double complex f = 0;
	//	x[i] = (B/(L-1))*i;                         // [0,B]
	//  x[i] = -B+(B/L)*i;                          // [-B,0]
	//  x[i] = -B+2*(B/L)*i;                        // [-B,B]
		for(int j=0; j<Amax; j++){
		    f += A[j]*p(j,2*x[i]/B-1);              // [0,B]
	//      f += A[j]*p(j,-2*x[i]/B-1);             // [-B,0]
	//      f += A[j]*p(j,x[i]/B);                  // [-B,B]
		}
		Y[i] = f;
	}
}
*/
//Legendre Transformation_Gaussian

void AN(double complex *func,double* x, double complex* A,double *ti,double *wi){
	
	double (*p)(int, double);
	p = P;
	
	for(int i=0; i<Amax; i++){
		double complex c=0;
		x[i] = (double)i;
		for(int j=0; j<L; j++){ 
			c += func[j]*p(i,2*ti[j]/B-1)*wi[j];
        }
	    A[i] = c*(2*i+1)/B;
	}
}

void LT(double complex *func,double *ti,double complex *A,double complex *Y){
	
	double (*p)(int, double);
	p = P;
	
	for(int i=0;i<L; i++){
		double complex f = 0;
		for(int j=0; j<Amax; j++){
			f += A[j]*p(j,2*ti[i]/B-1);
		}
		Y[i] = f;
	}

}
