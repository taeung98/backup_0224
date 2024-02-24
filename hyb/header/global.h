#define Nmax	8192
#define B		128
#define ETA		0.15
#define PI		M_PI
#define L		128
#define Amax	128
#define Wn(n, B)	((2*(n)+1)*M_PI/(B))
#define Wr(n, N) (-8. + (n)*16./(N-1)) // e_min=-8, e_max=8
#define N_EIGEN 1024
#define N_DOS 1024
#define lim		128

int my_rank, my_size, nearest, Nb;
double  ti[L], wi[L];	//Gaussian node, weight
double mu, *wn, *w;//, *nth;		
double complex *Giwn, *Gw, *Dw, *Gt, *Dt, *Dt_re, *Rhyb_i, *Rhyb_f, *Ihyb_i, *Ihyb_f, *An, *Gt_re;
char FILE_h5[512], Hamilton[512];

struct E_params { double n; double m; double h; int near;}; //N: nearest  n,m,h: 2,3,4 hopping term 
struct G_params { int n; double complex *w; struct E_params e_params;};

double bisection(double desired, double* arr_x, double* arr_y, int size){
	int left = 0;
	int right = size - 1;
	while( left<=right){
		int mid = (left+right)/2;
		double mid_value = (arr_y[left]+arr_y[right])/2.0;
		if( fabs(mid_value-desired) < 1e-6 ){
			return arr_x[mid];
		} else if ( mid_value-desired > 1e-6){
			right = mid - 1;
		} else {
			left = mid + 1;
		}
	}
	double closest_left = arr_x[right];
	double closest_right = arr_x[left];
	return (closest_left + closest_right)/2.0;
}

double NNH(double x, void* parmas) //Next Neighbor Hopping 
{
	struct E_params* p = (struct E_params*)parmas;
	double n,m,h;
	n = p->n; m = p->m; h = p->h; 		

	//normalization of bandwidth
	return -.5 * ( cos(x) + n*cos(2*x) + m*cos(3*x) + h*cos(4*x) ) / ( 1 + n + m + h);
}

void pathname(void* params){
	struct G_params* p = (struct G_params*)params;
	int N = p->e_params.near;
	int n = p->e_params.n * 20;
	int m = p->e_params.m * 20 * 20 / n; 
	int h = p->e_params.h * 20 * 20 * 20 / n / m;
	switch(N){
		case 1: 
			sprintf(Hamilton, "N%1d", N);								break;
		case 2:
			sprintf(Hamilton, "N%1d_n%02d", N, n);						break;
		case 3:
			sprintf(Hamilton, "N%1d_n%02d_m%02d", N, n, m);				break;
		case 4:
			sprintf(Hamilton, "N%1d_n%02d_m%02d_h%02d", N, n, m, h);	break;
		default:
			printf("N-Nearest Neighbor is larger than 4\n");	exit(1);
	}
}

// 함수 호출 많아

double Gre(double x, void *p){
	struct G_params* params = (struct G_params*)p;
	struct E_params* e_params  = &(params->e_params);
	int n = (params->n);
	double complex W = *(params->w + n);
	double (*E)(double, void*) = NNH;
		
	return creal( 1./( W +  mu - E(x,e_params)) );
}

double Gim(double x, void *p){
	struct G_params* params = (struct G_params*)p;
	struct E_params* e_params  = &(params->e_params);
	int n = (params->n);
	double complex W = *(params->w + n);
	double (*E)(double, void*) = NNH;
	
	return cimag( 1. / ( W +  mu - E(x,e_params)) );
}

void calculate_GD(double *W, double complex *G, double complex *D, void*p, char A){ 
	struct G_params* params = (struct G_params*)p;
//	struct E_params* e_params = &(params->e_params);
	double gre, gim, abserr; abserr=0;
	double epsabs = 1e-10; double epsrel = 1e-10;	
//	e_params->N = 1;
	
	if(A == 'r'){
		for(int n=0; n<Nmax; n++){
			W[n] = Wr(n,Nmax);
			params->w[n] = W[n]+(1.*I)*ETA;
		}
	}else{
		for(int n=0; n<Nmax; n++){
			W[n] = Wn(n,B);
			params->w[n] = (1.*I)*W[n];	
		}
	}
	
	gsl_integration_workspace *S = gsl_integration_workspace_alloc(N_EIGEN);
	gsl_function F1,F2;
	F1.function = &Gre;
	F1.params = params;
	F2.function = &Gim;
	F2.params = params;
	
	for(int n=0; n<Nmax; n++){
		gre = 0; gim =0;
		params->n = n;
					
		gsl_integration_qag(&F1,-PI,PI,epsabs,epsrel,N_EIGEN,6,S,&gre,&abserr);
		gsl_integration_qag(&F2,-PI,PI,epsabs,epsrel,N_EIGEN,6,S,&gim,&abserr);
		
		G[n] = (gre + (1.*I)*gim)/(2*PI);
		D[n] = params->w[n] + mu - 1./G[n];
	}
	
	gsl_integration_workspace_free(S);
}

void fermiNdos(void *params){
	mu = 0;
	struct G_params* p = (struct G_params*)params;
//	struct E_params* e_params = &(p->e_params);

	double epsabs = 1e-10; double epsrel = 1e-10;
	double DOS[N_DOS];	  memset(DOS, 0, sizeof(DOS));
	double energy[N_EIGEN];
	double sum_DOS[N_DOS-1]; memset(sum_DOS, 0, sizeof(sum_DOS)); // # of states(electorns)
	double dos, abserr;
	
	for(int n=0; n<N_EIGEN; n++)
		energy[n] = Wr(n, N_EIGEN);
			
	gsl_integration_workspace *S = gsl_integration_workspace_alloc(N_DOS);
	gsl_function F;
	F.function = &Gim;
	F.params = p; 
	for(int n=0;n<N_DOS;n++){
		dos = 0; abserr=0;
		p->n = n;
		p->w[n] = energy[n] + (1.*I)*ETA;
		gsl_integration_qag(&F, -PI, PI, epsabs, epsrel, N_EIGEN, 6, S, &dos, &abserr);
		DOS[n] = -dos/(2*PI)/PI; // (1/2PI): because k sum & (-1/PI): kk relation	
	}
	gsl_integration_workspace_free(S);
	
	double dE = (energy[N_EIGEN-1] - energy[0])/(N_EIGEN-1);
	sum_DOS[0] = (DOS[1]+DOS[0])/2.*dE;	
	for(int i=0;i<N_DOS-1-1;i++){
		sum_DOS[i+1] = sum_DOS[i] + (DOS[i+2]+DOS[i+1])/2.*dE;
	}
	
	mu = bisection(0.5, energy, sum_DOS, N_DOS-1);
}

/*
void fermiNdos(void *params){
	struct G_params* p = (struct G_params*)params;
//	struct E_params* e_params = &(p->e_params);
	double epsabs = 1e-10; double epsrel = 1e-10;
	double DOS[N_DOS];	  memset(DOS, 0, sizeof(DOS));
	double energy[N_EIGEN];
	double sum_DOS[N_DOS-1]; memset(sum_DOS, 0, sizeof(sum_DOS)); // # of states(electorns)
	double dos, abserr;
	
	for(int n=0; n<N_EIGEN; n++)
		energy[n] = Wr(n, N_EIGEN);
			
	gsl_integration_workspace *S = gsl_integration_workspace_alloc(N_DOS);
	gsl_function F;
	F.function = &Gim;
	F.params = p; 
	for(int n=0;n<N_DOS;n++){
		dos = 0; abserr=0;
		p->n = n;
		p->w[n] = energy[n] + (1.*I)*ETA;
		gsl_integration_qag(&F, -PI, PI, epsabs, epsrel, N_EIGEN, 6, S, &dos, &abserr);
		DOS[n] = -dos/(2*PI)/PI; // (1/2PI): because k sum & (-1/PI): kk relation	
	}
	gsl_integration_workspace_free(S);
	
	double dE = (energy[N_EIGEN-1] - energy[0])/(N_EIGEN-1);
	sum_DOS[0] = (DOS[1]+DOS[0])/2.*dE;	
	for(int i=0;i<N_DOS-1-1;i++){
		sum_DOS[i+1] = sum_DOS[i] + (DOS[i+2]+DOS[i+1])/2.*dE;
	}
	
	mu = bisection(0.5, energy, sum_DOS, N_DOS-1);
}
*/

void hyb(double complex *D, double *bath, char domain){
	if( domain == 'r' ){
		for(int n = 0; n < Nmax; n++){
			D[n] = 0.;
			for(int i = 0 ; i<Nb; i++){
				D[n] += bath[i]*bath[i]/(w[n] + (1.*I)*ETA - bath[i+Nb]);
			}
		}
	}else{
		for(int n = 0; n < Nmax; n++){
			D[n] = 0.;
			for(int i = 0; i < Nb; i++){
				D[n] += bath[i]*bath[i]/((1.*I)*wn[n] - bath[i+Nb]);
			}
		}
	}
}

void AN(double complex *func, double* x, double complex* A, double *ti, double *wi, char* fname){
// calculate Legendre coeff using Gaussian quad
    double (*p)(int, double);
    p = P;

    for(int i=0; i<Amax; i++){
        double complex c=0;
        x[i] = (double)i;
        for(int j=0; j<L; j++){
            c += func[j]*p(i, 2*ti[j]/B-1)*wi[j];
        }
        A[i] = c*(2*i+1)/B;
    }
	char path[50]; sprintf(path, "%s/%s", Hamilton, "An");
	//save_An(x, A, Amax, path, fname);
}

void LT( double *ti, double complex *A, double complex *Y, int num){
// Legendre transformation
    double (*p)(int, double);
    p = P;
	double complex f;
    for(int i=0;i<L; i++){
        f = 0;
        for(int j=0; j<num; j++){
            f += A[j]*p(j, 2*ti[i]/B-1);
        }
        Y[i] = f;
    }
}

void FT(double complex* b, double* node1, double complex* af, double* node2, double cycle){
// Fourier transformation using Gaussian quad
	int b_point_num = Nmax;
    int af_point_num = L;
    double complex f =0;
    for(int i=0; i<af_point_num; i++){  //tau index
        f = 0.;
        for(int j=0; j<b_point_num; j++){   // Wn index
            f += (cos(node1[j]*node2[i]) - (1.*I) * sin(node1[j]*node2[i])) * (b[j]-1 / ( (1.*I) * node1[j] ));
        }
        af[i] = f/cycle-0.5;
    }
}

void inverse( double* ti, double complex* A, double complex* y, char* folder){
    int num[5]; for(int i=0;i<3;i++){ num[i]=3+2*i; }; num[3]=30; num[4] = Amax;
	char path[50]="";	snprintf(path, sizeof(path), "%s/%s", Hamilton, folder);
    for(int i=0;i<5;i++){
        char fname[1024];	sprintf(fname, "%d", num[i]);
        LT(ti, A, y, num[i]);
        //save_re(ti, y, L, path, fname);
    }
}

double integral_imag(double *x, double complex *y, int size){
	double result = 0.;
	double dx = (x[size-1] - x[0]) / (size-1);
	for(int i = 1; i < size; i++){
		result += ( cimag(y[i-1]) + cimag(y[i]) ) * dx /2.;
	}
	return result;
}

double abs_sign(double x){
	if(x<0){
		return -1.;
	}
	else{
		return 1.;
	}
}
