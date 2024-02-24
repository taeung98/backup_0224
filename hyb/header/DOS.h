#define N_EIGEN 1000 // Number of eigenenergies
#define N_DOS 1000 // Number of energy interval
#define ETA 0.05 // Broadening rate
#define GREEN(e, eigen_e, fermi) ( ETA / ( pow((e) - ((eigen_e) - (fermi)), 2) + (ETA)*(ETA)) )  // Imaginary Green func
//#include <gsl/gsl_integration.h>
//#include <math.h>
//#include <stdio.h>
//#include <string.h>

void save(char*, double*,double*,int);
double state_num(double*, double);
void Gaussian(double*,double*);
double Fermi(double *eigen_e);

double Fermi(double *eigen_e) {
	double fermi, a, b, c;
	double t = 1;	// hopping term
	double count=0;
	double d = 2; // common diff
	double Ni, Nf,N1,N2;
	double xi[N_EIGEN],wi[N_EIGEN];// node & weight for Gaussian integration
	
	Gaussian(xi,wi); // Gaussian node & weight	

	// 실제 사용할 때는 자기가 구한 고유에너지, 페르미 레벨을 집어넣으면 됨
	for(int i=0; i<N_EIGEN; i++) 
		eigen_e[i] = -2.*t*cos(-M_PI+xi[i]) -2.*(t/2)*cos(2*(-M_PI+xi[i]));
	
	a = 3;// set initial fermi value
	b = a; c = 0; // set bisection node to find chemical potential
	
	// find fermi value with opposite signs to (Ni-0.5) 
	do{
		Ni = state_num(eigen_e,b);
		if (Ni < .5){
			b -= d;
			Nf = state_num(eigen_e,b);
			if(Nf > .5){
				//printf("Nf=%f\n",Nf);
				break;
			}
		}
		else{
			b += d;
			Nf = state_num(eigen_e,b);
			if(Nf < .5){
				//printf("Nf=%f\n",Nf);
				break;
			}
		}
	}while(1 == 1);
	
	// bisection method to find fermi value	with chemicla potential for N/2
	c = (a+b)/2;
	while ( fabs(state_num(eigen_e, c) - .5) > 1e-6 ) {
		if ( state_num(eigen_e, c) > .5 ) {
			a = c;
		} else {
			b = c;
			//printf("N<.5\n");
		}
		c = (a + b) / 2;
		//printf("fermi=%f,N=%.16f\n",c,state_num(eigen_e,c));
	}		

	printf("Fermi value for which N=.5 is %lf, N=%f\n",c,state_num(eigen_e,c));
	
	
	return 0;	
}

double state_num(double *eigen_e, double fermi) {
	double e_min = -8, e_max = 8; // 에너지 범위 지정
	double e[N_DOS],dos[N_DOS]; // DOS
	double xi[N_EIGEN],wi[N_EIGEN];
	//Gaussian Quadrature node & wieghts
	gsl_integration_glfixed_table *w;
	w = gsl_integration_glfixed_table_alloc(N_EIGEN);
	for(int i=0;i<N_EIGEN;i++){
		gsl_integration_glfixed_point(0,2*M_PI,i,&xi[i],&wi[i],w);
	}
	gsl_integration_glfixed_table_free(w);
	
	memset(dos, 0, sizeof(dos)); // DOS array 0으로 초기화

	for(int n=0; n<N_DOS; n++) {
		e[n] = e_min + (e_max - e_min) * n / N_DOS; // DOS 구할 e 계산

		for(int i=0; i<N_EIGEN; i++) {
			dos[n] += GREEN(e[n], eigen_e[i], fermi)/M_PI * wi[i]; // DOS 계산
		}

//		printf("%16.6f\t%16.6f\n", e, dos[n]); // 출력해보고 싶으면 주석 지우기 (왼쪽이 에너지, 오른쪽이 DOS)
	}
	
	double sum=0; // # of states(# of electorns)
	for(int i=0;i<N_DOS;i++)
		sum += dos[i]*(e_max-e_min)/N_DOS/(2*M_PI);
	//printf("%0.16f\n",sum);
	
	//save("DOS",e,dos,N_DOS);
	
	return sum;
}

//Gaussian Quadrature node & wieghts
void Gaussian(double*xi,double*wi){
	gsl_integration_glfixed_table *w;
	w = gsl_integration_glfixed_table_alloc(N_EIGEN);
	for(int i=0;i<N_EIGEN;i++){
		gsl_integration_glfixed_point(0,2*M_PI,i,&xi[i],&wi[i],w);
	}
	gsl_integration_glfixed_table_free(w);
}

void save(char *head, double* x, double* y, int LN){
    char fname[1024];
    sprintf(fname, "%s.txt", head);
    FILE* fp = fopen(fname,"w");
    for(int i=0; i<LN; i++){
        fprintf(fp, "%lf\t%lf\n", x[i], y[i] );
    }
    fclose(fp);

}

