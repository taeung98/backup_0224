#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <string.h>
#include <complex.h>
#include <time.h>
#include "header/gsl.h"
#include "header/legendre_poly.h"

#define Lmax 128
#define Nmax 8192
#define nquad 128
#define B 128

void Gaussian(double a,double b, int len, double *ni, double *wi)
{//Gaussian Quadrature node & weights
    gsl_integration_glfixed_table *w;
    w = gsl_integration_glfixed_table_alloc(len);
    for(int i=0;i<len;i++){
        gsl_integration_glfixed_point(a,b,i,&ni[i],&wi[i],w);
    } 
    gsl_integration_glfixed_table_free(w);
}

void L_coeff(double complex* func, double* A, double *ti, double *wi) 
{// calculate Legendre coeff using Gaussian quad
    double (*p)(int, double);
    p = P;

    for(int i=0; i<nquad; i++){
        double complex c=0;
        for(int j=0; j<Lmax; j++){
            c += func[j]*p(i, 2*ti[j]/B-1)*wi[j];
        }
        A[2*i] = creal(c*(2*i+1)/B);
		A[2*i+1] = cimag(c*(2*i+1)/B);
    }
}

void LT( double *ti, double complex *A, double complex *Y, int num)
{// Legendre transformation
    double (*p)(int, double);
    p = P;
	double complex f;
    for(int i=0;i<Lmax; i++){
        f = 0;
        for(int j=0; j<num; j++){
            f += A[j]*p(j, 2*ti[i]/B-1);
        }
        Y[i] = f;
    }
}

void FT(double complex* b, double* node1, double complex* af, double* node2, double cycle)
{// Fourier transformation using Gaussian quad
	int b_point_num = Nmax;
    int af_point_num = Lmax;
    double complex f =0;
    for(int i=0; i<af_point_num; i++){  //tau index
        f = 0.;
        for(int j=0; j<b_point_num; j++){   // Wn index
            f += (cos(node1[j]*node2[i]) - (1.*I) * sin(node1[j]*node2[i])) * (b[j]-1 / ( (1.*I) * node1[j] ));
        }
        af[i] = f/cycle-0.5;
    }
}

double ti[Lmax], wi[Lmax];
int main(int argc, char* argv[])
{
#ifdef DEBUG
	time_t init_time = time(NULL);	
#endif
	Gaussian(0, B, Lmax, ti, wi);

	H5Eset_auto(H5E_DEFAULT, NULL, NULL);
	hid_t file_id, omega_id, omgspace_id;
	hsize_t omega_dims[2];
	H5G_info_t group_info;	
	
	file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT); // argv[1]: file name
	double omega[2][Nmax];
	omega_id = H5Dopen2(file_id, "/omega", H5P_DEFAULT);
	omgspace_id = H5Dget_space(omega_id);
	H5Sget_simple_extent_dims(omgspace_id, omega_dims, NULL);
	H5Dread(omega_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL, H5P_DEFAULT, omega);

	H5Gget_info(file_id, &group_info);
	for(int i = 0; i < group_info.nlinks; i++){
		hid_t group_id, Gspace_id, Dspace_id, Giwn_id, Diwn_id;
		hsize_t dims[2];
		double Gdata[2][Nmax], Ddata[2][Nmax];
		char group_name[256];
		double complex Giwn[Nmax], Diwn[Nmax];
	
		H5Lget_name_by_idx(file_id, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, group_name, sizeof(group_name), H5P_DEFAULT);
		group_id = H5Gopen2(file_id, group_name, H5P_DEFAULT);
			
		Giwn_id = H5Dopen2(group_id, "Giwn", H5P_DEFAULT);		Diwn_id = H5Dopen2(group_id, "Diwn", H5P_DEFAULT);
        Gspace_id = H5Dget_space(Giwn_id);				        Dspace_id = H5Dget_space(Diwn_id);
		H5Sget_simple_extent_dims(Gspace_id, dims, NULL);		H5Sget_simple_extent_dims(Dspace_id, dims, NULL);
		
		H5Dread(Giwn_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL, H5P_DEFAULT, Gdata);
		H5Dread(Diwn_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL, H5P_DEFAULT, Ddata);
		for(int j = 0; j < Nmax; j++){
			Giwn[j] = Gdata[0][j] + 1i*Gdata[1][j];
			Diwn[j] = Ddata[0][j] + 1i*Ddata[1][j];
		}
			
		double complex Gt[Lmax], Dt[Lmax];
		double G_coeff[2*Lmax], D_coeff[2*Lmax];
		hsize_t coeff_dims[2] = {1, 2*Lmax};
		hid_t Gcff_id, Dcff_id, Gcffspace_id, Dcffspace_id;	

		FT(Giwn, omega[1], Gt, ti, B);	FT(Diwn, omega[1], Dt, ti, B);
		L_coeff(Gt, G_coeff, ti, wi);	L_coeff(Dt, D_coeff, ti, wi);

//		printf("%f, %f\n", creal(Dt[0]),cimag(Dt[0]));
		printf("%f, %f\n", creal(Gt[0]),cimag(Gt[0])); 

		H5Gunlink(group_id, "Gcoeff"); 		H5Gunlink(group_id, "Dcoeff");

		Gcffspace_id = H5Screate_simple(2, coeff_dims, NULL);	Dcffspace_id = H5Screate_simple(2, coeff_dims, NULL);
		Gcff_id = H5Dcreate2(group_id, "Gcoeff", H5T_NATIVE_DOUBLE, Gcffspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		Dcff_id = H5Dcreate2(group_id, "Dcoeff", H5T_NATIVE_DOUBLE, Dcffspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		H5Dwrite(Gcff_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, G_coeff);
		H5Dwrite(Dcff_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, D_coeff); 

		H5Dclose(Gcff_id);			H5Dclose(Dcff_id);
		H5Sclose(Gcffspace_id);		H5Sclose(Dcffspace_id);
		H5Sclose(Gspace_id);		H5Sclose(Dspace_id);
		H5Dclose(Giwn_id);			H5Dclose(Diwn_id);
		H5Gclose(group_id);

		printf("%d: sccuess\n", i);
	}
	
	H5Sclose(omgspace_id);
	H5Dclose(omega_id);
	H5Fclose(file_id);
#ifdef DEBUG
	time_t finl_time = time(NULL);
	printf("%f\n", difftime(finl_time, init_time) );
#endif
	return 0;
}
