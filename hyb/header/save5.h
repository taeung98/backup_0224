#ifndef SAVE_H
#include "header/sort.h"
#define SAVE_H

void sort_test(double *v, double *ep, double *V)
{
    int index[Nb];
    for(int i=0; i < Nb; i++){
        ep[i] = v[Nb+i]; // Nb ~ 2*Nb-1
        index[i] = i;
    }
    quickSort(ep, index, 0, Nb-1);
    for(int i = 0; i < Nb; i++)
		V[i] = v[index[i]];  
} 

double test_f(double *v, void *params){
    double V[Nb], ep[Nb], result = 0.; 
    double complex D;
    double complex* Diwn = (double complex*)params;

    sort_test(v, ep, V); 

    for(int n=0;n<lim;n++){
        D = 0.; 
        for(int j=0;j<Nb;j++){
            D += cabs(V[j])*cabs(V[j])/((1.*I)*wn[n]-ep[j]);
        }   
        result += (creal( Diwn[n] - D ) * creal( Diwn[n] - D ) + cimag( Diwn[n] - D ) * cimag( Diwn[n] - D )); // hyb diff
    }   
    return result/lim;
}

void save_freq(double* W, double* WN){
	hid_t file_id, dataset_id, dataspace_id, attribute_id, attr_dataspace_id;
	hsize_t dims[2]; 
	hsize_t attr_dims[1];
	
	dims[0] = 2; dims[1] = Nmax;
	double freq[2][Nmax];	// 0: real_domain, 1: iamg_domain 
	for(int i = 0; i<Nmax; i++){
		freq[0][i] = W[i];
		freq[1][i] = WN[i];
	}
	
	file_id = H5Fopen(FILE_h5, H5F_ACC_RDWR, H5P_DEFAULT); 
	if(file_id < 0)
		file_id	= H5Fcreate(FILE_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	char name[16]; sprintf(name, "/omega");	
	dataspace_id = H5Screate_simple(2, dims, NULL);
	dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
	if(dataset_id < 0){
		dataset_id = H5Dcreate2(file_id, name, H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}else{
		H5Sclose(dataspace_id);
		H5Fclose(file_id);	
		return;
	}
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, freq);	
		
	int attr_data = Nmax; 
	attr_dims[0] = 1;	
	attr_dataspace_id = H5Screate_simple(1, attr_dims, NULL);	
	attribute_id = H5Acreate2(dataset_id, "Nmax", H5T_STD_I32BE, attr_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute_id, H5T_NATIVE_INT, &attr_data);
	
	H5Aclose(attribute_id);
	H5Sclose(attr_dataspace_id);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);	
	H5Fclose(file_id);

	return;
}

void save_data(double complex* data, char* dataset_name){
	hid_t file_id, group1_id, dataset_id, dataspace_id;
	hsize_t func_dims[2]; func_dims[0] = 2; func_dims[1] = Nmax; 

	double func[2][Nmax];//0: real_part, 1: imag_part
	for(int i = 0; i < Nmax; i++){
		func[0][i] = creal(data[i]);
		func[1][i] = cimag(data[i]);
	}
	
	file_id = H5Fopen(FILE_h5, H5F_ACC_RDWR, H5P_DEFAULT); 
	if(file_id < 0)
		file_id = H5Fcreate(FILE_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	char Group1[1024]; 	sprintf(Group1, "/%s", Hamilton); 
	group1_id = H5Gopen2(file_id, Group1, H5P_DEFAULT);
	if( group1_id < 0)
		group1_id = H5Gcreate2(file_id, Group1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	dataspace_id = H5Screate_simple(2, func_dims, NULL);
	char Dname[16]; sprintf(Dname, "%s", dataset_name);
	dataset_id = H5Dopen2(group1_id, Dname, H5P_DEFAULT);
	if(dataset_id < 0){
		dataset_id = H5Dcreate2(group1_id, Dname, H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}else{	
		H5Sclose(dataspace_id);	
		H5Gclose(group1_id);
		H5Fclose(file_id);
		return;
	}
	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, func);	
		
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);	
	H5Gclose(group1_id);
	H5Fclose(file_id);
	
	return;
}

void load_chi(hid_t path_id, int chi_num, double *chi_val)
{
    hid_t dataset_id, attr_id;

	for (int i = 0; i < chi_num; i++) {
        char dataset_name[20];
        sprintf(dataset_name, "P%d", i);
        // Open dataset
        dataset_id = H5Dopen2(path_id, dataset_name, H5P_DEFAULT);
#ifdef DEBUG
		double dataset_data[2][2*Nb];
		hid_t dataspace_id;
		hsize_t dims[2];
        if (dataset_id < 0) {
            fprintf(stderr, "Unable to open dataset: %s\n", dataset_name);
            return;
        }
        // Open dataspace
        dataspace_id = H5Dget_space(dataset_id);
		H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
		// Read the dataset into the array
		if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataset_data) < 0) {
			fprintf(stderr, "Error reading dataset\n");
			H5Sclose(dataspace_id);
			H5Dclose(dataset_id);
			return;
		}
	    printf("Entire dataset:\n");
		for (int j = 0; j < dims[1]; j++) 
			printf("%e\t%e\n ", dataset_data[0][j], dataset_data[1][j]);
#endif
        // Read 'chi' attribute
        attr_id = H5Aopen(dataset_id, "chi", H5P_DEFAULT);
#ifdef DEBUG
        if (attr_id < 0) {
            fprintf(stderr, "Unable to open attribute 'chi' for dataset: %s\n", dataset_name);
            H5Dclose(dataset_id);
            return;
        }
        H5Sclose(dataspace_id);
#endif		

        H5Aread(attr_id, H5T_NATIVE_DOUBLE, &(chi_val[i]));
//		printf("chi = %e, test_f = %e\n", chi_val[i], test_f(dataset_data[1], NULL));
        // Close dataset and dataspace
        H5Aclose(attr_id);
        H5Dclose(dataset_id);
    }
}

void process_datasets(hid_t path_id) 
{
    hid_t dataset_id;
//	herr_t status;
	hsize_t num_objs;	H5Gget_num_objs(path_id, &num_objs);
	int dataset_count = (int)num_objs;	
    double chi_values[dataset_count];

	if(dataset_count == 10) 
		return;
#ifdef DEBUG
	printf("\n accumulate data\n");
	printf("\n dataset_count = %d\n", dataset_count);
#endif
	load_chi(path_id, dataset_count, chi_values);

	int chi_index[dataset_count];
	for(int i = 0; i < dataset_count; i++)
	{
		 chi_index[i] = i;
	//	 printf("load: chi[%d] = %e\n", i ,chi_values[i]);
	}
    // Sort chi_values in descending order
	quickSort(chi_values, chi_index, 0, dataset_count-1);	
    // Open datasets in descending order of chi_values
    for (int i = 0; i < dataset_count; i++)
	{
        char dataset_name[20];
        sprintf(dataset_name, "P%d", chi_index[i]);
        dataset_id = H5Dopen2(path_id, dataset_name, H5P_DEFAULT);
        // Remove datasets with higher chi_values
        if (i < 10) {
            // Rename datasets with lower chi_values
            char new_dataset_name[20];
            sprintf(new_dataset_name, "temp%d", i);
/*status=*/	H5Lmove(path_id, dataset_name, path_id, new_dataset_name, H5P_DEFAULT, H5P_DEFAULT);
/*
			if (status < 0) 
                fprintf(stderr, "Error renaming dataset: %s to %s\n", dataset_name, new_dataset_name);
*/
        }else {
/*status*/	H5Gunlink(path_id, dataset_name);
			H5Gget_num_objs(path_id, &num_objs);
/*
			if (status < 0) 
                fprintf(stderr, "Error removing dataset: %s\n", dataset_name);
*/
        }
		H5Dclose(dataset_id);
    }
	for (int i = 0; i < 10; i++){
		char dataset_name[20];	sprintf(dataset_name, "temp%d", i);
        dataset_id = H5Dopen2(path_id, dataset_name, H5P_DEFAULT);
		char new_dataset_name[20]; sprintf(new_dataset_name, "P%d", i);
/*status*/ H5Lmove(path_id, dataset_name, path_id, new_dataset_name, H5P_DEFAULT, H5P_DEFAULT);
		H5Dclose(dataset_id);
	}	
#ifdef DEBUG	
	printf("\n after sort: \n");	
	load_chi(path_id, dataset_count, chi_values);
#endif
}

void save_bath(double* sdata, double* chi, int snum, int* Iter, int TotIter){
	hid_t file_id, dataset_id, dataspace_id, group1_id, group2_id, chi_id, chi_dataspace_id, Iter_id, Iterspace_id, TotIter_id, TotIterspace_id;
	hsize_t dims[2], chi_dims[1]; chi_dims[0] = 1;	
	int	Iter_data;
	double chi_data,  bath[2][2*Nb];
	file_id = H5Fopen(FILE_h5, H5F_ACC_RDWR, H5P_DEFAULT);
	if(file_id < 0){
		file_id = H5Fcreate(FILE_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	}	
	dims[0] = 2; dims[1] = 2*Nb;
	
	char Group1[1024], Group2[1024], Data_num[128]; 
	sprintf(Group1, "/%s", Hamilton); sprintf(Group2, "%s%d", "bath", Nb);

	group1_id = H5Gopen2(file_id, Group1, H5P_DEFAULT);
	if( group1_id < 0)
		group1_id = H5Gcreate2(file_id, Group1, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


	group2_id = H5Gopen2(group1_id, Group2, H5P_DEFAULT);	
	if( group2_id < 0){
		group2_id = H5Gcreate2(group1_id, Group2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
	}

	TotIter_id = H5Aopen(group2_id, "Tot_Iter", H5P_DEFAULT);
	if ( TotIter_id < 0 ) {
		TotIterspace_id = H5Screate_simple(1, chi_dims, NULL);
		TotIter_id = H5Acreate2(group2_id, "Tot_Iter", H5T_NATIVE_INT, TotIterspace_id, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(TotIter_id, H5T_NATIVE_INT, &TotIter);
		H5Sclose(TotIterspace_id);
	}else {
		int TotIter_up;
		H5Aread(TotIter_id, H5T_NATIVE_INT, &TotIter_up);
		TotIter_up += TotIter;
		H5Awrite(TotIter_id, H5T_NATIVE_INT, &TotIter_up);
	}
	H5Aclose(TotIter_id);

	hsize_t num_objs; H5Gget_num_objs( group2_id, &num_objs);

	int snum1 = (int)num_objs; 
//	printf("snum1 = %d\n", snum1);
	int snum2 = snum1 + snum;
		
	for(int i = snum1; i < snum2; i++)
	{
		int k = i - snum1;
//		printf("k = %d\n", k);
		for(int j = 0; j < 2*Nb; j++)
		{
			bath[0][j] = sdata[k*4*Nb + j];
			bath[1][j] = sdata[k*4*Nb + 2*Nb + j];
//			printf("%e\t%e\n", bath[0][j], bath[1][j]);
		}

//		printf("chi = %e, test_f = %e\n\n", chi[k], test_f(&sdata[k*4*Nb + 2*Nb], NULL) );

		sprintf(Data_num,"%c%d",'P',i);
		dataspace_id = H5Screate_simple(2, dims, NULL);
		dataset_id = H5Dcreate2(group2_id, Data_num, H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, bath);	
		
		chi_data = chi[k];	chi_dims[0] = 1;	
		Iter_data = Iter[k];
		//printf("%d\n", Iter_data);
		chi_dataspace_id = H5Screate_simple(1, chi_dims, NULL);	
		Iterspace_id = H5Screate_simple(1, chi_dims, NULL);	
	
		chi_id = H5Acreate2(dataset_id, "chi", H5T_IEEE_F64BE, chi_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
		Iter_id = H5Acreate2(dataset_id, "Iter", H5T_NATIVE_INT, Iterspace_id, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(chi_id, H5T_NATIVE_DOUBLE, &chi_data);
		H5Awrite(Iter_id, H5T_NATIVE_INT, &Iter_data);
		
		H5Aclose(Iter_id);			H5Sclose(Iterspace_id);
		H5Sclose(chi_dataspace_id);	H5Aclose(chi_id);	
		H5Dclose(dataset_id);
	}
	process_datasets(group2_id);

	H5Gclose(group2_id);
	H5Gclose(group1_id);
	H5Fclose(file_id);
}

int compare_chi(const void *a, const void *b) {
	double chi_a = *((double *)a);
	double chi_b = *((double *)b);

	if (chi_a < chi_b) return -1;
	if (chi_a > chi_b) return 1;
	return 0;
}

#endif
