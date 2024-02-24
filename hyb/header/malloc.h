double* mdouble(int size){
	return (double*)malloc(sizeof(double)*size);
}

double complex* mdcomplex(int size){
	return (double complex*)malloc(sizeof(double complex)*size);
}

int* mint(int size){
	return (int*)malloc(sizeof(int)*size);
}

float* mfloat(int size){
	return (float*)malloc(sizeof(float)*size);
}

char* mchar(int size){
	return (char*)malloc(sizeof(char)*size);
}

long long int* mI64(int size){
	return (long long int*)malloc(sizeof(long long int)*size);
}
