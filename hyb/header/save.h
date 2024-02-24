void save_hyb(double *x, double complex *y, int LN, char* fname){
/*	char fname[512];
	snprintf(fname,sizeof(fname), "%s/%s.txt",path,name);
*/
	FILE *fp = fopen(fname, "w");
	if(fp != NULL){
		for(int i=0; i<LN; i++){
			fprintf(fp, "%lf\t%.16lf\t%.16lf\n", x[i], creal(y[i]),cimag(y[i]) );
		}
		fclose(fp);
	} else{
		printf("Unable to create or save the data at '%s'.\n", fname);
	}
}

/*void make_folder(char* path) {
//    printf("input path: %s\n",path);
	const mode_t mode = 0777;
    char* path_copy = strdup(path);  // 입력 경로를 복사하여 사용

    if (path_copy == NULL) {
        perror("Memory allocation failed");
        return;
    }

    char* token = strtok(path_copy, "/");
    char current_path[1024] = "";
    while (token != NULL) {
		strcat(current_path, token);
        // 폴더 생성 시 이미 존재하는 경우를 체크하여 처리
        int status = mkdir(current_path, mode);

        if (status != 0) {
            if (errno == EEXIST) {
//                printf("Folder '%s' already exists.\n", current_path);
            } else {
//                perror("Failed to create folder");
            }
        }

        strcat(current_path, "/");
        token = strtok(NULL, "/");
    }

    free(path_copy);
}

// 복소수 함수 txt로 저장
void save_data(double *wn, double complex *arr, int LN,char *path,char* head){
	make_folder(path);
	char fname[512];
    char filepath[512];
	sprintf(fname, "%s/%s.txt",path, head);
    if( access(fname,F_OK) != -1 ){
	}else{
		FILE *fp = fopen(fname, "w");
		if (fp == NULL){
			perror("Error opening file");
			exit(EXIT_FAILURE);
		}
		for(int i=0; i<LN; i++){
			fprintf(fp, "%lf\t%.16lf\t%.16lf\n", wn[i], creal(arr[i]), cimag(arr[i]) );
		}
		fclose(fp);
		printf("%s created and data saved.\n",fname);
	}
}

void save_An(double *wn, double complex *arr, int LN,char *path,char* head){
	make_folder(path);	
	char fname[512];
	sprintf(fname, "%s/%s.txt",path, head);
	if( access(fname,F_OK) != -1){
	}else{
		FILE *fp = fopen(fname, "w");
		if(fp == NULL){
			perror("Error opening file");
			exit(EXIT_FAILURE);
		}
		for(int i=0; i<LN; i++){
			fprintf(fp, "%lf\t%.16lf\t%.16lf\n", wn[i], creal(arr[i]), cimag(arr[i]) );
		}
		fclose(fp);
		printf("%s created and data saved.\n",fname);
	}
}

void save_re(double *wn, double complex *arr, int LN,char *path,char* head){
	make_folder(path);	
	char fname[512];
	sprintf(fname, "%s/%s.txt",path, head);
	if( access(fname,F_OK) != -1){
	}else{
		FILE *fp = fopen(fname, "w");
		if(fp == NULL){
			perror("Error opening file");
			exit(EXIT_FAILURE);
		}
		for(int i=0; i<LN; i++){
			fprintf(fp, "%lf\t%.16lf\t%.16lf\n", wn[i], creal(arr[i]), cimag(arr[i]) );
		}
		fclose(fp);
		printf("%s created and data saved.\n",fname);
	}
}

void save(double* x, double* y, int LN, char *path, char* head){
	char fname[1024];
	sprintf(fname, "%s/%s.txt", path, head);
	if( access(fname,F_OK) != -1){
	}else{
		FILE *fp = fopen(fname, "w");
		if(fp == NULL){
			perror("Error opening file");
			exit(EXIT_FAILURE);
		}
		for(int i=0; i<LN; i++){
			fprintf(fp, "%0.16lf\t%0.16lf\n", x[i], y[i] );
		}
		fclose(fp);
		printf("%s created and data saved.\n",fname);
	}
}

//소수점 16자리까지 비교
void save_err(double *wn, double complex *arr, int LN, char *path,char* head){
	make_folder(path); 
	char fname[1024];
    sprintf(fname, "%s/%s.txt",path,head);
    FILE *fp = fopen(fname, "w");
    for(int i=0; i<LN; i++){
        fprintf(fp, "%lf\t%.16lf\t%.16lf\n", wn[i], fabs(creal(arr[i])),-fabs(cimag(arr[i])) );
    }
    fclose(fp);
}

void save_bath(double* bath, int Nb, char* path,char* name){
	make_folder(path);
	char fname[1024];	
	snprintf(fname,sizeof(fname), "%s/%s.txt",path,name);
	FILE *fp = fopen(fname,"w");
	if(fp != NULL){
		for(int i=0; i<Nb; i++){
			fprintf(fp, "%0.16lf\t%0.16lf\n", bath[Nb+i], bath[i] );
		}
		fclose(fp);
	} else{
		printf("Unable to create or save the data at '%s'.\n", fname);
	}
}


//	 간단하게 비교

	save_hyb(wn, Diwn, Nmax, "o.txt"); 
	hyb(Ihyb_i, bath_init, 'i'); save_hyb(wn, Ihyb_i, Nmax, "b.txt");
	hyb(Ihyb_f, bath, 'i');		save_hyb(wn ,Ihyb_f, Nmax, "f.txt");
	/
	save_hyb(w, Dw, Nmax, "O.txt");
	hyb(Rhyb_i, bath_init, 'r'); save_hyb(w, Rhyb_i, Nmax, "B.txt");
	hyb(Rhyb_f, bath, 'r');		save_hyb(w ,Rhyb_f, Nmax, "F.txt");
	printf("f = %e\n", *chi);
*/
/*
void save_add(char* path, char* filename, void* column1,void* column2,int s){
	char filepath[256];
	sprintf(filepath, "%s/%s.txt",path,filename);

	FILE* file = fopen(filepath, "a");
	if( file == NULL){
		printf("Error: Unable to open the file %s\n",filepath);	
	}
	if(s){
		char* c1 = (char*)column1;
		double* c2 = (double*)column2;
		fprintf(file,"%16s\t%.16lf\n",c1, *c2);
	}else{
		double* c1 = (double*)column1;
		double* c2 = (double*)column2;
		fprintf(file,"%.16lf\t%.16lf\n",*c1, *c2);
	}
	fclose(file);
}

// 특정 시간에 대한 주파수 함수
void save_d(char *head,double *wn, double complex *arr, int LN, int v, char* subfolder){
    char fname[1024];
    sprintf(fname, "%s/%s%d.txt", subfolder, head,v);
    FILE *fp = fopen(fname, "w");
    for(int i=0; i<LN; i++){
        fprintf(fp, "%lf\t%lf\t%lf\n", wn[i], creal(arr[i]), cimag(arr[i]) );
    }
    fclose(fp);
}
int out_data(char *head, double complex *y, int LN, char* subfolder) {
    char *index[2];
	index[0] = "_Real.txt"; index[1] = "_Imag.txt";
    char data[32*LN];
    FILE *fp = NULL;
	for(int i=0;i<LN;i++)
		y[i] = 0;
	for(int i=0;i<2;i++){
		char name[1024]="";
		sprintf(name,"%s/%s%s", subfolder, head , index[i]);
		//strcat(name,head); strcat(name,index[i]);
		fp = fopen(name,"r");
		if(fp == NULL){
			printf("Error: Cannot open file %s\n",name);
			return -1;
		}
		for(int j=0;j<LN;j++){
		    if(fgets(data,sizeof(data), fp) == NULL){
				(void)data;
				return -1;
			}
			y[j] += (atof(data))*(1-i) + (atof(data))*i*1I;
		}
		fclose(fp);
    }
    return 0;
}

int save_data(char *head, double *wn, double complex *arr, int LN, char* subfolder,double t){
    char fname[1024];
    sprintf(fname, "%s/%s_t%f.txt", subfolder, head ,t);
    FILE *fp = fopen(fname, "w");
    for(int i=0; i<LN; i++){
        fprintf(fp, "%lf\t%.16lf\t%.16lf\n", wn[i], creal(arr[i]), cimag(arr[i]) );
    }
    fclose(fp);

    return 0;
}
 */
