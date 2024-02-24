#ifndef SORT_H
#define SORT_H

void swap_int(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void swap_double(double* a, double* b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

// Partition the array into two sub-arrays and return the pivot index
int partition(double* arr, int* indices, int low, int high) {
    double pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap_double(&arr[i], &arr[j]);
			swap_int(&indices[i], &indices[j]);
        }
    }

    swap_double(&arr[i + 1], &arr[high]);
	swap_int(&indices[i + 1], &indices[high]);
    return (i + 1);
}

// Main quick sort function
void quickSort(double* arr, int* indices, int low, int high) {
    if (low < high) {
        int pi = partition(arr, indices, low, high);
        quickSort(arr, indices,  low,  pi - 1);
        quickSort(arr, indices, pi + 1, high);
    }
}

#endif
