#include "spmv.h"
#include "stdio.h"


//Estructura para matriz
typedef struct{
	int* fila_inicio;
	int* indices_columnas;
	double* val;
	int num_filas;
	int num_columnas;
	int num_sin_ceros;
} MatrizCSR;

void my_sparse(){

}
