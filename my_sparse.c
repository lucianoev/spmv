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


//Función para el producto entre matriz y vector 
void prod_matriz_vector(const MatrizCSR* matriz_csr, const double* vector, double* resultado){

	//Inicializo vector
	for (int fila = 0; fila < matriz_csr->num_filas; fila++) {
	        resultado[fila] = 0.0;
	}

	//Producto de matriz CSR con vector
	for (int fila = 0; fila < matriz_csr.num_filas; fila++) {
		        for (int idx = matriz_csr.fila_inicio[fila]; idx < matriz_csr.fila_inicio[fila + 1]; idx++) {
				            resultado[fila] += matriz_csr.valores[idx] * vector[matriz_csr.indices_columnas[idx]];    
			}
	}
}
