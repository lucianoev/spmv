#include "spmv.h"
#include "stdio.h"
#include <stdlib.h>

/*typedef struct {
    int* fila_inicio;       // Índices de inicio de cada fila
    int* indices_columnas;  // Índices de las columnas de los elementos no cero
    double* val;            // Valores de los elementos no cero
    int num_filas;          // Número de filas
    int num_columnas;       // Número de columnas
    int num_sin_ceros;      // Número de elementos no cero
} MatrizCSR;*/

MatrizCSR convert_to_csr(const double* mat, int num_filas, int num_columnas) {
    MatrizCSR csr;
    csr.num_filas = num_filas;
    csr.num_columnas = num_columnas;

    // Contar números no cero
    int nnz = 0; // número de elementos no cero
    for (int i = 0; i < num_filas; i++) {
        for (int j = 0; j < num_columnas; j++) {
            if (mat[i * num_columnas + j] != 0) {
                nnz++;
            }
        }
    }

    // Asignar memoria para la estructura CSR
    csr.fila_inicio = (int*)malloc((num_filas + 1) * sizeof(int));
    csr.indices_columnas = (int*)malloc(nnz * sizeof(int));
    csr.val = (double*)malloc(nnz * sizeof(double));

    // Llenar la estructura CSR
    int current_nnz = 0;
    csr.fila_inicio[0] = 0;
    for (int i = 0; i < num_filas; i++) {
        for (int j = 0; j < num_columnas; j++) {
            if (mat[i * num_columnas + j] != 0) {
                csr.indices_columnas[current_nnz] = j; // Índice de la columna
                csr.val[current_nnz] = mat[i * num_columnas + j]; // Valor
                current_nnz++;
            }
        }
        csr.fila_inicio[i + 1] = current_nnz; // Número de elementos no cero hasta la fila siguiente
    }

    csr.num_sin_ceros = nnz; // Guardar el número total de elementos no cero
    return csr; // Retorna la matriz CSR
}

//Función para el producto entre matriz y vector 
int my_sparse(const MatrizCSR* matriz_csr, double vector[], double resultado[]){

	//Inicializo vector
	for (int fila = 0; fila < matriz_csr->num_filas; fila++) {
	        resultado[fila] = 0.0;
	}

	//Producto de matriz CSR con vector
	for (int fila = 0; fila < matriz_csr->num_filas; fila++) {
		        for (int idx = matriz_csr->fila_inicio[fila]; idx < matriz_csr->fila_inicio[fila + 1]; idx++) {
				            resultado[fila] += matriz_csr->val[idx] * vector[matriz_csr->indices_columnas[idx]];    
			}
	}
	return 0; // Retorna 0 si todo salió bien
}
