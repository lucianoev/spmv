#ifndef SPMV_H
#define SPMV_H

// Estructura para matriz CSR
typedef struct {
    int* fila_inicio;       // Índices de inicio de cada fila
    int* indices_columnas;  // Índices de las columnas de los elementos no cero
    double* val;            // Valores de los elementos no cero
    int num_filas;          // Número de filas
    int num_columnas;       // Número de columnas
    int num_sin_ceros;      // Número de elementos no cero
} MatrizCSR;

// Declaraciones de funciones
int my_dense(const unsigned int n, const double mat[], double vec[], double result[]);
MatrizCSR convert_to_csr(const double* mat, int num_filas, int num_columnas); // Nueva función de conversión
int my_sparse(const MatrizCSR* matriz_csr, double vec[], double result[]); // Producto de matriz CSR con vector

#endif // SPMV_H

