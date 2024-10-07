#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>      // CBLAS in GSL (the GNU Scientific Library)
#include "timer.h"
#include "spmv.h"

#define DEFAULT_SIZE 1024
#define DEFAULT_DENSITY 0.25

unsigned int populate_sparse_matrix(double mat[], unsigned int n, double density, unsigned int seed)
{
    unsigned int nnz = 0;

    srand(seed);

    for (unsigned int i = 0; i < n * n; i++) {
        if ((rand() % 100) / 100.0 < density) {
            // Get a pseudorandom value between -9.99 and 9.99
            mat[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
            nnz++;
        } else {
            mat[i] = 0;
        }
    }

    return nnz;
}

unsigned int populate_vector(double vec[], unsigned int size, unsigned int seed)
{
    srand(seed);

    for (unsigned int i = 0; i < size; i++) {
        vec[i] = ((double)(rand() % 10) + (double)rand() / RAND_MAX) * (rand() % 2 == 0 ? 1 : -1);
    }

    return size;
}

int is_nearly_equal(double x, double y)
{
    const double epsilon = 1e-5; // some small number
    return fabs(x - y) <= epsilon * fabs(x);
}

unsigned int check_result(double ref[], double result[], unsigned int size)
{
    for(unsigned int i = 0; i < size; i++) {
        if (!is_nearly_equal(ref[i], result[i]))
            return 0;
    }
    return 1;
}

/*
// Función para convertir matriz densa a CSR
MatrizCSR convert_to_csr(const double* mat, int num_filas, int num_columnas) {
    // Contar elementos no cero
    int nnz = 0;
    for (int i = 0; i < num_filas; i++) {
        for (int j = 0; j < num_columnas; j++) {
            if (mat[i * num_columnas + j] != 0) {
                nnz++;
            }
        }
    }

    // Asignar memoria para CSR
    MatrizCSR matriz_csr;
    matriz_csr.num_filas = num_filas;
    matriz_csr.num_columnas = num_columnas;
    matriz_csr.num_sin_ceros = nnz;
    matriz_csr.fila_inicio = (int*)malloc((num_filas + 1) * sizeof(int));
    matriz_csr.indices_columnas = (int*)malloc(nnz * sizeof(int));
    matriz_csr.val = (double*)malloc(nnz * sizeof(double));

    // Llenar la estructura CSR
    int idx = 0; // Índice para elementos no cero
    matriz_csr.fila_inicio[0] = 0; // Inicio de la primera fila
    for (int i = 0; i < num_filas; i++) {
        for (int j = 0; j < num_columnas; j++) {
            if (mat[i * num_columnas + j] != 0) {
                matriz_csr.indices_columnas[idx] = j;
                matriz_csr.val[idx] = mat[i * num_columnas + j];
                idx++;
            }
        }
        matriz_csr.fila_inicio[i + 1] = idx; // Guardar índice de inicio de la siguiente fila
    }

    return matriz_csr;
}

// Función para producto de matriz CSR con vector
void my_sparse(const MatrizCSR* matriz_csr, double vec[], double result[]) {
    // Inicializar el vector resultado
    for (int fila = 0; fila < matriz_csr->num_filas; fila++) {
        result[fila] = 0.0;
    }

    // Producto de matriz CSR con vector
    for (int fila = 0; fila < matriz_csr->num_filas; fila++) {
        for (int idx = matriz_csr->fila_inicio[fila]; idx < matriz_csr->fila_inicio[fila + 1]; idx++) {
            result[fila] += matriz_csr->val[idx] * vec[matriz_csr->indices_columnas[idx]];
        }
    }
}*/

// Función principal
int main(int argc, char *argv[])
{
    int size;        // number of rows and cols (size x size matrix)
    double density;  // aprox. ratio of non-zero values

    if (argc < 2) {
        size = DEFAULT_SIZE;
        density = DEFAULT_DENSITY;
    } else if (argc < 3) {
        size = atoi(argv[1]);
        density = DEFAULT_DENSITY;
    } else {
        size = atoi(argv[1]);
        density = atof(argv[2]);
    }

    double *mat, *vec, *refsol, *mysol;
    MatrizCSR matriz_csr;  // Variable para la matriz CSR

    mat = (double *) malloc(size * size * sizeof(double));
    vec = (double *) malloc(size * sizeof(double));
    refsol = (double *) malloc(size * sizeof(double));
    mysol = (double *) malloc(size * sizeof(double));

    unsigned int nnz = populate_sparse_matrix(mat, size, density, 1);
    populate_vector(vec, size, 2);

    printf("Matriz size: %d x %d (%d elements)\n", size, size, size*size);
    printf("%d non-zero elements (%.2lf%%)\n\n", nnz, (double) nnz / (size*size) * 100.0);

    // Cálculo denso usando CBLAS
    printf("Dense computation\n----------------\n");
    timeinfo start, now;
    timestamp(&start);

    cblas_dgemv(CblasRowMajor, CblasNoTrans, size, size, 1.0, mat, size, vec, 1, 0.0, refsol, 1);

    timestamp(&now);
    printf("Time taken by CBLAS dense computation: %ld ms\n", diff_milli(&start, &now));

    // Usando tu propia implementación densa
    timestamp(&start);
    my_dense(size, mat, vec, mysol);
    timestamp(&now);
    printf("Time taken by my dense matrix-vector product: %ld ms\n", diff_milli(&start, &now));

    if (check_result(refsol, mysol, size) == 1)
        printf("Result is ok!\n");
    else
        printf("Result is wrong!\n");

    // Ahora probemos SpMV: Producto matriz dispersa - vector denso
    matriz_csr = convert_to_csr(mat, size, size); // Convertir mat a formato CSR

    // Usar tu propia implementación dispersa
    timestamp(&start);
    my_sparse(&matriz_csr, vec, mysol); // Multiplicación CSR
    timestamp(&now);
    printf("Time taken by my sparse matrix-vector product: %ld ms\n", diff_milli(&start, &now));

    // Comprobar el resultado
    if (check_result(refsol, mysol, size) == 1)
        printf("Sparse result is ok!\n");
    else
        printf("Sparse result is wrong!\n");

    // Liberar recursos
    free(mat);
    free(vec);
    free(refsol);
    free(mysol);
    free(matriz_csr.fila_inicio);
    free(matriz_csr.indices_columnas);
    free(matriz_csr.val);

    return 0;
}
