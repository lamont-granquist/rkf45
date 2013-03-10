#include <stdlib.h>
#include <stdio.h>

/*Allocate Matrix */
double** allocate_double_matrix(int m, int n) {
  /* Allocate memory for the elements */
  double *mem = malloc(m * n * sizeof(double));
  /* Allocate memory for the matrix array */
  double **mat = malloc(m * sizeof(double *));
  /* Setup array */
  if (mem != NULL && mat != NULL) {
    int i;
    for (i = 0; i < m; ++i) {
      mat[i] = &mem[i * n];
    }
  } else {
    printf("Out of memory!\n"); exit(-1);
  }
  return mat;
}

/* Print Matrix */
void print_matrix(int m,int n,double** mat) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      printf("%.16lf, ",mat[i][j]);
    }
    printf("\n");
  }
}

/* Free Matrix */
void free_double_matrix(double **mat) {
  /* Free elements */
  free(*mat);
  /* Free matrix */
  free(mat);
}
