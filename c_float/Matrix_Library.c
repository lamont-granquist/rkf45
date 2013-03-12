#include <stdlib.h>
#include <stdio.h>

/*Allocate Matrix */
float** allocate_float_matrix(int m, int n) {
  /* Allocate memory for the elements */
  float *mem = malloc(m * n * sizeof(float));
  /* Allocate memory for the matrix array */
  float **mat = malloc(m * sizeof(float *));
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
void print_matrix(int m,int n,float** mat) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      printf("%.16lf, ",mat[i][j]);
    }
    printf("\n");
  }
}

/* Free Matrix */
void free_float_matrix(float **mat) {
  /* Free elements */
  free(*mat);
  /* Free matrix */
  free(mat);
}
