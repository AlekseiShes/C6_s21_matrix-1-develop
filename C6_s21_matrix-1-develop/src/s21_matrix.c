#include "s21_matrix.h"

double s21_get_determinant(matrix_t *A, int *status);
int s21_get_minor(int row, int col, matrix_t *A, matrix_t *result);

/* return code error: OK 0; ERROR_MATRIX 1 */
int CheckCorrectMatrix(matrix_t *A) {
  return ((A != NULL) && (A->columns > 0) && (A->rows > 0) && (A->matrix) &&
          *(A->matrix))
             ? OK
             : ERROR_MATRIX;
}

/* return code error: OK 0; ERROR_MATRIX 1; ERROR_CALC 2 */
int CheckSameMatrix(matrix_t *A, matrix_t *B) {
  int status = CheckCorrectMatrix(A);
  status = (status == OK) ? CheckCorrectMatrix(B) : status;
  if (status == OK && !((A->columns == B->columns) && (A->rows == B->rows)))
    status = ERROR_CALC;
  return status;
}

/* return code error: OK 0; ERROR_MATRIX 1; ERROR_CALC 2 */
int checkMultAbleMatrix(matrix_t *A, matrix_t *B) {
  int status = CheckCorrectMatrix(A);
  status = (status == OK) ? CheckCorrectMatrix(B) : status;
  if (status == OK && !(A->columns == B->rows)) status = ERROR_CALC;
  return status;
}

/* return code error: OK 0; ERROR_MATRIX 1; ERROR_CALC 2 */
int CheckSquerMatrix(matrix_t *A) {
  int status = CheckCorrectMatrix(A);
  if (status == OK && !(A->columns == A->rows)) status = ERROR_CALC;
  return status;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = result ? OK : ERROR_MATRIX;
  if ((status == OK) && rows > 0 && columns > 0) {
    result->columns = columns;
    result->rows = rows;
    result->matrix = (double **)malloc(rows * sizeof(double *));
    status = result->matrix != NULL ? OK : ERROR_MATRIX;
    for (int r = 0; (r < rows) && (status == OK); r++)
      status =
          (!(*(result->matrix + r) = (double *)calloc(columns, sizeof(double))))
              ? ERROR_MATRIX
              : OK;
  } else {
    status = ERROR_MATRIX;
  }
  return status;
}

void s21_remove_matrix(matrix_t *A) {
  if (CheckCorrectMatrix(A) == OK) {
    for (int r = 0; r < A->rows; r++) free(*(A->matrix + r));
    free(A->matrix);
    A->matrix = NULL;
    A->columns = A->rows = 0;
  }
}

/* SUCCESS 1 FAILURE 0 */
int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int isEqual = (CheckSameMatrix(A, B) == OK);
  for (int r = 0; isEqual && (r < A->rows); r++)
    for (int c = 0; isEqual && (c < A->columns); c++)
      isEqual = (fabs(A->matrix[r][c] - B->matrix[r][c]) < PRECISION);
  return isEqual;
}

int sumsub_matrix(int sign, matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = CheckSameMatrix(A, B);
  if (status == OK && (result != NULL) &&
      (status = s21_create_matrix(A->rows, A->columns, result)) == OK) {
    sign = sign ? -1 : 1;
    for (int r = 0; r < A->rows; r++)
      for (int c = 0; c < A->columns; c++)
        result->matrix[r][c] = A->matrix[r][c] + sign * B->matrix[r][c];
  }
  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  return sumsub_matrix(0, A, B, result);
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  return sumsub_matrix(1, A, B, result);
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = CheckCorrectMatrix(A);
  if (status == OK &&
      (status = s21_create_matrix(A->rows, A->columns, result)) == OK)
    for (int r = 0; r < A->rows; r++)
      for (int c = 0; c < A->columns; c++)
        result->matrix[r][c] = A->matrix[r][c] * number;
  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = checkMultAbleMatrix(A, B);
  if (status == OK &&
      (status = s21_create_matrix(A->rows, B->columns, result)) == OK)
    for (int r = 0; r < A->rows; r++)
      for (int c = 0; c < B->columns; c++) {
        result->matrix[r][c] = 0;
        for (int n = 0; n < A->columns; n++)
          result->matrix[r][c] += (A->matrix[r][n] * B->matrix[n][c]);
      }
  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = CheckCorrectMatrix(A);
  if (status == OK &&
      (status = s21_create_matrix(A->columns, A->rows, result)) == OK)
    for (int r = 0; r < A->rows; r++)
      for (int c = 0; c < A->columns; c++)
        result->matrix[c][r] = A->matrix[r][c];
  return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = CheckSquerMatrix(A);
  if (status == OK &&
      (status = s21_create_matrix(A->rows, A->columns, result)) == OK)
    for (int r = 0; r < A->rows; r++)
      for (int c = 0; c < A->columns; c++) {
        matrix_t T = {NULL, 0, 0}; /* T - temperary matrix */
        if (s21_create_matrix(A->rows - 1, A->columns - 1, &T) == OK) {
          if ((status = s21_get_minor(r, c, A, &T)) == OK) {
            result->matrix[r][c] = s21_get_determinant(&T, &status);
            if ((r + c) % 2) result->matrix[r][c] *= (-1);
          }
          s21_remove_matrix(&T);
        } else {
          result->matrix[r][c] = 1;
        }
      }
  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  double det;
  int status = s21_determinant(A, &det);
  if (status == OK) {
    if (fabs(det) > 0) {
      matrix_t Complem = {NULL, 0, 0};
      if ((status = s21_calc_complements(A, &Complem)) == OK) {
        matrix_t Trans = {NULL, 0, 0};
        s21_transpose(&Complem, &Trans);
        s21_mult_number(&Trans, 1 / det, result);
        s21_remove_matrix(&Complem);
        s21_remove_matrix(&Trans);
      }
    } else {
      status = ERROR_CALC;
    }
  }
  return status;
}

int s21_determinant(matrix_t *A, double *result) {
  int status = CheckSquerMatrix(A);
  if (status == OK) *result = s21_get_determinant(A, &status);
  return status;
}

double s21_get_determinant(matrix_t *A, int *status) {
  double det = 0;
  if (*status == OK) {
    if (A->rows == 1) {
      det = A->matrix[0][0];
    } else {
      matrix_t T = {NULL, 0, 0};
      if (s21_create_matrix(A->rows - 1, A->columns - 1, &T) == OK)
        for (int c = 0; (c < A->columns) && (*status == OK); c++)
          if ((*status = s21_get_minor(0, c, A, &T)) == OK) {
            if (c % 2)
              det -= A->matrix[0][c] * s21_get_determinant(&T, status);
            else
              det += A->matrix[0][c] * s21_get_determinant(&T, status);
          }
      s21_remove_matrix(&T);
    }
  }
  return det;
}

int s21_get_minor(int row, int col, matrix_t *A, matrix_t *result) {
  int status = CheckCorrectMatrix(A);
  status = (status == OK) ? CheckCorrectMatrix(result) : status;
  status = (status == OK && !((A->columns - 1 == result->columns) &&
                              (A->rows - 1 == result->rows)))
               ? ERROR_CALC
               : status;
  int minor_row = 0, minor_col;
  if (status == OK)
    for (int r = 0; r < A->rows; r++)
      if (r != row) {
        minor_col = 0;
        for (int c = 0; c < A->columns; c++)
          if (c != col) {
            result->matrix[minor_row][minor_col] = A->matrix[r][c];
            minor_col++;
          }
        minor_row++;
      }
  return status;
}