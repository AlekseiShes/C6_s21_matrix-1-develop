#include "s21_matrix.h"

#include <check.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define Mode_SUM 0
#define Mode_SUB 1
#define Mode_Equ 2
#define Mode_MuN 3
#define Mode_Tra 4
#define Mode_ONE 5

#define MAX 3
#define parMAX 4
#define numbMAX 1

/*
typedef struct CollectionTests {
  char *name;
#ifdef Q_OS_MAC
  const TTest *func;
#else
  void (*func)(int _i);
#endif
} CollectionTests;
*/

// from -1 to +1
double getrand() { return (rand() % 2 ? 1 : -1) * (float)rand() / RAND_MAX; }

void remove_matrix(int n, ...) {
  va_list argptr;
  va_start(argptr, n);
  for (int i = 0; i < n; i++) s21_remove_matrix(va_arg(argptr, matrix_t *));
  va_end(argptr);
}

int s21_copy_matrix(matrix_t *A, matrix_t *result) {
  int status = (A != NULL) && (A->columns > 0) && (A->rows > 0) &&
                       (A->matrix) && *(A->matrix)
                   ? OK
                   : ERROR_MATRIX;
  if (status == OK) {
    status = s21_create_matrix(A->rows, A->columns, result);
    for (int r = 0; (status == OK) && (r < A->rows); r++)
      for (int c = 0; (status == OK) && (c < A->columns); c++)
        result->matrix[r][c] = A->matrix[r][c];
  }
  return status;
}

void addnoise(double *element, int *eisEq, int par) {
  if (par > 1) *element += (rand() % 2 ? 1 : -1) * PRECISION * 0.1;
  if (par > 2) {
    *element += (rand() % 2 ? 1 : -1) * PRECISION * 10;
    *eisEq = FAILURE;
  }
}

int Init(int mode, matrix_t *A, matrix_t *B, matrix_t *Res, matrix_t *Expec,
         int par, double number, int *estatus, int *eisEq) {
  int status = FAILURE;
  double v1, v2;
  if (mode == Mode_SUM || mode == Mode_SUB || mode == Mode_Equ ||
      mode == Mode_MuN || mode == Mode_Tra || mode == Mode_ONE) {
    *estatus =
        (A && A->matrix && A->rows > 0 && A->columns > 0) ? OK : ERROR_MATRIX;
    if (*estatus == OK && mode == Mode_ONE)
      *estatus = (A->rows == A->columns) ? OK : ERROR_CALC;
    if (*estatus == OK && mode != Mode_MuN && mode != Mode_ONE)
      *estatus =
          (B && B->matrix && B->rows > 0 && B->columns > 0) ? OK : ERROR_MATRIX;
    if (*estatus == OK && mode != Mode_MuN && mode != Mode_ONE) {
      if (mode == Mode_Tra)
        *estatus =
            (A->rows == B->columns && A->columns == B->rows) ? OK : ERROR_CALC;
      else
        *estatus =
            (A->rows == B->rows && A->columns == B->columns) ? OK : ERROR_CALC;
    }
    s21_copy_matrix(A, Expec);
    if (*estatus == OK && par > 0) {
      for (int i = 0; i < A->rows; i++)
        for (int j = 0; j < A->columns; j++) {
          A->matrix[i][j] = v1 = getrand() * 100;
          if (mode != Mode_MuN && mode != Mode_ONE) {
            if (mode == Mode_Tra)
              B->matrix[j][i] = v1;
            else
              B->matrix[i][j] =
                  (mode == Mode_Equ) ? v1 : (v2 = getrand() * 100);
          }
          if (mode == Mode_SUM) Expec->matrix[i][j] = v1 + v2;
          if (mode == Mode_SUB) Expec->matrix[i][j] = v1 - v2;
          if (mode == Mode_Equ || mode == Mode_Tra) Expec->matrix[i][j] = v1;
          if (mode == Mode_MuN) Expec->matrix[i][j] = v1 * number;
          if (mode == Mode_ONE) Expec->matrix[i][j] = (i == j);
          addnoise(&(Expec->matrix[i][j]), eisEq, par);
        }
      if (A->rows == 1 && A->columns == 1 && fabs(A->matrix[0][0]) < 1E-15 &&
          mode == Mode_ONE)
        *estatus = ERROR_CALC;
    } else {
      if (mode == Mode_ONE) *estatus = ERROR_CALC;
    }
    if (mode == Mode_SUM) status = s21_sum_matrix(A, B, Res);
    if (mode == Mode_SUB) status = s21_sub_matrix(A, B, Res);
    if (mode == Mode_Equ) {
      s21_copy_matrix(B, Res);
      status = *estatus;
    }
    if (mode == Mode_MuN) status = s21_mult_number(A, number, Res);
    if (mode == Mode_Tra) {
      s21_transpose(B, Res);
      status = *estatus;
    }
    if (mode == Mode_ONE) {
      s21_remove_matrix(B);
      status = s21_inverse_matrix(A, B);
      status = (status != OK) ? status : s21_mult_matrix(A, B, Res);
    }
  }
  return status;
}

START_TEST(test_s21_create_matrix) {
  matrix_t Matrix;
  matrix_t *Pointer = NULL;
  int status, expected_status;
  for (int r = -MAX; r < MAX; r++)
    for (int c = -MAX; c < MAX; c++) {
      Pointer = rand() % 4 ? &Matrix : NULL;
      status = s21_create_matrix(r, c, Pointer);
      expected_status = (r > 0 && c > 0 && Pointer != NULL) ? OK : ERROR_MATRIX;
      ck_assert_int_eq(status, expected_status);
      if (status == OK) s21_remove_matrix(Pointer);
    }
}
END_TEST

void test4matrix(int mode) {
  int expected_status = OK, status;
  int expected_isEq = SUCCESS, isEq;
  matrix_t matrix_1 = {NULL, 0, 0}, matrix_2 = {NULL, 0, 0},
           matrix_Res = {NULL, 0, 0}, matrix_Expected = {NULL, 0, 0};
  double num;
  for (int r1 = 1; r1 < MAX; r1++)
    for (int c1 = 1; c1 < MAX; c1++)
      for (int r2 = 1; r2 < MAX; r2++)
        for (int c2 = 1; c2 < MAX; c2++)
          for (int par = 0; par < parMAX; par++)
            for (int numb = -numbMAX; numb < numbMAX + 1; numb++) {
              num = numb + getrand() * 0.5;
              s21_create_matrix(r1, c1, &matrix_1);
              s21_create_matrix(r2, c2, &matrix_2);
              expected_status = OK;
              expected_isEq = SUCCESS;
              status = Init(mode, &matrix_1, &matrix_2, &matrix_Res,
                            &matrix_Expected, par, num, &expected_status,
                            &expected_isEq);

              ck_assert_int_eq(expected_status, status);
              if (expected_status == OK && status == OK) {
                isEq = s21_eq_matrix(&matrix_Res, &matrix_Expected);
                ck_assert_int_eq(isEq, expected_isEq);
              }
              remove_matrix(4, &matrix_1, &matrix_2, &matrix_Res,
                            &matrix_Expected);
            }
}

START_TEST(test_s21_eq_matrix) {
  test4matrix(Mode_Equ);

  matrix_t matrix_1 = {NULL, 0, 0}, matrix_2 = {NULL, 0, 0};
  int res;
  s21_create_matrix(4, 4, &matrix_1);
  s21_create_matrix(4, 4, &matrix_2);
  matrix_1.matrix[2][3] = -1E-2;
  res = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(res, 0);
  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&matrix_2);

  s21_create_matrix(4, 4, &matrix_1);
  s21_create_matrix(4, 4, &matrix_2);
  matrix_1.matrix[2][3] = 21;
  matrix_2.matrix[2][3] = 21 + 1E-15;
  res = s21_eq_matrix(&matrix_1, &matrix_2);
  ck_assert_int_eq(res, 1);
  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&matrix_2);
}
END_TEST

START_TEST(test_s21_sum_matrix) { test4matrix(Mode_SUM); }
END_TEST

START_TEST(test_s21_sub_matrix) {
  test4matrix(Mode_SUB);

  matrix_t matrix_1 = {NULL, 0, 0}, matrix_2 = {NULL, 0, 0},
           matrix_R = {NULL, 0, 0}, matrix_E = {NULL, 0, 0};
  int res;
  s21_create_matrix(4, 4, &matrix_1);
  s21_create_matrix(4, 4, &matrix_2);
  s21_create_matrix(4, 4, &matrix_E);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      matrix_1.matrix[i][j] = i + j;
      matrix_2.matrix[i][j] = j;
      matrix_E.matrix[i][j] = i;
    }
  }
  res = s21_sub_matrix(&matrix_1, &matrix_2, &matrix_R);
  ck_assert_int_eq(res, 0);
  if (res == OK) res = s21_eq_matrix(&matrix_R, &matrix_E);
  ck_assert_int_eq(res, SUCCESS);
  s21_remove_matrix(&matrix_1);
  s21_remove_matrix(&matrix_2);
  s21_remove_matrix(&matrix_R);
  s21_remove_matrix(&matrix_E);
}
END_TEST

START_TEST(test_s21_mult_number) { test4matrix(Mode_MuN); }
END_TEST

START_TEST(test_s21_transpose) { test4matrix(Mode_Tra); }
END_TEST

START_TEST(test_s21_inverse_mult_matrix) { test4matrix(Mode_ONE); }
END_TEST

START_TEST(test_s21_determinant) {
  matrix_t matrix = {NULL, 0, 0};
  s21_create_matrix(3, 3, &matrix);
  int count = 1;
  for (int i = 0; i < matrix.rows; i++)
    for (int j = 0; j < matrix.columns; j++) matrix.matrix[i][j] = count++;
  matrix.matrix[2][2] = 10;
  double res;
  int status = s21_determinant(&matrix, &res);
  ck_assert_uint_eq(res, -3);
  ck_assert_int_eq(status, 0);
  s21_remove_matrix(&matrix);
}
END_TEST

START_TEST(test_s21_calc_complements) {
  double m[3][3] = {{1, 2, 3}, {0, 4, 2}, {5, 2, 1}};
  double r[3][3] = {{0, 10, -20}, {4, -14, 8}, {-8, -2, 4}};
  matrix_t matrix = {NULL, 0, 0}, res_org = {NULL, 0, 0},
           res_s21 = {NULL, 0, 0};
  s21_create_matrix(3, 3, &matrix);
  s21_create_matrix(3, 3, &res_org);
  for (int i = 0; i < matrix.rows; i++)
    for (int j = 0; j < matrix.columns; j++) matrix.matrix[i][j] = m[i][j];
  for (int i = 0; i < res_org.rows; i++)
    for (int j = 0; j < res_org.columns; j++) res_org.matrix[i][j] = r[i][j];
  int status = s21_calc_complements(&matrix, &res_s21);
  ck_assert_int_eq(status, OK);
  if (status == OK) {
    status = s21_eq_matrix(&res_org, &res_s21);
    ck_assert_int_eq(status, SUCCESS);
  }
  s21_remove_matrix(&matrix);
  s21_remove_matrix(&res_s21);
  s21_remove_matrix(&res_org);
}
END_TEST

START_TEST(test_s21_inverse_matrix) {
  double m[3][3] = {{2, 5, 7}, {6, 3, 4}, {5, -2, -3}};
  double r[3][3] = {{1, -1, 1}, {-38, 41, -34}, {27, -29, 24}};
  matrix_t matrix = {NULL, 0, 0}, res_org = {NULL, 0, 0},
           res_s21 = {NULL, 0, 0};
  s21_create_matrix(3, 3, &matrix);
  s21_create_matrix(3, 3, &res_org);
  for (int i = 0; i < matrix.rows; i++)
    for (int j = 0; j < matrix.columns; j++) matrix.matrix[i][j] = m[i][j];
  for (int i = 0; i < matrix.rows; i++)
    for (int j = 0; j < matrix.columns; j++) res_org.matrix[i][j] = r[i][j];
  int status = s21_inverse_matrix(&matrix, &res_s21);
  ck_assert_int_eq(status, OK);
  if (status == OK) {
    status = s21_eq_matrix(&res_org, &res_s21);
    ck_assert_int_eq(status, SUCCESS);
  }
  s21_remove_matrix(&matrix);
  s21_remove_matrix(&res_s21);
  s21_remove_matrix(&res_org);
}
END_TEST

START_TEST(test_s21_mult_matrix) {
  matrix_t mat_1 = {NULL, 0, 0}, mat_2 = {NULL, 0, 0}, res_mult = {NULL, 0, 0},
           res_org = {NULL, 0, 0};
  s21_create_matrix(2, 2, &mat_1);
  s21_create_matrix(2, 2, &mat_2);
  s21_create_matrix(2, 2, &res_org);
  for (int i = 0; i < mat_1.rows; i++)
    for (int j = 0; j < mat_1.columns; j++) {
      mat_1.matrix[i][j] = i + j + 1;
      mat_2.matrix[i][j] = i + j + 1;
    }
  res_org.matrix[0][0] = 5;
  res_org.matrix[0][1] = 8;
  res_org.matrix[1][0] = 8;
  res_org.matrix[1][1] = 13;
  int status = s21_mult_matrix(&mat_1, &mat_2, &res_mult);
  ck_assert_int_eq(status, OK);
  if (status == OK) {
    status = s21_eq_matrix(&res_mult, &res_org);
    ck_assert_int_eq(status, SUCCESS);
  }
  s21_remove_matrix(&mat_1);
  s21_remove_matrix(&mat_2);
  s21_remove_matrix(&res_mult);
  s21_remove_matrix(&res_org);
}
END_TEST

Suite *s21_matrix_suite() {
  Suite *s = suite_create("s21_matrix_test");
  /*  CollectionTests Test[] = {
      {"test_s21_create_matrix", test_s21_create_matrix},
      {"test_s21_eq_matrix", test_s21_eq_matrix},
      {"test_s21_sum_matrix", test_s21_sum_matrix},
      {"test_s21_sub_matrix", test_s21_sub_matrix},
      {"test_s21_mult_number", test_s21_mult_number},
      {"test_s21_transpose", test_s21_transpose},
      {"test_s21_inverse_mult_matrix", test_s21_inverse_mult_matrix}
      //,{"test_s21_mult_matrix", test_s21_mult_matrix}
      //,{"test_s21_calc_complements", test_s21_calc_complements}
      //,{"test_s21_determinant", test_s21_determinant}
      //,{"test_s21_inverse_matrix", test_s21_inverse_matrix}
  };
  for (size_t i = 0; i < sizeof(Test) / sizeof(Test[0]); i++) {
    TCase *tc_s21_matrix = tcase_create(Test[i].name);
    suite_add_tcase(s, tc_s21_matrix);
    tcase_add_test(tc_s21_matrix, Test[i].func);
  }
*/
  TCase *tc_s21_create_matrix = tcase_create("test_s21_create_matrix");
  suite_add_tcase(s, tc_s21_create_matrix);
  tcase_add_test(tc_s21_create_matrix, test_s21_create_matrix);
  TCase *tc_test_s21_eq_matrix = tcase_create("test_s21_eq_matrix");
  suite_add_tcase(s, tc_test_s21_eq_matrix);
  tcase_add_test(tc_test_s21_eq_matrix, test_s21_eq_matrix);
  TCase *tc_test_s21_sum_matrix = tcase_create("test_s21_sum_matrix");
  suite_add_tcase(s, tc_test_s21_sum_matrix);
  tcase_add_test(tc_test_s21_sum_matrix, test_s21_sum_matrix);
  TCase *tc_test_s21_sub_matrix = tcase_create("test_s21_sub_matrix");
  suite_add_tcase(s, tc_test_s21_sub_matrix);
  tcase_add_test(tc_test_s21_sub_matrix, test_s21_sub_matrix);
  TCase *tc_test_s21_mult_number = tcase_create("test_s21_mult_number");
  suite_add_tcase(s, tc_test_s21_mult_number);
  tcase_add_test(tc_test_s21_mult_number, test_s21_mult_number);
  TCase *tc_test_s21_transpose = tcase_create("test_s21_transpose");
  suite_add_tcase(s, tc_test_s21_transpose);
  tcase_add_test(tc_test_s21_transpose, test_s21_transpose);
  TCase *tc_test_s21_inverse_mult_matrix =
      tcase_create("test_s21_inverse_mult_matrix");
  suite_add_tcase(s, tc_test_s21_inverse_mult_matrix);
  tcase_add_test(tc_test_s21_inverse_mult_matrix, test_s21_inverse_mult_matrix);

  TCase *tc_test_s21_determinant = tcase_create("test_s21_determinant");
  suite_add_tcase(s, tc_test_s21_determinant);
  tcase_add_test(tc_test_s21_determinant, test_s21_determinant);
  TCase *tc_test_s21_calc_complements =
      tcase_create("test_s21_calc_complements");
  suite_add_tcase(s, tc_test_s21_calc_complements);
  tcase_add_test(tc_test_s21_calc_complements, test_s21_calc_complements);
  TCase *tc_test_s21_inverse_matrix = tcase_create("test_s21_inverse_matrix");
  suite_add_tcase(s, tc_test_s21_inverse_matrix);
  tcase_add_test(tc_test_s21_inverse_matrix, test_s21_inverse_matrix);
  TCase *tc_s21_mult_matrix = tcase_create("test_s21_mult_matrix");
  suite_add_tcase(s, tc_s21_mult_matrix);
  tcase_add_test(tc_s21_mult_matrix, test_s21_mult_matrix);

  return s;
}

int main() {
  int status = OK;
  srand(time(NULL));
  Suite *s = s21_matrix_suite();
  SRunner *runner = srunner_create(s);
  srunner_run_all(runner, CK_NORMAL);
  status = srunner_ntests_failed(runner);
  srunner_free(runner);
  return (status == OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}