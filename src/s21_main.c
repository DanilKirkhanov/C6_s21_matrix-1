#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int err = OK;
  if (rows < 1 || columns < 1)
    err = ERR_MATRIX;
  else {
    result->matrix = (double **)calloc(rows, sizeof(double *));
    if (result->matrix == NULL)
      err = 1;
    else {
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = (double *)calloc(columns, sizeof(double));
        if (result->matrix[i] == NULL) {
          err = 1;
          for (int j = 0; j < i; j++) free(result->matrix[j]);
          free(result->matrix);
          break;
        }
      }
    }
    result->columns = columns;
    result->rows = rows;
  }
  return err;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) free(A->matrix[i]);
    free(A->matrix);
    A->matrix = NULL;
  }
  A->rows = 0;
  A->columns = 0;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int eq = 1;
  if (A->columns != B->columns || A->rows != B->rows || s21_valid_matrix(A) ||
      s21_valid_matrix(B))
    eq = 0;
  else {
    for (int i = 0; i < A->rows && eq; i++)
      for (int j = 0; j < A->columns && eq; j++)
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > EPS) eq = 0;
  }
  return eq;
}

int s21_valid_matrix(matrix_t *A) {
  int err = 0;
  if (A == NULL || A->rows < 1 || A->columns < 1 || A->matrix == NULL) err = 1;
  return err;
}

int s21_calculate(matrix_t *A, matrix_t *B, matrix_t *result, int sign) {
  int err = 0;
  if (s21_valid_matrix(A) || s21_valid_matrix(B))
    err = ERR_MATRIX;
  else if (A->rows != B->rows || A->columns != B->columns)
    err = ERR_CALC;
  else if (!s21_create_matrix(A->rows, A->columns, result))
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] + sign * B->matrix[i][j];
  else
    err = ERR_MATRIX;
  return err;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  return s21_calculate(A, B, result, 1);
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  return s21_calculate(A, B, result, -1);
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int err = 0;
  if (s21_valid_matrix(A) || s21_valid_matrix(B))
    err = ERR_MATRIX;
  else if (A->columns != B->rows)
    err = ERR_CALC;
  else if (!s21_create_matrix(A->rows, B->columns, result))
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < B->columns; j++) {
        result->matrix[i][j] = 0;
        for (int k = 0; k < B->rows; k++)
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
      }
  else
    err = ERR_MATRIX;
  return err;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int err = OK;
  if (s21_valid_matrix(A) || s21_create_matrix(A->rows, A->columns, result))
    err = ERR_MATRIX;
  else {
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[i][j] = A->matrix[i][j] * number;
  }
  return err;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int err = OK;
  if (s21_valid_matrix(A) || s21_create_matrix(A->columns, A->rows, result))
    err = ERR_MATRIX;
  else {
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[j][i] = A->matrix[i][j];
  }
  return err;
}
int s21_determinant(matrix_t *A, double *result) {
  int err = OK;
  if (s21_valid_matrix(A))
    err = ERR_MATRIX;
  else if (A->rows != A->columns)
    err = ERR_CALC;
  else if (A->rows == 1)
    *result = A->matrix[0][0];
  else if (A->rows == 2)
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  else if (A->rows > 2) {
    *result = 0;
    matrix_t temp;
    int sign = 1;
    for (int i = 0; i < A->rows; i++) {
      s21_create_matrix(A->rows - 1, A->rows - 1, &temp);
      s21_minor(A, &temp, i, 0);
      double res = 0;
      s21_determinant(&temp, &res);
      *result += sign * A->matrix[i][0] * res;
      sign *= -1;
      s21_remove_matrix(&temp);
    }
  }
  return err;
}

void s21_minor(matrix_t *src, matrix_t *temp, int null_row, int null_column) {
  int i = 0, j = 0;
  for (int ki = 0; ki < src->rows; ki++)
    for (int kj = 0; kj < src->columns; kj++) {
      if (ki != null_row && null_column != kj) {
        temp->matrix[i][j] = src->matrix[ki][kj];
        j++;
        if (j == src->columns - 1) {
          j = 0;
          i++;
        }
      }
    }
}
int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int err = OK;
  if (s21_valid_matrix(A))
    err = ERR_MATRIX;
  else if (A->rows != A->columns || A->rows == 1 || A->columns == 1)
    err = ERR_CALC;
  else {
    s21_create_matrix(A->rows, A->rows, result);
    matrix_t temp;
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++) {
        double res_minor = 0;
        s21_create_matrix(A->rows - 1, A->rows - 1, &temp);
        s21_minor(A, &temp, i, j);
        s21_determinant(&temp, &res_minor);
        result->matrix[i][j] = res_minor * pow(-1, i + j);

        s21_remove_matrix(&temp);
      }
  }
  return err;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  double det = 0;
  int err = s21_determinant(A, &det);
  if (!err) {
    if (det == 0)
      err = ERR_CALC;
    else {
      if (A->rows == 1) {
        s21_create_matrix(A->rows, A->rows, result);
        result->matrix[0][0] = 1.00 / A->matrix[0][0];
      } else {
        matrix_t transpouse, complem;
        int err_compl;
        err_compl = s21_calc_complements(A, &complem);
        if (!err_compl) {
          s21_transpose(&complem, &transpouse);
          s21_create_matrix(A->rows, A->rows, result);
          for (int i = 0; i < A->rows; i++)
            for (int j = 0; j < A->columns; j++)
              result->matrix[i][j] = transpouse.matrix[i][j] / det;
        } else
          err = ERR_CALC;
        if (transpouse.matrix != NULL) s21_remove_matrix(&transpouse);
        if (complem.matrix != NULL) s21_remove_matrix(&complem);
      }
    }
  }

  return err;
}