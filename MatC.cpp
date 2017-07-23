#pragma warning(disable : 4996)
#include "MatC.h"

void reorderOutput(double** output, int dimension);
void swap(double** input, double**output, int dimension, int first_row, int second_row);


Matrix * malloc_mat(int rows, int cols)
{
	Matrix * tot = (Matrix*)malloc(sizeof(Matrix));
	tot->data = (element**)malloc(sizeof(element*) * rows);
	for (int i = 0;i < rows;i++)
	{
		tot->data[i] = (element*)malloc(sizeof(element) * cols);
		memset(tot->data[i], 0, sizeof(element) * cols);
	}
	tot->rows = rows;
	tot->cols = cols;
	return tot;
}
Matrix * eyes(int dim)
{
	Matrix * tot = malloc_mat(dim, dim);
	for (int i = 0;i < dim;i++)tot->data[i][i] = 1;
	return tot;
}
void mat_sum(Matrix * M1, Matrix * M2)
{
	for (int i = 0;i < M1->rows;i++)
		for (int j = 0;j < M1->cols;j++)
			M1->data[i][j] += M2->data[i][j];
}
Matrix * copy_mat(Matrix * from)
{
	Matrix * tot = malloc_mat(from->rows, from->cols);
	for (int i = 0;i < from->rows;i++)
		memcpy(tot->data[i], from->data[i], sizeof(element) * from->cols);
	return tot;
}
void mat_sum(Matrix * M1, Matrix * M2, Matrix *& out)
{
	if (!out) out = malloc_mat(M1->rows, M1->cols);
	else if (!out->data) out = malloc_mat(M1->rows, M1->cols);
	for (int i = 0;i < M1->rows;i++)
		for (int j = 0;j < M1->cols;j++)
			out->data[i][j] = M1->data[i][j] + M2->data[i][j];
}
void mat_minus(Matrix * M1, Matrix * M2, Matrix *& out)
{
	if (!out) out = malloc_mat(M1->rows, M1->cols);
	else if (!out->data) out = malloc_mat(M1->rows, M1->cols);
	for (int i = 0;i < M1->rows;i++)
		for (int j = 0;j < M1->cols;j++)
			out->data[i][j] = M1->data[i][j] - M2->data[i][j];
}
void free_mat(Matrix *& item)
{
	if (item == NULL)return;
	for (int i = 0;i < item->rows;i++)free(item->data[i]);
	free(item->data);
	free(item);
	item = NULL;
}
void mat_multiply(Matrix * M1, double num, Matrix *& out)
{
	if (!out) out = malloc_mat(M1->rows, M1->cols);
	else if (!out->data) out = malloc_mat(M1->rows, M1->cols);
	for (int i = 0;i < M1->rows;i++)
		for (int j = 0;j < M1->cols;j++)
			out->data[i][j] = M1->data[i][j] * num;
}
void mat_multiply(Matrix * M1, Matrix * M2, Matrix *& out)
{
	if (out)free_mat(out);
	out = malloc_mat(M1->rows, M2->cols);
	for (int i = 0;i < M1->rows;i++)
	{
		for (int j = 0;j < M2->cols;j++)
		{
			for (int k = 0;k < M1->cols;k++)
			{
				element num = M1->data[i][k] * M2->data[k][j];
				out->data[i][j] += num;
			}
		}
	}
}
void mat_trans(Matrix * M1, Matrix *& tot)
{
	if (tot)free_mat(tot);
	tot = malloc_mat(M1->cols, M1->rows);
	for (int i = 0; i < M1->rows;i++)
	{
		for (int j = 0;j < M1->cols;j++)
		{
			tot->data[j][i] = M1->data[i][j];
		}
	}
}
void mat_inv(Matrix * M1, Matrix *& out)
{
	Matrix * M1_copy = copy_mat(M1);
	int dimension = M1->cols;
	out = eyes(dimension);
	element ** input = M1_copy->data;
	element ** output = out->data;
	int i, j, k;
	for (i = 0;i<dimension;++i)
	{
		for (j = 0;j<dimension;++j)
		{
			if (input[j][i] != 0)
			{
				swap(input, output, dimension, 0, j);
				break;
			}
		}
		for (j = 0;j<dimension;++j)
		{
			if (j == 0)
			{
				for (k = dimension - 1;k >= 0;--k)
					output[j][k] /= input[j][i];
				for (k = dimension - 1;k >= i;--k)
					input[j][k] /= input[j][i];
			}
			else
			{
				for (k = dimension - 1;k >= 0;--k)
					output[j][k] = output[j][k] - input[j][i] / input[0][i] * output[0][k];
				for (k = dimension - 1;k >= i;--k)
					input[j][k] = input[j][k] - input[j][i] / input[0][i] * input[0][k];
			}
		}
		swap(input, output, dimension, 0, (i + 1) % dimension);
	}
	reorderOutput(output, dimension);
	free_mat(M1_copy);
}

void swap(double** input, double**output, int dimension, int first_row, int second_row)
{
	double* temp_row1 = new double[dimension];
	double* temp_row2 = new double[dimension];
	int i;
	for (i = 0;i<dimension;i++)
	{
		temp_row1[i] = input[first_row][i];
		temp_row2[i] = output[first_row][i];
	}
	for (i = 0;i<dimension;i++)
	{
		input[first_row][i] = input[second_row][i];
		output[first_row][i] = output[second_row][i];
		input[second_row][i] = temp_row1[i];
		output[second_row][i] = temp_row2[i];
	}
	delete[] temp_row1;
	delete[] temp_row2;
}
void reorderOutput(double** output, int dimension)
{
	double**temp = new double*[dimension];
	for (int i = 0;i<dimension;++i)
		temp[i] = new double[dimension];

	for (int i = 1;i<dimension;++i)
		memcpy(temp[i - 1], output[i], sizeof(double)*dimension);
	memcpy(temp[dimension - 1], output[0], sizeof(double)*dimension);

	for (int i = 0;i<dimension;++i)
		memcpy(output[i], temp[i], sizeof(double)*dimension);

	for (int i = 0;i<dimension;++i)
		delete[] temp[i];
	delete[] temp;
}
Matrix * mat_read(const char * filename)
{
	int rows, cols;
	FILE * fp = fopen(filename, "r");
	fscanf(fp, "%d%d", &rows, &cols);
	Matrix * tot = malloc_mat(rows, cols);
	for (int i = 0;i < rows;i++)
		for (int j = 0;j < cols;j++)
			fscanf(fp, "%lf", &tot->data[i][j]);
	fclose(fp);
	return tot;
}
void mat_output(Matrix * m, const char * title)
{
	printf("Matrix %s : rows=%d\tcols=%d\n", title, m->rows, m->cols);
	for (int i = 0;i < m->rows;i++)
	{
		for (int j = 0;j < m->cols;j++)
			printf("%.2lf\t", m->data[i][j]);
		printf("\n");
	}
}
void mat_save(Matrix * m, const char * filename)
{
	FILE * fp = fopen(filename, "w");
	fprintf(fp, "%d\t%d\n", m->rows, m->cols);
	for (int i = 0;i < m->rows;i++)
	{
		for (int j = 0;j < m->cols;j++)
			fprintf(fp, "%lf\t", m->data[i][j]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

void LMS(Matrix * Z, Matrix * H, Matrix * D, Matrix *& X, Matrix *& S, Matrix *& V)
{
	Matrix * Ht = NULL, *Dinv = NULL;
	Matrix * temp1 = NULL, *temp2 = NULL, *temp3 = NULL, *temp4 = NULL, *temp5 = NULL;
	mat_trans(H, Ht);
	mat_inv(D, Dinv);
	mat_multiply(Ht, Dinv, temp1);
	mat_multiply(temp1, H, temp2);
	mat_inv(temp2, temp3);
	mat_multiply(temp3, Ht, temp4);
	mat_multiply(temp4, Dinv, temp5);
	mat_multiply(temp5, Z, X);

	Matrix * temp6 = NULL, *temp7 = NULL;
	Matrix * I = eyes(Z->rows);
	mat_multiply(H, temp5, temp6);
	mat_minus(temp6, I, temp7);
	mat_multiply(temp7, Z, V);

	Matrix * theta = NULL, *vt = NULL, *vts = NULL;
	mat_trans(V, vt);
	mat_multiply(vt, Dinv, vts);
	mat_multiply(vts, V, theta);
	double tt = 1.0 / pow(theta->data[0][0] / (Z->rows - X->rows), 2);
	mat_multiply(temp3, tt, S);

	//mat_save(Z, "Z.txt");
	//mat_save(H, "H.txt");
	//mat_save(V, "V.txt");
	//mat_save(S, "S.txt");
	//mat_save(D, "D.txt");


	free_mat(theta);
	free_mat(temp6);free_mat(I);free_mat(temp7);
	free_mat(Ht);free_mat(Dinv);free_mat(temp1);
	free_mat(temp2);free_mat(temp3);free_mat(temp4);
	free_mat(temp5);
}
void mat_resize(Matrix * m, int n_rows, int n_cols, Matrix *& tot)
{
	if (tot)free_mat(tot);
	tot = malloc_mat(n_rows, n_cols);
	int n = m->rows * m->cols;
	for (int i = 0;i < n;i++)
	{
		int a_rows = i / m->cols;
		int a_cols = i % m->cols;

		int b_rows = i / n_cols;
		int b_cols = i % n_cols;

		tot->data[b_rows][b_cols] = m->data[a_rows][a_cols];
	}
}

void mat_merge_horizontal(Matrix * M1, Matrix * M2, Matrix *& out)
{
	if (out)free_mat(out);
	out = malloc_mat(M1->rows, M1->cols + M2->cols);
	for (int i = 0;i < M1->rows;i++)
	{
		for (int j = 0;j < M1->cols;j++)
		{
			out->data[i][j] = M1->data[i][j];
		}
		for (int j = 0;j < M2->cols;j++)
		{
			out->data[i][j + M1->cols] = M2->data[i][j];
		}
	}
}
void sub_mat(Matrix * m, int start_row, int end_row, int start_col, int end_col, Matrix *& out)
{
	if (out)free_mat(out);
	out = malloc_mat(end_row - start_row + 1, end_col - start_col + 1);
	for (int i = start_row;i <= end_row;i++)
	{
		for (int j = start_col;j <= end_col;j++)
		{
			out->data[i - start_row][j - start_col] = m->data[i][j];
		}
	}
}
