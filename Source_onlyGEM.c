typedef struct LLmatrix
{
	int col;
	double value;
	struct LLmatrix *next;
	struct LLmatrix *back;
}LLmatrix;



//create matrix from file
LLmatrix** readFromFile(char *filename, int dim, double* diagarray)
{
	LLmatrix** mtx = (LLmatrix**)calloc(dim, sizeof(LLmatrix*));
	LLmatrix *temp = NULL;
	int ret_code, M, N, nonzero, i, x, y;
	double val;
	FILE *fmtx = fopen(filename, "r");

	memcount_GEM_inputs = memcount_GEM_inputs + dim * sizeof(LLmatrix*)+sizeof(LLmatrix*)+7 * sizeof(int)+sizeof(double);
	if (fmtx == NULL)
	{
		printf("file not found\n");
		exit(1);
	}
	MM_typecode matcode;

	if (mm_read_banner(fmtx, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}
	if (mm_is_complex(matcode) || mm_is_pattern(matcode))
	{
		printf("This application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}
	if ((ret_code = mm_read_mtx_crd_size(fmtx, &M, &N, &nonzero)) != 0)
	{
		printf("Unable to read size: exiting\n");
		exit(1);
	}

	if (M != N)
	{
		printf("Matrix not symmetric. Exiting!");
		exit(1);
	}

	for (i = 0; i < nonzero; i++)
	{
		fscanf(fmtx, "%d %d %lf\n", &x, &y, &val);

		if (mm_is_symmetric(matcode))
		{
			insertsortedLL(mtx, x - 1, y - 1, val);
			if (x != y)
				insertsortedLL(mtx, y - 1, x - 1, val);
			else
				diagarray[x - 1] = val;
		}
		else
		{
			insertsortedLLasym(mtx, x - 1, y - 1, val);
			if (x != y)
				insertsortedLLasym(mtx, y - 1, x - 1, val);
			else
				diagarray[x - 1] = val;
		}
	}
	return mtx;
}

//performs gaussian elimination
double* GaussianElim(LLmatrix **matrix, double* b, int dim, fillinLL **fill, int* order, double* diagarray, int* fillincount_gauss, fillinLL **fillintrack, int* proc)//, int* vertexToOrder) //check later: delete fillin
{
	int ithrow, i, j, index = 0;
	int* processed = (int*)calloc(dim, sizeof(int));
	LLmatrix* currrow = NULL;
	double factor, diag_row = 0.0;
	fillinLL *fillinlist = NULL;
	double *soln = NULL, tempabs = 0.0;
	*proc = dim - 1;

	memcount_GEM = memcount_GEM + (dim + 4) * sizeof(int)+4 * sizeof(double)+sizeof(fillinLL*)+sizeof(LLmatrix*);

	for (i = 0; i < dim; i++)
	{
		ithrow = order[i];
		diag_row = 0.0;
		currrow = matrix[ithrow];
		//get diagonal element of this current row
		diag_row = diagarray[ithrow];
		tempabs = diag_row > 0 ? diag_row : (-1.0 * diag_row);

		if (tempabs > 0.000001) //very small numbers are artifacts of floating point numbers, so treating them as zero
		{
			//traverse adjacency list
			currrow = matrix[ithrow];
			while (currrow != NULL)
			{
				if (currrow->col != ithrow && processed[currrow->col] == 0)  //actually we dont need  processed[currrow->col] == 0... check later
				{
					//attack matrix[currrow->col] row and make matrix[currrow->col, ithrow] element 0
					//factor = M[row2, row1] / M[row1, row1]
					factor = currrow->value / diag_row;
					subtractrows(matrix, ithrow, currrow->col, factor, diagarray, fillincount_gauss, fillintrack);
					b[currrow->col] = b[currrow->col] - factor*b[ithrow];
				}
				currrow = currrow->next;
			}
		}
		else //diagonal element is 0, so terminate
		{
			printf("immature termination\n");
			*proc = i;  //number of rows processed			
			break;
		}
		processed[ithrow] = 1;
	}

	//if all rows have been procesed, attempt to find solutions
	*proc = i;
	if (i == dim)
	{
		int *found = (int*)calloc(dim, sizeof(int)); //initially no variables are solved for
		LLmatrix* backsubstraverse = NULL;
		double sum, coeff = 0;
		soln = (double*)malloc(dim*sizeof(double));
		memcount_GEM = memcount_GEM + dim*sizeof(int)+(dim + 2)*sizeof(double)+sizeof(LLmatrix*);
		for (i = dim - 1; i >= 0; i--)
		{
			ithrow = order[i];
			backsubstraverse = matrix[ithrow];
			sum = 0.0;
			while (backsubstraverse != NULL)
			{
				if (backsubstraverse->col != ithrow)
				{
					sum += soln[backsubstraverse->col] * backsubstraverse->value;
				}
				else
				{
					coeff = backsubstraverse->value;
				}
				backsubstraverse = backsubstraverse->next;
			}
			soln[ithrow] = (b[ithrow] - sum) / coeff;
		}
	}
	return soln;
}

//row2 = row2 - factor*row1. This function subtracts two rows when they are represented sparsely (linked list form)
void subtractrows(LLmatrix **matrix, int row1, int row2, double factor, double* diagarray, int* fillincount_gauss, fillinLL **fillintrack)//, int* vertexToOrder)
{
	LLmatrix* rowLL1 = matrix[row1];
	LLmatrix* rowLL2 = matrix[row2];
	LLmatrix* temp = NULL, *lastrow2 = NULL;
	int freeflag;
	double tempabs = 0.0;
	fillinLL *tempfillin;

	while (rowLL1 != NULL && rowLL2 != NULL)
	{
		freeflag = 0;
		if (rowLL1->col < rowLL2->col)
		{
			//insert temp before rowLL2
			LLmatrix* temp = (LLmatrix*)malloc(sizeof(LLmatrix));
			memcount_GEM = memcount_GEM + sizeof(LLmatrix);
			temp->col = rowLL1->col;
			temp->value = -factor * rowLL1->value;
			temp->next = rowLL2;
			temp->back = rowLL2->back;
			if (temp->back != NULL)
				temp->back->next = temp;
			else
				matrix[row2] = temp;
			rowLL2->back = temp;
			tempabs = temp->value > 0 ? temp->value : (-1.0 * temp->value);

			//diagonal elements that turn non zero are not fillins (in the graph theoritic sense)
			if (row2 != rowLL1->col)
			{
				*fillincount_gauss = *fillincount_gauss + 1;
				tempfillin = (fillinLL*)malloc(sizeof(fillinLL));
				memcount_GEM = memcount_GEM + sizeof(fillinLL);
				tempfillin->varnum = rowLL1->col;
				tempfillin->next = fillintrack[row2];
				tempfillin->val = temp->value;
				fillintrack[row2] = tempfillin;
			}
			if (rowLL1->col == row2)
				diagarray[row2] = temp->value;

			rowLL1 = rowLL1->next;

			continue;
		}
		if (rowLL1->col > rowLL2->col)
		{
			lastrow2 = rowLL2;
			rowLL2 = rowLL2->next;
			continue;
		}
		if (rowLL1->col == rowLL2->col)
		{
			//the 2 values cancel out so we need to free a LL member as value is 0			
			tempabs = rowLL1->value*factor - rowLL2->value;
			tempabs = tempabs > 0 ? tempabs : (-1.0 * tempabs);
			//if(tempabs < 0.000001)
			if (rowLL2->col == row1)
			{
				temp = rowLL2;
				if (rowLL2->back != NULL)
					rowLL2->back->next = rowLL2->next;
				else
					matrix[row2] = rowLL2->next;

				if (rowLL2->next != NULL)
					rowLL2->next->back = rowLL2->back;
				freeflag = 1;

				if (rowLL2->col == row2)
					diagarray[row2] = 0;
			}
			else
			{
				rowLL2->value = rowLL2->value - factor * rowLL1->value;
				if (rowLL2->col == row2)
					diagarray[row2] = rowLL2->value;
			}

			lastrow2 = rowLL2;
			rowLL1 = rowLL1->next;
			rowLL2 = rowLL2->next;
			if (freeflag == 1)
			{
				if (lastrow2->back != NULL)
					lastrow2 = lastrow2->back;
				else
					lastrow2 = NULL;
				free(temp);
			}
		}
	}
	//append rest of row1 to row2
	if (rowLL2 == NULL)
	{
		while (rowLL1 != NULL)
		{
			LLmatrix* temp = (LLmatrix*)malloc(sizeof(LLmatrix));
			memcount_GEM = memcount_GEM + sizeof(LLmatrix);

			temp->col = rowLL1->col;
			temp->value = -factor * rowLL1->value;
			tempabs = temp->value > 0 ? temp->value : (-1.0 * temp->value);

			//diagonal elements that turn non zero are not fillins (in the graph theoritic sense)
			if (row2 != rowLL1->col)
			{
				*fillincount_gauss = *fillincount_gauss + 1;
				tempfillin = (fillinLL*)malloc(sizeof(fillinLL));
				memcount_GEM = memcount_GEM + sizeof(fillinLL);
				tempfillin->varnum = rowLL1->col;
				tempfillin->next = fillintrack[row2];
				tempfillin->val = temp->value;
				fillintrack[row2] = tempfillin;
			}

			temp->back = lastrow2;

			if (rowLL1->col == row2)
				diagarray[row2] = temp->value;

			if (lastrow2 != NULL)
				lastrow2->next = temp;
			else
				matrix[row2] = temp;
			temp->next = NULL;
			lastrow2 = temp;
			rowLL1 = rowLL1->next;
		}
	}
}



//auxilliary function to compare fillins from graph and GEM

//compares fillins generated by graph an fillins generated by GEM
int compareFillins(fillinLL **fill, fillinLL **fillin_gauss, int dim, int *vertexToOrder, int processed)
{
	int i, flag = 0;
	fillinLL *temp, *tempgauss, *temppp, *temp1;

	int delfromhead = 0;
	int successfulmatch = 0;
	for (i = 0; i < dim; i++)
	{
		temppp = fill[i];
		temp = fill[i];
		while (temp != NULL)
		{
			flag = 0;
			delfromhead = 0;
			temp1 = temp->next;
			if (vertexToOrder[temp->varnum] <= processed)
			{
				tempgauss = fillin_gauss[i];
				if (tempgauss != NULL)
				{
					flag = deleteFromLL(fillin_gauss, i, temp->varnum, &delfromhead);
					successfulmatch = successfulmatch + flag;
					if (flag == 0)
						printf("error1");
				}
				else
				{
					printf("error2\n");
				}
				flag = 0;


				tempgauss = fillin_gauss[temp->varnum];
				if (tempgauss != NULL)
				{
					flag = deleteFromLL(fillin_gauss, temp->varnum, i, &delfromhead);
					successfulmatch = successfulmatch + flag;
					if (flag == 0)
						printf("error3");
				}
				else
				{
					printf("error4\n");
				}

				deleteFromLL(fill, i, temp->varnum, &delfromhead);
			}
			if (delfromhead == 1)
				temp = fill[i];
			else
				temp = temp1;
		}
	}
	return successfulmatch;
}

