#define _CRT_SECURE_NO_WARNINGS
#define COMPARE_FILLINS 1
#define ALLOWLEXM 1
#define ALLOWLEXP 1

#include<stdio.h>
#include<stdlib.h>
#include "mmio.h"
#include<time.h>
#include<limits.h>

//structure for header cell of doubly double linked list
typedef struct headercell
{
	int flag;
	struct headercell *head;
	struct vertexcell *vset;
	struct headercell *back;
}headercell;

//structure for adjacency list
typedef struct vertex
{
	int num; //final ordering
	int varnum;  //variable number
	struct vertexcell* vc;
	struct vertexLL* neighbours;
}vertex;

//structure for adjacency list's linked list that shows neighbours
typedef struct vertexLL
{
	struct vertex* v;
	struct vertexLL *next;
}vertexLL;

//structure for body cell of doubly double linked list
typedef struct vertexcell
{
	struct headercell *flag;
	//struct vertex *vert;
	int vertexLocation; //Location of vertex in the vertex array (adj list array)
	struct vertexcell *next;
	struct vertexcell *back;
}vertexcell;

//auxiliary structure for implementing fixlist as a linked list
typedef struct newSetLL
{
	struct headercell* hcell;
	struct newSetLL* next;
}newSetLL;  //points to newly created sets, so that we can later set their header->flag to 0

//structure containing pointers to structures required by lexp (doubly double linked list and adjacency list)
typedef struct init
{
	struct headercell* eqgraph;
	struct vertex* adjlist;
}init;

//structure to contain fillin edges
typedef struct fillinLL
{
	double val;
	int varnum;
	struct fillinLL *next;
}fillinLL;

//structure used for keeping track of labels in lexm
typedef struct lexMstruct
{
	int varnum;
	int label;  //even labels are old, odd represent +1/2
	int reached;  //0 shows unreached, 1 shows reached
	struct lexMstruct *next;
}lexMstruct;

//simple linked list for implementing 'reach' in lexm
typedef struct simpleLL
{
	struct vertex* v;
	struct simpleLL *next;
}simpleLL;

typedef struct LLmatrix
{
	int col;
	double value;
	struct LLmatrix *next;
	struct LLmatrix *back;
}LLmatrix;

//function prototypes
init initialize(float* symmatrix, int dim);
void insertIngraph(vertex *adjlist, int i, int j);
void printAdjList(vertex *adjlist, int dim);
void printGraph(headercell *eqgraph);
int* lexP(init *initial, int dim);
void checkcells(init *initial, int dim);
void fillin(init* initial, int *lexpOrder, int dim);
int* lexM(vertex* adjlist, int dim);
fillinLL **fillin2(vertex* initial, int *lexpOrder, int dim, int *fillincount);
void printFillin(fillinLL **fillin, int dim);
init initializelexP(char *filename, int *dim);
vertex* initFromFile(char *filename, int* dim);
void printreached(lexMstruct **labeltrack, int dim);
void printreach(simpleLL** reach, int dim);
void printlabel(lexMstruct **labeltrack, int dim);
void printcountingsort(lexMstruct **countingSort, int dim2);
void writeFillinToFile(FILE* f, fillinLL **fillin, int dim, int fillincount);
void freeadjlist(vertex* adjlist, int dim);

void insertsortedLL(LLmatrix** row, int x, int y, double temp);
LLmatrix** readFromFile(char *filename, int dim, double* diagarray);
void printMtx(LLmatrix **matrix, int dim, int mode);
double* GaussianElim(LLmatrix **matrix, double* b, int dim, fillinLL **fill, int* order, double* diagarray, int* fillincount_gauss, fillinLL **fillin_gauss, int* processed);//, int* vertexToOrder);

void subtractrows(LLmatrix **matrix, int row1, int row2, double factor, double* diagarray, int* fillincount_gauss, fillinLL **fillintrack);
void printMtxRow(LLmatrix* m);
int compareFillins(fillinLL **fill, fillinLL **fillin_gauss, int dim, int* vertexToOrder, int processed);
int* getOrder(vertex* adjlist, int dim);
int deleteFromLL(fillinLL **fill, int row, int varnum, int *delfromhead);

void printToFile(FILE* f, int fillincount_gauss, fillinLL** fillin_graph, fillinLL** fillin_gauss, double *soln, int successfulmatch, int dim, int processed, int zerodiagonal);
void printFillinsToFile(FILE* f, fillinLL** fillin_gauss, fillinLL** fill, int dim);

void insertsortedLLasym(LLmatrix** mtx, int row, int col, double val);
void insertIngraphasym(vertex *adjlist, int i, int j);
void filloutGraph(vertex* adjlist, fillinLL** fillinlexp, int dim);
fillinLL** makeFillinlistCopy(fillinLL** fillin, int dim);

//global variables to count memory usage
long maxmemcount = 0;
long memcount = 0;

long memcount_GEM_inputs = 0;
long memcount_GEM = 0;

//	memcount_GEM += dim*sizeof(int);
//	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;

int main(int argc, char *argv[])
{
	int dim;
	int i = 0, j, fillincount = 0;
	char *filename = argv[1];
	char *outfile = (char*)calloc(20, sizeof(char));
	int *originalOrder = NULL;
	clock_t start, stop, temp, stop1;
	int *order;
	fillinLL **fill;
	int fcount = 0;
	char type;
	int *vertexToOrder = NULL;
	int fillincount_gauss = 0;
	int processed = 0;
	fillinLL **fillin_gauss_lexp, **fillin_gauss_lexm, **fillin_gauss_orig;

	int lexmfillincount = INT_MAX, origfillincount = INT_MAX;
	fillinLL **fillinOrigOrder = NULL;
	fillinLL **fillinlexm = NULL;
	int *lexmOrder = NULL;
	int *vertexToOrder_lexm = NULL;
	int *vertexToOrder_orig = NULL;

	int gaussleftover, graphleftover, successfulmatch;
	start = clock();

	float startfloat, stopfloat, timediff;

	srand(time(NULL));

	//generate output filename
	while (1)
	{
		if (argv[1][i] == '.')
			break;
		outfile[i] = argv[1][i];
		i++;
	}
	outfile[i++] = '.';
	outfile[i++] = 'o'; outfile[i++] = 'u'; outfile[i++] = 't';


	FILE* f = fopen(outfile, "w");

	//call function to create adjacency list and doubly double linked list for lexp to use
	init initial = initializelexP(filename, &dim);

	fprintf(f, "File name: %s. Number of variables = %d\n", argv[1], dim);

	//generate b matrix
	double *matrixB_lexp = (double*)malloc(dim * sizeof(double));
	double *matrixB_lexm = (double*)malloc(dim * sizeof(double));
	double *matrixB_orig = (double*)malloc(dim * sizeof(double));
	memcount_GEM_inputs = memcount_GEM_inputs + dim * sizeof(double);
	FILE *fmtx = fopen(filename, "r");
	MM_typecode matcode;
	if (mm_read_banner(fmtx, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}
	fprintf(f, "B matrix: ");
	if (mm_is_real(matcode) == 0)  //integers
	{
		double r;
		for (i = 0; i < dim; i++)
		{
			r = (rand() >> 8) % 50;
			fprintf(f, "%lf ", r);
			matrixB_lexp[i] = r;
			matrixB_lexm[i] = r;
			matrixB_orig[i] = r;
		}
	}
	else
	{
		double r1, r2;
		for (i = 0; i < dim; i++)
		{
			r1 = (rand() >> 8) % 50;
			r2 = (rand() >> 8) % 50;
			if (r2 != 0)
			{
				matrixB_lexp[i] = r1 / r2;
				matrixB_lexm[i] = r1 / r2;
				matrixB_orig[i] = r1 / r2;
				fprintf(f, "%lf ", r1 / r2);
			}
			else
			{
				matrixB_lexp[i] = r1 / (r2 + 1.0);
				matrixB_lexm[i] = r1 / (r2 + 1.0);
				matrixB_orig[i] = r1 / (r2 + 1.0);
				fprintf(f, "%lf ", r1 / (r2 + 1.0));
			}
		}
	}
	fprintf(f, "\n\n");
	fclose(fmtx);

	printf("Matrix: %s, Nodes: %d", argv[1], dim);

	temp = clock();
	fillin_gauss_lexp = (fillinLL**)calloc(dim, sizeof(fillinLL*));
	double *diagarray_lexp = (double*)calloc(dim, sizeof(double));
	memcount_GEM_inputs = memcount_GEM_inputs + dim * sizeof(double);

	printf("\nMemory used for inputs to GEM: %d\n", memcount_GEM_inputs);
	LLmatrix **matrixA = readFromFile(filename, dim, diagarray_lexp);  //read A matrix from file
	stop = clock();
	startfloat = temp / 1.0;
	stopfloat = stop / 1.0;
	timediff = (stopfloat - startfloat) / CLOCKS_PER_SEC;
	printf("\n\n Time taken for reading and storing matrix for GEM is: %f seconds  (%lu, %lu)\n\n", timediff, temp, stop);


#if ALLOWLEXP
	//run lexp
	int *lexpOrder = lexP(&initial, dim);
	int *vertexToOrder_lexp = getOrder(initial.adjlist, dim);

	//calculate fillin edges
	fillinLL **fillinlexp = fillin2(initial.adjlist, lexpOrder, dim, &fillincount);
	int lexpfillincount = fillincount;
	printf("LexP gave %d fillin edges\n", fillincount);
	fprintf(f, "Lexp order: (listed in order of elimination)\n");
	for (i = 0; i < dim; i++)
		fprintf(f, "%d ", lexpOrder[i] + 1);
	fprintf(f, "\nLexP gave %d fillin edges\n", fillincount);

	//Run GEM for lexp ordering
	order = lexpOrder;
	fill = fillinlexp;
	fcount = lexpfillincount;
	vertexToOrder = vertexToOrder_lexp;

	temp = clock();
	memcount_GEM = dim * sizeof(fillinLL*); //accounts for size of fillin array which contains fill in edges generated by GEM
	double *soln_lexp = GaussianElim(matrixA, matrixB_lexp, dim, fill, order, diagarray_lexp, &fillincount_gauss, fillin_gauss_lexp, &processed);//, vertexToOrder);
	printf("Memory used for Gaussian Elimination in lexp ordering = %d", memcount_GEM);
	stop = clock();
	startfloat = temp / 1.0;
	stopfloat = stop / 1.0;
	timediff = (stopfloat - startfloat) / CLOCKS_PER_SEC;
	printf("\nTime taken for running GEM on lexp ordering: %f seconds  (%lu, %lu)", timediff, temp, stop);
	printFillinsToFile(f, fillin_gauss_lexp, fill, dim);

#if COMPARE_FILLINS
	fillinLL** fillinlexp_copy = makeFillinlistCopy(fillinlexp, dim);


	temp = clock();
	successfulmatch = compareFillins(fill, fillin_gauss_lexp, dim, vertexToOrder, processed);  //compare fillins generated by graph and by lexp
	stop = clock();
	startfloat = temp / 1.0;
	stopfloat = stop / 1.0;
	timediff = (stopfloat - startfloat) / CLOCKS_PER_SEC;
	printf("\nTime taken for comparing fillins of GEM and fillins from running Lexp is: %f seconds  (%lu, %lu)\n\n", timediff, temp, stop);

	printToFile(f, fillincount_gauss, fill, fillin_gauss_lexp, soln_lexp, successfulmatch, dim, processed, order[processed] + 1);
#endif
	if (soln_lexp != NULL)
		free(soln_lexp);
	if (diagarray_lexp != NULL)
		free(diagarray_lexp);
	if (lexpOrder != NULL)
		free(lexpOrder);
	if (fillin_gauss_lexp != NULL)
		free(fillin_gauss_lexp);
	//free structures used in lexp
	freeadjlist(initial.adjlist, dim);
#endif


#if ALLOWLEXM
	//create new adjacency list
	vertex* adjlist1 = initFromFile(filename, &dim);
	//run lexm
	lexmOrder = lexM(adjlist1, dim);
	vertexToOrder_lexm = getOrder(adjlist1, dim);

	fillincount = 0;
	//calculate fillin for lexm
	fillinlexm = fillin2(adjlist1, lexmOrder, dim, &fillincount);

#if COMPARE_FILLINS
	fillinLL** fillinlexm_copy = makeFillinlistCopy(fillinlexm, dim);
#endif

	lexmfillincount = fillincount;

	printf("LexM gave %d fillin edges\n", fillincount);
	fprintf(f, "Lexm order: (listed in order of elimination)\n");
	for (i = 0; i < dim; i++)
		fprintf(f, "%d ", lexmOrder[i] + 1);
	fprintf(f, "\nLexM gave %d fillin edges", fillincount);

	fprintf(f, "\n");

	//Run GEM on lexm ordering
	order = lexmOrder;
	fill = fillinlexm;
	fcount = lexmfillincount;
	vertexToOrder = vertexToOrder_lexm;
	fillincount_gauss = 0;

	double *diagarray_lexm = (double*)calloc(dim, sizeof(double));
	matrixA = readFromFile(filename, dim, diagarray_lexm);
	fillin_gauss_lexm = (fillinLL**)calloc(dim, sizeof(fillinLL*));

	temp = clock();
	memcount_GEM = dim * sizeof(fillinLL*); //accounts for size of fillin array which contains fill in edges generated by GEM
	double *soln_lexm = GaussianElim(matrixA, matrixB_lexm, dim, fill, order, diagarray_lexm, &fillincount_gauss, fillin_gauss_lexm, &processed);//, vertexToOrder);
	printf("Memory used for Gaussian Elimination in lexm ordering = %d\n", memcount_GEM);
	stop = clock();
	startfloat = temp / 1.0;
	stopfloat = stop / 1.0;
	timediff = (stopfloat - startfloat) / CLOCKS_PER_SEC;
	printf("Time taken for running GEM on lexm ordering: %f seconds  (%lu, %lu)\n", timediff, temp, stop);

	printFillinsToFile(f, fillin_gauss_lexm, fill, dim);

#if COMPARE_FILLINS
	temp = clock();
	successfulmatch = compareFillins(fill, fillin_gauss_lexm, dim, vertexToOrder, processed); //compare fillins generated by graph and by lexm
	stop = clock();
	startfloat = temp / 1.0;
	stopfloat = stop / 1.0;
	timediff = (stopfloat - startfloat) / CLOCKS_PER_SEC;
	printf("Time taken for comparing fillins of GEM and fillins from running Lexm is: %f seconds  (%lu, %lu)\n", timediff, temp, stop);

	printToFile(f, fillincount_gauss, fill, fillin_gauss_lexm, soln_lexm, successfulmatch, dim, processed, order[processed] + 1);
#endif
	if (soln_lexm != NULL)
		free(soln_lexm);
	if (diagarray_lexm != NULL)
		free(diagarray_lexm);
	if (lexmOrder != NULL)
		free(lexmOrder);
	if (fillin_gauss_lexm != NULL)
		free(fillin_gauss_lexm);

	freeadjlist(adjlist1, dim);

#endif


	//original order
	vertex* adjlist2 = initFromFile(filename, &dim);
	originalOrder = (int*)calloc(dim, sizeof(int));
	memcount += dim*sizeof(int);
	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
	printf("init finished\n");
	for (i = 0; i < dim; i++)
	{
		originalOrder[i] = i;
		adjlist2[i].num = i;
	}
	fillincount = 0;

	vertexToOrder_orig = (int*)malloc(dim * sizeof(int));
	printf("init arrays\n");
	for (i = 0; i < dim; i++)
		vertexToOrder_orig[i] = i;
	printf("start fillin\n");
	fillinOrigOrder = fillin2(adjlist2, originalOrder, dim, &fillincount);

#if COMPARE_FILLINS
	fillinLL** fillinorig_copy = makeFillinlistCopy(fillinOrigOrder, dim);
#endif

	origfillincount = fillincount;
	printf("\nOriginal order gave %d fillin edges\n", fillincount);
	fprintf(f, "Original order: ");
	for (i = 0; i < dim; i++)
		fprintf(f, "%d ", i + 1);
	fprintf(f, "\nOriginal Order gave %d fillin edges\n", fillincount);

	//Run GEM on original ordering
	order = originalOrder;
	fill = fillinOrigOrder;
	fcount = origfillincount;
	vertexToOrder = vertexToOrder_orig;
	fillincount_gauss = 0;

	double *diagarray_orig = (double*)calloc(dim, sizeof(double));
	matrixA = readFromFile(filename, dim, diagarray_orig);
	fillin_gauss_orig = (fillinLL**)calloc(dim, sizeof(fillinLL*));

	temp = clock();
	memcount_GEM = dim * sizeof(fillinLL*); //accounts for size of fillin array which contains fill in edges generated by GEM
	double *soln_orig = GaussianElim(matrixA, matrixB_orig, dim, fill, order, diagarray_orig, &fillincount_gauss, fillin_gauss_orig, &processed);//, vertexToOrder);
	printf("Memory used for Gaussian Elimination in original ordering = %d\n", memcount_GEM);
	stop = clock();
	startfloat = temp / 1.0;
	stopfloat = stop / 1.0;
	timediff = (stopfloat - startfloat) / CLOCKS_PER_SEC;
	printf("Time taken for running GEM on original ordering: %f seconds  (%lu, %lu)\n", timediff, temp, stop);

	printFillinsToFile(f, fillin_gauss_orig, fill, dim);
#if COMPARE_FILLINS
	temp = clock();
	successfulmatch = compareFillins(fill, fillin_gauss_orig, dim, vertexToOrder, processed); //compare fillins generated by graph and by original order
	stop = clock();
	startfloat = temp / 1.0;
	stopfloat = stop / 1.0;
	timediff = (stopfloat - startfloat) / CLOCKS_PER_SEC;
	printf("Time taken for comparing fillins of GEM and fillins from running original order is: %f seconds  (%lu, %lu)\n", timediff, temp, stop);

	printToFile(f, fillincount_gauss, fill, fillin_gauss_orig, soln_orig, successfulmatch, dim, processed, order[processed] + 1);
#endif


#if COMPARE_FILLINS
	init initial_filled ;
#if ALLOWLEXP
	//this section fills out the graph with fillins generated by lexp and runs lexp on it again
	initial_filled = initializelexP(filename, &dim);
	filloutGraph(initial_filled.adjlist, fillinlexp_copy, dim);
	int *lexpOrderfilled = lexP(&initial_filled, dim);
	fprintf(f, "LexP ordering (when run on graph filled out by lexp fillin edges): ");
	for (i = 0; i < dim; i++)
		fprintf(f, "%d ", lexpOrderfilled[i] + 1);
	fillincount = 0;
	fillin2(initial_filled.adjlist, lexpOrderfilled, dim, &fillincount);
	fprintf(f, "\nNo of fillin edges when lexp is run on a graph filled out by lexp fillin edges = %d\n", fillincount);
	if (lexpOrderfilled != NULL)
		free(lexpOrderfilled);

	//this section fills out the graph with fillins generated by lexp and runs lexm on it again
	{
		vertex* adj_filled = initFromFile(filename, &dim);
		filloutGraph(adj_filled, fillinlexp_copy, dim);
		lexpOrderfilled = lexM(adj_filled, dim);
		fprintf(f, "LexM ordering (when run on graph filled out by lexp fillin edges): ");
		for (i = 0; i < dim; i++)
			fprintf(f, "%d ", lexpOrderfilled[i] + 1);
		fillincount = 0;
		fillin2(adj_filled, lexpOrderfilled, dim, &fillincount);
		fprintf(f, "\nNo of fillin edges when lexm is run on a graph filled out by lexp fillin edges = %d\n", fillincount);
		if (lexpOrderfilled != NULL)
			free(lexpOrderfilled);
	}

#endif

#if ALLOWLEXM
	//this section fills out the graph with fillins generated by lexm and runs lexp on it again
	initial_filled = initializelexP(filename, &dim);
	filloutGraph(initial_filled.adjlist, fillinlexm_copy, dim);
	int *lexmOrderfilled = lexP(&initial_filled, dim);
	fprintf(f, "LexP ordering (when run on graph filled out by lexm fillin edges): ");
	for (i = 0; i < dim; i++)
		fprintf(f, "%d ", lexmOrderfilled[i] + 1);
	fillincount = 0;
	fillin2(initial_filled.adjlist, lexmOrderfilled, dim, &fillincount);
	fprintf(f, "\nNo of fillin edges when lexp is run on a graph filled out by lexm fillin edges = %d\n", fillincount);
	if (lexmOrderfilled != NULL)
		free(lexmOrderfilled);

	//this section fills out the graph with fillins generated by lexm and runs lexm on it again
	{
		vertex* adj_filled = initFromFile(filename, &dim);
		filloutGraph(adj_filled, fillinlexm_copy, dim);
		lexmOrderfilled = lexM(adj_filled, dim);
		fprintf(f, "LexM ordering (when run on graph filled out by lexm fillin edges):");
		for (i = 0; i < dim; i++)
			fprintf(f, "%d ", lexmOrderfilled[i] + 1);
		fillincount = 0;
		fillin2(adj_filled, lexmOrderfilled, dim, &fillincount);
		fprintf(f, "\nNo of fillin edges when lexm is run on a graph filled out by lexm fillin edges = %d\n", fillincount);
		if (lexmOrderfilled != NULL)
			free(lexmOrderfilled);
	}
#endif
	//this section fills out the graph with fillins generated by original order and runs lexp on it again
	initial_filled = initializelexP(filename, &dim);
	filloutGraph(initial_filled.adjlist, fillinorig_copy, dim);
	int *origOrderfilled = lexP(&initial_filled, dim);
	fprintf(f, "LexP ordering (when run on graph filled out by original order fillin edges): ");
	for (i = 0; i < dim; i++)
		fprintf(f, "%d ", origOrderfilled[i] + 1);
	fillincount = 0;
	fillin2(initial_filled.adjlist, origOrderfilled, dim, &fillincount);
	fprintf(f, "\nNo of fillin edges when lexp is run on a graph filled out by original order fillin edges = %d\n", fillincount);
	if (origOrderfilled != NULL)
		free(origOrderfilled);

	//this section fills out the graph with fillins generated by original order and runs lexm on it again
	{
		vertex* adj_filled = initFromFile(filename, &dim);
		filloutGraph(adj_filled, fillinorig_copy, dim);
		origOrderfilled = lexM(adj_filled, dim);
		fprintf(f, "LexM ordering (when run on graph filled out by lexm fillin edges):");
		for (i = 0; i < dim; i++)
			fprintf(f, "%d ", origOrderfilled[i] + 1);
		fillincount = 0;
		fillin2(adj_filled, origOrderfilled, dim, &fillincount);
		fprintf(f, "\nNo of fillin edges when lexm is run on a graph filled out by lexm fillin edges = %d\n", fillincount);
	}
#endif

	stop1 = clock();
	startfloat = start / 1.0;
	stopfloat = stop1 / 1.0;
	timediff = (stopfloat - startfloat) / CLOCKS_PER_SEC;
	printf("\n\n Total time taken is: %f seconds (%lu, %lu)\n", timediff, start, stop1);

	fclose(f);
	printf("bye!");
	return 0;
}

//this function takes an adjacency list and a list of fill ins and creates a 'perfect' graph by adding the fillin edges
void filloutGraph(vertex* adjlist, fillinLL** fillin, int dim)
{
	fillinLL* temp = NULL;
	int i;

	for (i = 0; i < dim; i++)
	{
		temp = fillin[i];
		while (temp != NULL)
		{
			insertIngraph(adjlist, i, temp->varnum);
			temp = temp->next;
		}
	}
}

//makes a copy of fillin list
fillinLL** makeFillinlistCopy(fillinLL** fillin, int dim)
{
	int i;
	fillinLL* temp;
	fillinLL** copy = (fillinLL**)calloc(dim, sizeof(fillinLL*));
	for (i = 0; i < dim; i++)
	{
		temp = fillin[i];
		while (temp != NULL)
		{
			fillinLL* newnode = (fillinLL*)calloc(1, sizeof(fillinLL));
			newnode->varnum = temp->varnum;
			newnode->next = copy[i];
			copy[i] = newnode;
			temp = temp->next;
		}
	}
	return copy;
}

void printFillinsToFile(FILE* f, fillinLL** fillin_gauss, fillinLL** fill, int dim)
{
	int i;
	fillinLL* temp = NULL;
	fprintf(f, "Fillins generated by graph\n");
	for (i = 0; i < dim; i++)
	{
		temp = fill[i];
		while (temp != NULL)
		{
			if (i < temp->varnum)
				fprintf(f, "%d,%d  ", i + 1, temp->varnum + 1);
			else
				fprintf(f, "%d,%d  ", temp->varnum + 1, i + 1);
			temp = temp->next;
		}
	}
	fprintf(f, "\n");
	fprintf(f, "Fillins generated by GEM\n");
	for (i = 0; i < dim; i++)
	{
		temp = fillin_gauss[i];
		while (temp != NULL)
		{
			fprintf(f, "%d %d %lf, ", i + 1, temp->varnum + 1, temp->val);
			temp = temp->next;
		}
	}
	fprintf(f, "\n");
}

void printToFile(FILE* f, int fillincount_gauss, fillinLL** fillin_graph, fillinLL** fillin_gauss, double *soln, int successfulmatch, int dim, int processed, int zerodiagonal)
{
	int i;
	fillinLL* temp;
	int gaussleftover = 0, graphleftover = 0;

	fprintf(f, "Number of non zero elements generated by GEM = %d\n", fillincount_gauss);

	if (soln != NULL)
	{
		fprintf(f, "SOLUTION:  ");
		for (i = 0; i < dim; i++)
		{
			fprintf(f, "x%d = %lf ", i + 1, soln[i]);
		}
		fprintf(f, "\n");
	}
	else
	{
		fprintf(f, "Immature Termination:  ");
		fprintf(f, "Number of rows processed = %d, ", processed);
		fprintf(f, "Zero encountered at row number %d \n", zerodiagonal);
	}

	fprintf(f, "Number of successful matches between fillins in graphs and non zero elements Generated by GEM: %d\n", successfulmatch);

	//left over edges
	fprintf(f, "fillin edges in graph not found in GEM: ");
	for (i = 0; i < dim; i++)
	{
		temp = fillin_graph[i];
		while (temp != NULL)
		{
			fprintf(f, "%d,%d, ", i + 1, temp->varnum + 1);
			graphleftover++;
			temp = temp->next;
		}
	}
	fprintf(f, "\n");
	fprintf(f, "Number of fillin edges from graph not found in fillins from GEM: %d\n", graphleftover);

	//left over edges
	fprintf(f, "fillin edges in GEM not found in graph: ");
	for (i = 0; i < dim; i++)
	{
		temp = fillin_gauss[i];
		while (temp != NULL)
		{
			fprintf(f, "%d %d %lf, ", i + 1, temp->varnum + 1, temp->val);
			gaussleftover++;
			temp = temp->next;
		}
	}
	fprintf(f, "\n");
	fprintf(f, "Number of fillin edges from graph not found in fillins from GEM: %d\n\n", gaussleftover);
}

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

//a helper function used by compareFillins. deletes elements from linked lists
int deleteFromLL(fillinLL **fill, int row, int varnum, int *delfromhead)
{
	//delete varnum from row
	fillinLL* temp = fill[row];
	fillinLL* temp1 = NULL;
	//delete from head
	if (temp == NULL)
	{
		*delfromhead = 0;
		return 0;
	}
	if (temp->varnum == varnum)
	{
		fill[row] = temp->next;
		free(temp);
		*delfromhead = 1;
		return 1;
	}
	//delete from body
	while (temp->next != NULL)
	{
		if (temp->next->varnum == varnum)
		{
			temp1 = temp->next;
			temp->next = temp->next->next;
			free(temp1);
			*delfromhead = 0;
			return 1;
		}
		temp = temp->next;
	}
	*delfromhead = 0;
	return 0;
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

int* getOrder(vertex* adjlist, int dim)
{
	int *vertexToOrder = (int*)malloc(dim * sizeof(int));
	int i;
	for (i = 0; i < dim; i++)
	{
		vertexToOrder[i] = (adjlist + i)->num;
	}
	return vertexToOrder;
}

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

//function for inserting into asymmetric matrices
void insertsortedLLasym(LLmatrix** mtx, int row, int col, double val)
{
	LLmatrix *traverse, *last = NULL;
	LLmatrix* temp;
	double abstemp1 = 0.0, abstemp2 = 0.0;
	int flag = 0, looped = 0;

	temp = (LLmatrix*)malloc(sizeof(LLmatrix));
	memcount_GEM_inputs = memcount_GEM_inputs + sizeof(LLmatrix);
	temp->value = val;
	temp->col = col;

	traverse = mtx[row];

	if (traverse == NULL)
	{
		mtx[row] = temp;
		temp->next = NULL;
		temp->back = NULL;
		flag = 1;
	}
	else
	{
		while (traverse != NULL)
		{
			if (traverse->col == col) //element already present, so check values and modify
			{
				abstemp1 = traverse->value > 0 ? traverse->value : (-1.0 * traverse->value);
				abstemp2 = val > 0 ? val : (-1.0 * val);
				if (abstemp1 > abstemp2)  //current value is greater than the value that was already present, so replace
				{
					traverse->value = val;
				}
				return;
			}
			if (traverse->col > col)
			{
				temp->back = traverse->back;
				if (temp->back != NULL)
					temp->back->next = temp;
				else
					mtx[row] = temp; //insert at head;
				temp->next = traverse;
				traverse->back = temp;
				flag = 1;
				break;
			}
			last = traverse;
			traverse = traverse->next;
		}
	}

	if (flag == 0) //insert at the end
	{
		last->next = temp;
		temp->back = last;
		temp->next = NULL;
	}
}

//insert into matrix structure in a sorted fashion
void insertsortedLL(LLmatrix** mtx, int row, int col, double val)
{
	LLmatrix *traverse, *last = NULL;
	LLmatrix* temp;
	int flag = 0, looped = 0;

	temp = (LLmatrix*)malloc(sizeof(LLmatrix));
	memcount_GEM_inputs = memcount_GEM_inputs + sizeof(LLmatrix);
	temp->value = val;
	temp->col = col;

	traverse = mtx[row];

	if (traverse == NULL)
	{
		mtx[row] = temp;
		temp->next = NULL;
		temp->back = NULL;
		flag = 1;
	}
	else
	{
		while (traverse != NULL)
		{
			if (traverse->col > col)
			{
				temp->back = traverse->back;
				if (temp->back != NULL)
					temp->back->next = temp;
				else
					mtx[row] = temp; //insert at head;
				temp->next = traverse;
				traverse->back = temp;
				flag = 1;
				break;
			}
			last = traverse;
			traverse = traverse->next;
		}
	}

	if (flag == 0) //insert at the end
	{
		last->next = temp;
		temp->back = last;
		temp->next = NULL;
	}
}

//print a row of the matrix
void printMtxRow(LLmatrix* m)
{
	printf("printing row:\n");
	while (m != NULL)
	{
		printf(" %d:%lf ", m->col, m->value);
		m = m->next;
	}
	printf("\n");
}

void printMtx(LLmatrix **matrix, int dim, int mode)
{
	int i, j, lastcol, flag;
	printf("printing matrix\n");
	LLmatrix *currrow;
	for (i = 0; i < dim; i++)
	{
		currrow = matrix[i];
		printf("%d:: ", i);
		lastcol = 0;
		flag = 0;
		while (currrow != NULL)
		{
			if (mode == 1)
			{
				if (flag == 0)
				{
					for (j = 0; j < currrow->col - lastcol; j++)
						printf(" xx%lf", 0);
				}
				else
				{
					for (j = 0; j < currrow->col - lastcol - 1; j++)
						printf(" xx%lf", 0);
				}
			}
			flag = 1;
			printf(" %lf:%d", currrow->value, currrow->col);
			lastcol = currrow->col;
			currrow = currrow->next;
		}
		if (mode == 1)
		{
			for (j = 0; j < dim - lastcol - 1; j++)
				printf(" xx%lf", 0);
		}
		printf("\n");
	}
}


//frees adjacency list structure's memory
void freeadjlist(vertex* adjlist, int dim)
{
	int i;
	for (i = 0; i < dim; i++)
	{
		vertexLL* temp = adjlist[i].neighbours;
		while (temp != NULL)
		{
			adjlist[i].neighbours = adjlist[i].neighbours->next;
			free(temp);
			temp = adjlist[i].neighbours;
			memcount -= sizeof(vertexLL);
			maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
		}
	}
	free(adjlist);
	memcount -= dim*sizeof(vertexLL);
	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
}

//write fillin edges to file
void writeFillinToFile(FILE* f, fillinLL **fillin, int dim, int fillincount)
{
	fillinLL* temp;
	int i = 0;

	for (i = 0; i < dim; i++)
	{
		temp = fillin[i];
		while (temp != NULL)
		{
			if (i < temp->varnum)
				fprintf(f, "%d,%d ", i + 1, temp->varnum + 1);
			else
				fprintf(f, "%d,%d ", temp->varnum + 1, i + 1);
			temp = temp->next;
		}
	}
	fprintf(f, "\n");
}

//read .mtx file using mmio.c and mmio.h APIs and create adjacency list
vertex* initFromFile(char *filename, int* dim)
{
	init initial;
	int ret_code, M, N, nonzero, i, x, y;
	double val;
	memcount = memcount + sizeof(init)+(7 * sizeof(int)) + sizeof(double);
	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
	FILE *fmtx = fopen(filename, "r");
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

	vertex* adjlist = (vertex*)calloc(N, sizeof(vertex)); //adjacency list.
	memcount = memcount + N*sizeof(vertex*);
	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
	for (i = 0; i < N; i++)
	{
		adjlist[i].num = 0;
		adjlist[i].varnum = i;
	}

	for (i = 0; i < nonzero; i++)
	{
		fscanf(fmtx, "%d %d %f\n", &x, &y, &val);
		//assuming val is always non-zero
		if (x != y) //ignore diagonal elements
		{
			if (mm_is_symmetric(matcode))
				insertIngraph(adjlist, x - 1, y - 1); //mtx market format starts from 1 not 0
			else
				insertIngraphasym(adjlist, x - 1, y - 1);
		}
	}
	*dim = N;
	fclose(fmtx);
	return adjlist;
}

//create adjacency list using initFromFile and create doubly double linked list and connect the 2 structures
init initializelexP(char *filename, int *dim)
{
	headercell *eqgraph, *labelphi;
	int i, j;
	vertexcell *tempvcell;
	init initial;

	vertex* adjlist = initFromFile(filename, dim);
	eqgraph = (headercell*)calloc(1, sizeof(headercell));
	labelphi = (headercell*)calloc(1, sizeof(headercell));
	memcount = memcount + sizeof(vertex*)+2 * sizeof(int)+2 * sizeof(headercell*)+sizeof(init)+2 * sizeof(headercell);
	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
	eqgraph->flag = 0;
	eqgraph->back = NULL;
	eqgraph->head = labelphi;
	eqgraph->vset = NULL;
	labelphi->back = eqgraph;
	labelphi->flag = 0;
	labelphi->head = NULL;
	labelphi->vset = NULL;
	tempvcell = labelphi->vset;

	for (i = 0; i < *dim; i++) //add all vertexcells to philabel
	{
		vertexcell* v = (vertexcell*)calloc(1, sizeof(vertexcell));
		memcount = memcount + sizeof(vertexcell);
		maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
		v->vertexLocation = i;
		v->flag = labelphi;
		v->next = tempvcell;
		v->back = NULL;   //Note in paper implementation it points to headercell. But we can get to headercell using flag.
		if (tempvcell != NULL)
			tempvcell->back = v;
		labelphi->vset = v;
		tempvcell = labelphi->vset;
	}

	//make pointer assignments from adlist's vertices to corresponding cells
	tempvcell = eqgraph->head->vset;
	for (i = *dim - 1; i >= 0; i--)
	{
		adjlist[i].vc = tempvcell;
		tempvcell = tempvcell->next;
	}

	initial.adjlist = adjlist;
	initial.eqgraph = eqgraph;
	return initial;
}

int* lexM(vertex* adjlist, int dim)
{
	int *lexmOrder = (int*)calloc(dim, sizeof(int)); //output of this function
	lexMstruct** labeltrack = (lexMstruct**)calloc(dim, sizeof(lexMstruct*));
	lexMstruct** vertexToLabel = (lexMstruct**)calloc(dim, sizeof(lexMstruct*));
	lexMstruct *temp, *tempCS;
	simpleLL* tempSimpleLL;
	int i, k, m, n;
	simpleLL** reach = (simpleLL**)calloc(dim, sizeof(simpleLL*));
	lexMstruct** countingSort = (lexMstruct**)calloc(2 * dim, sizeof(lexMstruct*));

	memcount = memcount + dim * sizeof(int)+2 * dim*sizeof(lexMstruct*)+4 * sizeof(int)+dim*sizeof(simpleLL*)+2 * dim*sizeof(lexMstruct*);
	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
	//initialize 
	for (i = 0; i < dim; i++)
	{
		temp = (lexMstruct*)calloc(1, sizeof(lexMstruct));
		memcount = memcount + sizeof(lexMstruct);
		maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
		temp->varnum = i;
		temp->label = 0; //check
		temp->next = NULL;
		temp->reached = 0;
		labeltrack[i] = temp;
		vertexToLabel[i] = temp;
	}
	k = 0;
	for (i = dim - 1; i >= 0; i--)
	{
		//at iteration i, ith position of labeltrack contains the largest label
		//assign alpha and alpha inverse mappings
		adjlist[labeltrack[i]->varnum].num = i;
		lexmOrder[i] = labeltrack[i]->varnum;
		labeltrack[i]->reached = 1;

		//traverse adj list of the vertex.
		vertexLL* neighbours = adjlist[labeltrack[i]->varnum].neighbours;
		//update labels of immediate neighbours
		while (neighbours != NULL)
		{
			//consider only unreached vertices
			if (vertexToLabel[neighbours->v->varnum]->reached == 0)
			{
				tempSimpleLL = (simpleLL*)calloc(1, sizeof(simpleLL));
				memcount = memcount + sizeof(simpleLL);
				maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
				tempSimpleLL->v = neighbours->v;

				//add w to reach(l(w))
				//l(w) label of w is vertexToLabel[neighbours->v->varnum]->label, which is going to be even for sure
				vertexToLabel[neighbours->v->varnum]->reached = 1;

				tempSimpleLL->next = reach[(vertexToLabel[neighbours->v->varnum]->label) >> 1];
				reach[(vertexToLabel[neighbours->v->varnum]->label) >> 1] = tempSimpleLL;
				vertexToLabel[neighbours->v->varnum]->label += 1;
			}
			neighbours = neighbours->next;
		}
		simpleLL* tempreach;
		//now try to update labels of nodes that are not immediate neighbours
		for (m = 0; m <= k >> 1; m++)
		{
			tempreach = reach[m];
			//traverse through mth linked list of reach
			while (tempreach != NULL)
			{
				neighbours = tempreach->v->neighbours;
				reach[m] = reach[m]->next;

				while (neighbours != NULL)
				{
					lexMstruct *tempLM = vertexToLabel[neighbours->v->varnum];
					if (tempLM->reached == 0) //find an unreached node
					{
						tempLM->reached = 1; //mark it reached

						tempSimpleLL = (simpleLL*)calloc(1, sizeof(simpleLL));
						memcount = memcount + sizeof(simpleLL);
						maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
						tempSimpleLL->v = neighbours->v;
						if (tempLM->label > 2 * m)
						{
							//add to reach
							tempSimpleLL->next = reach[(vertexToLabel[neighbours->v->varnum]->label) >> 1];
							reach[(vertexToLabel[neighbours->v->varnum]->label) >> 1] = tempSimpleLL;
							vertexToLabel[neighbours->v->varnum]->label += 1;
						}
						else
						{
							//add to reach
							tempSimpleLL->next = reach[m];
							reach[m] = tempSimpleLL;
						}
					}

					neighbours = neighbours->next;
				}

				free(tempreach);
				memcount -= sizeof(simpleLL);
				maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
				tempreach = reach[m];
			}
		}

		//sort by label and readjust labels
		//place in counting sort array
		for (m = 0; m <= i - 1; m++)
		{
			labeltrack[m]->next = countingSort[labeltrack[m]->label];
			countingSort[labeltrack[m]->label] = labeltrack[m];
			labeltrack[m] = NULL;
		}

		//retrieve from counting sort array
		m = 0;
		n = 0;
		int lb = 0, flag = 0;;
		while (n <= i - 1)
		{
			tempCS = countingSort[m];

			while (tempCS != NULL)
			{
				tempCS->label = 2 * lb;
				labeltrack[n] = tempCS;
				countingSort[m] = tempCS->next;
				labeltrack[n]->next = NULL;
				n++;
				tempCS = countingSort[m];
				flag = 1;
			}
			if (flag == 1)
				lb++;
			flag = 0;
			m++;
		}
		k = 2 * (lb - 1);
		for (n = 0; n < i - 1; n++)
			labeltrack[n]->reached = 0;
	}

	return lexmOrder;
}


fillinLL** fillin2(vertex* adjlist, int *lexpOrder, int dim, int *fillincount)
{
	int flg = 1;
	int *test = (int*)calloc(dim, sizeof(int)); //calloc sets it to 0 (ie false)
	int i, k, mv;
	vertex *vert, *nextElimination, *tempvert;
	vertexLL *vLL, *temp, *temp1;
	fillinLL **fillin = (fillinLL**)calloc(dim, sizeof(fillinLL**));
	fillinLL *tempfillinLL, *temp1fillinLL;
	memcount = memcount + (dim + 3)*sizeof(int)+3 * sizeof(vertex*);
	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
	int present = 1;
	vertexLL* check = NULL;
	for (i = 0; i < dim - 1; i++)
	{
		k = dim - 1;
		vert = adjlist + lexpOrder[i];

		vLL = vert->neighbours;
		if (vLL == NULL)
			continue;

		//for first node (wich cant be a duplicate as its only 1 node)
		if (vLL->v->num > i)
		{
			test[vLL->v->num] = 1;
			k = vLL->v->num;
		}
		else
			k = dim - 1;

		while (vLL->next != NULL) //eliminate duplicates
		{
			if (test[vLL->next->v->num] == 1) //delete duplicate
			{
				tempfillinLL = fillin[lexpOrder[i]];
				temp = vLL->next;
				vLL->next = temp->next;
				free(temp);
				memcount = memcount - sizeof(vertexLL);
				maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
				continue;
			}
			else
			{
				test[vLL->next->v->num] = 1;
				if (vLL->next->v->num > i)
					k = k < vLL->next->v->num ? k : vLL->next->v->num;
			}
			vLL = vLL->next;
			printf("");
		}
		if (vLL->v->num > i)
			k = k < vLL->v->num ? k : vLL->v->num;

		tempvert = adjlist + lexpOrder[k];

		mv = tempvert->varnum;
		nextElimination = adjlist + lexpOrder[k];
		vLL = vert->neighbours;
		while (vLL != NULL)
		{
			present = 0;
			temp1 = vLL->next;
			test[vLL->v->num] = 0;  //reset test[]
			//add to the neighbour that is going to be eliminated next
			if (vLL->v->varnum != mv && vLL->v->num > i)
			{
				check = (adjlist + lexpOrder[k])->neighbours;
				while (check != NULL)
				{
					if (vLL->v->varnum == check->v->varnum)
					{
						present = 1;
						break;
					}
					check = check->next;
				}
				tempvert = adjlist + lexpOrder[k];
				temp = tempvert->neighbours;
				tempvert->neighbours = vLL;
				vLL->next = temp;

				if (present == 0)
				{
					fillinLL* newnode = (fillinLL*)calloc(1, sizeof(fillinLL));
					memcount = memcount + sizeof(fillinLL);
					maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
					newnode->varnum = vLL->v->varnum;
					newnode->next = fillin[mv];
					fillin[mv] = newnode;
					*fillincount += 1;
				}
			}
			vert->neighbours = temp1;
			vLL = temp1;
		}
	}
	return fillin;
}

int* lexP(init *initial, int dim)
{
	int *lexpOrder = (int*)calloc(dim, sizeof(int));
	int i;
	headercell *eqgraph = initial->eqgraph;
	headercell *eqgraphtemp = eqgraph;
	headercell *tempheadercell;
	vertexcell *vcell, *tempcell;
	vertexLL *neighbours;
	vertexcell *tempvset;
	headercell *newlabel;
	headercell *tempheader;
	newSetLL *fixlist, *tempFL;
	memcount = memcount + 5 * sizeof(headercell*)+sizeof(vertexLL*)+sizeof(vertexcell*)+2 * sizeof(vertexcell*)+2 * sizeof(newSetLL*);
	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
	for (i = dim - 1; i >= 0; i--)
	{
		eqgraphtemp = eqgraph;
		while ((eqgraphtemp->head != NULL) && (eqgraphtemp->head->vset == NULL)) //this while loop deletes empty labels
		{
			tempheadercell = eqgraphtemp->head;
			if (tempheadercell->head != NULL)
			{
				eqgraphtemp->head = tempheadercell->head;
				eqgraphtemp->head->back = eqgraphtemp;
			}
			else
			{
				eqgraphtemp->head = NULL;
			}
			free(tempheadercell);
			memcount = memcount - sizeof(headercell);
			maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
		}

		vcell = eqgraph->head->vset;

		eqgraph->head->vset = vcell->next;
		if (eqgraph->head->vset != NULL)
		{
			eqgraph->head->vset->back = NULL;
		}
		initial->adjlist[vcell->vertexLocation].num = i;
		initial->adjlist[vcell->vertexLocation].vc = NULL;
		lexpOrder[i] = vcell->vertexLocation;

		fixlist = NULL;
		vertexLL *neighbours = initial->adjlist[vcell->vertexLocation].neighbours;
		while (neighbours != NULL)
		{
			tempcell = neighbours->v->vc; //locate cell corresponding to this neighbour

			if (tempcell != NULL)
			{
				//delete cell from current set:
				//delete only if it is unnumbered (but its present in double double LL oly if its unnumbered so looks safe)
				if (tempcell->back != NULL)
				{
					tempcell->back->next = tempcell->next;
				}
				else//tempcell is the first cell in this list
				{
					tempcell->flag->vset = tempcell->next;
				}
				if (tempcell->next != NULL)
				{
					tempcell->next->back = tempcell->back;
				}

				//add cell to new set
				if (tempcell->flag->back->flag == 0)  //create a new label set if needed
				{
					newlabel = (headercell*)calloc(1, sizeof(headercell));
					memcount = memcount + sizeof(headercell);
					maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
					newlabel->flag = 1;
					tempheader = tempcell->flag->back;
					newlabel->head = tempcell->flag;
					newlabel->back = tempheader;
					tempheader->head = newlabel;
					tempcell->flag->back = newlabel;

					//add to fixlist
					newSetLL *tempFixlistElement = (newSetLL*)calloc(1, sizeof(newSetLL));
					memcount = memcount + sizeof(newSetLL);
					maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
					tempFixlistElement->hcell = newlabel;
					tempFixlistElement->next = fixlist;
					fixlist = tempFixlistElement;
				}
				tempvset = tempcell->flag->back->vset;
				tempcell->flag->back->vset = tempcell;
				tempcell->next = tempvset;
				tempcell->flag = tempcell->flag->back;
				tempcell->back = NULL;
				if (tempcell->next != NULL)
					tempcell->next->back = tempcell;
			}

			neighbours = neighbours->next;
		}
		free(vcell);
		memcount -= sizeof(vertexcell);
		maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;

		//update fixlist elements's flags
		tempFL = fixlist;
		while (fixlist != NULL)
		{
			tempFL = fixlist->next;
			fixlist->hcell->flag = 0;
			free(fixlist);
			fixlist = tempFL;
			memcount -= sizeof(newSetLL);
			maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;
		}
	}
	return lexpOrder;
}

void insertIngraphasym(vertex *adjlist, int i, int j)
{
	int flag = 0;
	vertexLL* temp;
	temp = adjlist[i].neighbours;
	while (temp != NULL)
	{
		if (temp->v->varnum == j) //edge has already been added, so dont reinsert
			return;
		temp = temp->next;
	}
	insertIngraph(adjlist, i, j);
}

void insertIngraph(vertex *adjlist, int i, int j)
{
	//insert at adjlist[i] and adjlist[j]
	vertexLL* temp;
	vertexLL* v1 = (vertexLL*)calloc(1, sizeof(vertexLL));
	v1->v = (adjlist + j);
	temp = adjlist[i].neighbours;
	adjlist[i].neighbours = v1;
	v1->next = temp;

	vertexLL* v2 = (vertexLL*)calloc(1, sizeof(vertexLL));
	v2->v = (adjlist + i);
	temp = adjlist[j].neighbours;
	adjlist[j].neighbours = v2;
	v2->next = temp;

	memcount += 2 * sizeof(vertexLL);
	maxmemcount = maxmemcount > memcount ? maxmemcount : memcount;

	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//auxiliary functions (not used in main code)
void printreach(simpleLL** reach, int dim)
{
	int i;
	simpleLL *temp;
	printf("in printreach\n");
	for (i = 0; i < dim; i++)
	{
		printf("%d:  ", i);
		temp = reach[i];
		while (temp != NULL)
		{
			printf("%d ", temp->v->varnum);
			temp = temp->next;
		}
		printf("\n");
	}
	printf("exit printreach\n");
}

void printcountingsort(lexMstruct **countingSort, int dim2)
{
	int i;
	lexMstruct *temp;
	printf("printcountingsort\n");
	for (i = 0; i < dim2; i++)
	{
		temp = countingSort[i];
		printf("\n%d: ", i);
		while (temp != NULL)
		{
			printf("%d ", temp->varnum);
			temp = temp->next;
		}
	}
}

void printlabel(lexMstruct **labeltrack, int dim)
{
	int i;
	printf("printlabels\n");
	for (i = 0; i < dim; i++)
	{
		printf("%d: %d\n", labeltrack[i]->varnum, labeltrack[i]->label);
	}
	printf("\n");
}

void printreached(lexMstruct **labeltrack, int dim)
{
	int i;
	printf("printreached\n");
	for (i = 0; i < dim; i++)
	{
		printf("%d: %d %d\n", i, labeltrack[i]->reached, labeltrack[i]->label);
	}
	printf("\n");
}


void printFillin(fillinLL **fillin, int dim)
{
	fillinLL* temp;
	int i = 0;
	printf("\nenter printFillin\n");
	int count = 0;

	for (i = 0; i < dim; i++)
	{
		temp = fillin[i];
		printf("%d: ", i);
		while (temp != NULL)
		{
			printf(" %d ", temp->varnum);
			temp = temp->next;
			count++;
		}
		printf("\n");
	}
	printf("\nexit printFillin...count = %d\n", count);
}


void printGraph(headercell *eqgraph)
{
	headercell *temp = eqgraph;
	vertexcell *vtemp;
	printf("printgraph\n");
	temp = temp->head;
	printf(" wqwqwqwqwhead \n");
	while (temp != NULL)
	{
		vtemp = temp->vset;
		printf(" wqwqwqwqw ");
		while (vtemp != NULL)
		{
			printf("%d (%d)  ", vtemp->vertexLocation, vtemp->flag->flag);
			vtemp = vtemp->next;
		}
		printf("\n");
		temp = temp->head;
	}
}

void printAdjList(vertex *adjlist, int dim)
{
	int i;
	vertexLL* temp;
	printf("printAdjList\n");
	for (i = 0; i < dim; i++)
	{
		temp = adjlist[i].neighbours;
		printf("%d: ", i);
		while (temp != NULL)
		{
			printf("%d ", temp->v->varnum);
			temp = temp->next;
		}
		printf("\n");
	}
}
