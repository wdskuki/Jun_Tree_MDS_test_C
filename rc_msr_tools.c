#include "rc_msr_tools.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "main.h"
#include "jerasure.h"
int * generateRandomMatrix(int rows, int cols, int w){
	int w1 = 1;
	int w2 = (int)pow(2, w) - 1;
	int *matrix;  
	int i, j;
	matrix = Malloc(int, rows * cols);
	if(!matrix) {
		printf("Not enough memory!\n");
		exit(1);
	}
	for(i = 0; i < rows * cols; i++)
		matrix[i] = rand()%(w2 +1 - w1) + w1;
	return matrix;
}

/*从a中随机选出n_c个放到c中，但不能是b中的数。其中n_a, n_b, n_c为分别对应数组的个数*/
int random_select(int *a,  int *b, int *c, int n_a , int n_b, int n_c ){
	int *d;
	int n_d;
	int i, j, t, flag;

	if(n_b + n_c> n_a)
		return -1;

	n_d = n_a - n_b;
	d = Malloc(int, n_d);
	if(b != NULL && n_b != 0){
		t = 0;
		for(i = 0; i < n_a; i++){
			flag = 0;
			for(j = 0; j < n_b; j ++){
				if(a[i] == b[j]){
					flag = 1;
					break;
				}
			}
			if(flag == 0){
				d[t] = a[i];
				t++;
			}
		}
	}
	else
	{
		t = 0;
		for(i = 0; i < n_a; i++){
			d[t] = a[i];
			t++;
		}
	}

	if(n_d == n_c){
		for(i = 0; i < n_d; i++)
			c[i] = d[i];
		//free(d);
		return 0;
	}

	i = 0;
	while(i < n_c){
		t = rand()%n_d;
		c[i] = d[t];
		d[t] = d[n_d-1];
		n_d--;
		i++;
	}
	//	free(d);
	return 0;
}

int *store_coding_forward( int *buf, 
										  int rev_n,
										  int r1, 
										  int c1,
										  int *eVector,
										  int r2,
										  int c2,
										  int *p,
										  int r3,
										  int c3,
										  int sed_n_max,
										  int *sed_rows,
										  int w)
{
	//int sed_rows;
	int * random_matrix;
	int * temp_matrix;
	int * local_coding_vector;
	int * send_matrix;
	int i, j, t;
	if(buf == NULL || rev_n == 0 || r1 ==0 || c1 == 0){
		local_coding_vector = jerasure_matrix_multiply(p, eVector, r3, c3, r2, c2,w);
		*sed_rows = 1;
		return local_coding_vector;
		//free(local_coding_vector);
	}

	if(c3 != r2 || r3 != r1 || c2 != c1){
		printf("Error in store_coding_forward()!\n");
		sed_rows = NULL;
		exit(1);
	}
	if(rev_n + 1 >= sed_n_max) *sed_rows = sed_n_max;
	else *sed_rows = rev_n + 1;

	local_coding_vector = jerasure_matrix_multiply(p, eVector, r3, c3, r2, c2,w);
	temp_matrix = Malloc(int, (rev_n + 1)*r1*c1);
	t = rev_n*r1*c1;
	for(i = 0; i < t; i++) temp_matrix[i] = buf[i];
	for(j = t; j < t+r1*c1; j++) temp_matrix[j] = local_coding_vector[j -t];
	free(local_coding_vector);
	random_matrix = generateRandomMatrix((*sed_rows)*r1, (rev_n + 1)*r1, w);
	send_matrix = jerasure_matrix_multiply(random_matrix, temp_matrix, (*sed_rows)*r1, (rev_n + 1)*r1, (rev_n + 1)*r1, c1, w);
	free(temp_matrix);
	free(random_matrix);

	return send_matrix;
}
int * topo_store_coding_forward(int r, 
												  int ** linkGraph,
												  int nodeNum,
												  int r1,  //(D-K+1)	  eVector
												  int c1,		//K*(D-K+1)	eVector
												  int r2,		//1   P:
												  int c2,		//(D-K+1)	P:
												  int sed_n_max,
												  int *eMatrix,
												  int *buf,	//buf
												  int *buf_rows,		// buffer size
												  int *sed_rows,	//send matrix size
												  int w){
	int flag;
	int i, j, k, t;
	int *eVector;
	int *p;
	int *temp;
	int * sub_buf;
	int  r0, c0;
	int rev_n = 0;
	int sub_buf_rows;
	int sub_sed_rows;
	if(c2 != r1){
		printf("Error in topo_store_coding_forward()!\n");
		buf_rows = NULL;
		sed_rows = NULL;
		exit(1);
	}
	eVector = generateRandomMatrix(r1, c1, w);
	for(i = 0; i < r1*c1; i++) 
		eMatrix[r*r1*c1+i] = eVector[i];
	p = generateRandomMatrix(r2, c2, w);

	flag = 0;
	for(i = 0; i < nodeNum; i++)
		if(linkGraph[r][i] != 0){
			flag = 1;break;
		}

	if(flag == 0){
		temp = store_coding_forward(NULL, 0, 0, 0, eVector, r1, c1, p, r2, c2, sed_n_max, &sub_sed_rows, w);
		*buf_rows = 0;
		*sed_rows = sub_sed_rows;
		free(eVector);
		free(p);
		return temp;
	}
	else{
		r0 = r2;
		c0 = c1;
		//buf = Malloc(int, nodeNum*r0*c0);
		k = 0;
		rev_n = 0;
		for(i = 0; i < nodeNum; i++){
			if(linkGraph[r][i] != 0){
				t = i;
				//sub_buf_rows = 0;
				//Malloc(int, nodeNum*K*(D-K+1));
				sub_buf = Malloc(int, nodeNum*c1);
				temp = topo_store_coding_forward(t,linkGraph, nodeNum, r1, c1, r2, c2, sed_n_max,eMatrix, sub_buf, &sub_buf_rows, &sub_sed_rows, w);
				free(sub_buf);
				rev_n += sub_sed_rows;
				for(j = 0; j < sub_sed_rows*r2*c1; j++)
					buf[k++] = temp[j];
			}
		}
//		printf("\nnodeID = %d, buf: \n", r);jerasure_print_matrix(buf, rev_n, r2*c1, W);
		*buf_rows = rev_n;
		temp = store_coding_forward(buf, rev_n, r0, c0, eVector, r1, c1, p, r2, c2, sed_n_max, &sub_sed_rows, w);
		*sed_rows = sub_sed_rows;
		free(eVector);
		free(p);
		//free(buf);
		return temp;
	}
}

int * topo_store_coding_forward2(int r, 
													int ** linkGraph,
													int nodeNum,
													int r1,  //(D-K+1)	  eVector
													int c1,		//K*(D-K+1)	eVector
													int r2,		//1   P:
													int c2,		//(D-K+1)	P:
													int sed_n_max,
													int *eMatrix,
													int *buf,	//buf
													int *buf_rows,		// buffer size
													int *sed_rows,	//send matrix size
													int w){
		int flag;
		int i, j, k, t;
		int *eVector;
		int *p;
		int *temp;
		int * sub_buf;
		int  r0, c0;
		int rev_n = 0;
		int sub_buf_rows;
		int sub_sed_rows;
		if(c2 != r1){
			printf("Error in topo_store_coding_forward()!\n");
			buf_rows = NULL;
			sed_rows = NULL;
			exit(1);
		}
		eVector = Malloc(int, r1*c1);
		//eVector = generateRandomMatrix(r1, c1, w);
		for(i = 0; i < r1*c1; i++) 
				eVector[i] = eMatrix[r*r1*c1+i];
		p = generateRandomMatrix(r2, c2, w);

//		jerasure_print_matrix(eMatrix, N*r1, c1, w);printf("\n");
		flag = 0;
		for(i = 0; i < nodeNum; i++)
			if(linkGraph[r][i] != 0){
				flag = 1;break;
			}

			if(flag == 0){
				temp = store_coding_forward(NULL, 0, 0, 0, eVector, r1, c1, p, r2, c2, sed_n_max, &sub_sed_rows, w);
				*buf_rows = 0;
				*sed_rows = sub_sed_rows;
				free(eVector);
				free(p);
//				jerasure_print_matrix(eMatrix, N*r1, c1, w);printf("\n");
				return temp;
			}
			else{
				r0 = r2;
				c0 = c1;
				//buf = Malloc(int, nodeNum*r0*c0);
				k = 0;
				rev_n = 0;
				for(i = 0; i < nodeNum; i++){
					if(linkGraph[r][i] != 0){
						t = i;
						//sub_buf_rows = 0;
						//Malloc(int, nodeNum*K*(D-K+1));
						sub_buf = Malloc(int, nodeNum*c1);
//						jerasure_print_matrix(eMatrix, N*r1, c1, w);printf("\n");
						temp = topo_store_coding_forward2(t,linkGraph, nodeNum, r1, c1, r2, c2, sed_n_max,eMatrix, sub_buf, &sub_buf_rows, &sub_sed_rows, w);
//						jerasure_print_matrix(eMatrix, N*r1, c1, w);printf("\n");
						free(sub_buf);
						rev_n += sub_sed_rows;
						for(j = 0; j < sub_sed_rows*r2*c1; j++)
							buf[k++] = temp[j];
					}
				}
				//		printf("\nnodeID = %d, buf: \n", r);jerasure_print_matrix(buf, rev_n, r2*c1, W);
				*buf_rows = rev_n;
				temp = store_coding_forward(buf, rev_n, r0, c0, eVector, r1, c1, p, r2, c2, sed_n_max, &sub_sed_rows, w);
				*sed_rows = sub_sed_rows;
				free(eVector);
				free(p);
				//free(buf);
				//free(buf);
			//	jerasure_print_matrix(eMatrix, N*r1, c1, w);
				return temp;
			}
}

int judge_linear_independent(int *eMatrix,
											  int t,
											  int *newComerVector,
											  int testTimes,
											  int n,
											  int d,
											  int k,
											  int w){
	/*test if Linear Independent */
	int *stest, *all, *failed;
	int *stestMatrix;
	int i, j;

	printf("\nTest Repair Successful?: ");
	while(testTimes--){
		stest = Malloc(int, k-1);
		all = Malloc(int, n);
		failed = Malloc(int, 1);
		for(i = 0; i < n; i++) 
			all[i] = i;
		failed[0] = t;
		random_select(all, failed, stest, n, 1, k-1);
		free(all);free(failed);
		stestMatrix = Malloc(int, k*(d-k+1)*k*(d-k+1));
		for(i = 0; i < k; i++){
			if(i < k-1)
				assignMatrixToMatrix(eMatrix, 
				n*(d-k+1), 
				k*(d-k+1), 
				stest[i]*(d-k+1), 
				(stest[i]+1)*(d-k+1)-1, 
				0,
				k*(d-k+1)-1, 
				stestMatrix, 
				i*(d-k+1),
				(i+1)*(d-k+1)-1,
				0, 
				k*(d-k+1)-1);
			else
				assignMatrixToMatrix(newComerVector,
				(d-k+1),
				k*(d-k+1),
				0,
				(d-k+1) -1,
				0,
				k*(d-k+1)-1,
				stestMatrix,
				(k-1)*(d-k+1),
				k*(d-k+1)-1,
				0,
				k*(d-k+1)-1);
		}
//		printf("\nSelected test if Independent Matrix:\n"); jerasure_print_matrix(stestMatrix, k*(d-k+1), k*(d-k+1), w);
		i = jerasure_invertible_matrix(stestMatrix, k*(d-k+1), w);
		if(i != 1){
			//printf("failed repair, the generated matrix is linear dependent with others\n");
			return -1;
		}
		else
			//printf("Repair Successful!\n");
		
		free(stestMatrix);
		free(stest);
	}
	return 0;
}

int judge_retrievable_original_file(int *eMatrix,
									     int t,
										 int *newComerVector,
										 int testTimes,
										 int n,
										 int d,
										 int k,
										 int w)
{
	/*test if Linear Independent */
	int *stest, *all;
	int *stestMatrix;
	int i, j;

	int sucRevNum = 0;
	if(!(t == -1 && newComerVector == NULL))
		assignMatrixToMatrix(newComerVector,(d-k+1),k*(d-k+1),0,(d-k+1)-1,0,k*(d-k+1)-1,eMatrix,t*(d-k+1),(t+1)*(d-k+1)-1,0,k*(d-k+1)-1);

	//printf("\nTest Repair Successful?: ");
	while(testTimes--){
		stest = Malloc(int, k);
		all = Malloc(int, n);
		for(i = 0; i < n; i++) 
			all[i] = i;
		random_select(all, NULL, stest, n, 0, k);
		free(all);
		stestMatrix = Malloc(int, k*(d-k+1)*k*(d-k+1));
		for(i = 0; i < k; i++)
				assignMatrixToMatrix(eMatrix, 
				n*(d-k+1), 
				k*(d-k+1), 
				stest[i]*(d-k+1), 
				(stest[i]+1)*(d-k+1)-1, 
				0,
				k*(d-k+1)-1, 
				stestMatrix, 
				i*(d-k+1),
				(i+1)*(d-k+1)-1,
				0, 
				k*(d-k+1)-1);
		//		printf("\nSelected test if Independent Matrix:\n"); jerasure_print_matrix(stestMatrix, k*(d-k+1), k*(d-k+1), w);
		i = jerasure_invertible_matrix(stestMatrix, k*(d-k+1), w);
		if(i == 1)
			sucRevNum++;
			free(stestMatrix);
		free(stest);
	}

	return sucRevNum;
}

int linear_independent_num(int *eMatrix,
											int n,
											int d,
											int k,
											int avgTimes,
//											int *select,
//											int select_num,
											int w
											)
{
	int * selectMatrix;
	int * selectMatrix2;
	int * all, *select, *select2;
	int i,j,t, i2,j2;
	all = Malloc(int, n);
	for(i = 0; i < n; i++)  all[i] = i;
	for(i = k; i <=n; i++){
		select = Malloc(int, i);
		random_select(all, NULL, select, n, 0, i);
		selectMatrix = Malloc(int, i*(d-k+1)*k*(d-k+1));
		for(j = 0; j < i; j++){
			assignMatrixToMatrix(eMatrix, 
											  n*(d-k+1),
											  k*(d-k+1),
											  select[j]*(d-k+1),
											  (select[j]+1)*(d-k+1)-1,
											  0,
											  k*(d-k+1)-1,
											  selectMatrix,
											  j*(d-k+1),
											  (j+1)*(d-k+1)-1,
											  0,
											  k*(d-k+1)-1);
		}
		
		for(t = 0; t < avgTimes; t++){
			select2 = Malloc(int, k);
			random_select(select, NULL, select2, i, 0, k);
			selectMatrix2 = Malloc(int,  k*(d-k+1)*k*(d-k+1));
			for(i2 = 0; i2 < k; i2++){
				assignMatrixToMatrix(eMatrix, 
 				    i*(d-k+1),
					k*(d-k+1),
					select2[j]*(d-k+1),
					(select2[j]+1)*(d-k+1)-1,
					0,
					k*(d-k+1)-1,
					selectMatrix2,
					i2*(d-k+1),
					(i2+1)*(d-k+1)-1,
					0,
					k*(d-k+1)-1);
					}
			j2 = jerasure_invertible_matrix(selectMatrix2, k*(d-k+1),w);
			if(j2 == 1) return i;
			free(selectMatrix2);
			free(select2);
		}
		free(select);
		free(selectMatrix);
	}
	free(all);
	return (n+1);
}

int linear_independent_num2(int *eMatrix,
	int t,
	int *newComerVector,
	int n,
	int d,
	int k,
	int testTimes,
	//											int *select,
	//											int select_num,
	int w
	)
{
	int * selectMatrix;
	int * selectMatrix2;
	int * all, *select, *select2;
	int i,j, i2,j2;
	
	if(!(t == -1 && newComerVector == NULL))
		assignMatrixToMatrix(newComerVector,(d-k+1),k*(d-k+1),0,(d-k+1)-1,0,k*(d-k+1)-1,eMatrix,t*(d-k+1),(t+1)*(d-k+1)-1,0,k*(d-k+1)-1);
	all = Malloc(int, n);
	for(i = 0; i < n; i++)  all[i] = i;
	for(i = k; i <=n; i++){
		for(i2 = 0; i2 < testTimes*(i-k+1); i2++){
			select = Malloc(int, k);
			random_select(all, NULL, select, n, 0, k);
			selectMatrix = Malloc(int, k*(d-k+1)*k*(d-k+1));
			for(j = 0; j < k; j++){
				assignMatrixToMatrix(eMatrix, 
					n*(d-k+1),
					k*(d-k+1),
					select[j]*(d-k+1),
					(select[j]+1)*(d-k+1)-1,
					0,
					k*(d-k+1)-1,
					selectMatrix,
					j*(d-k+1),
					(j+1)*(d-k+1)-1,
					0,
					k*(d-k+1)-1);
			}
			j2 = jerasure_invertible_matrix(selectMatrix, k*(d-k+1),w);
			if(j2 == 1) return i;
			free(select);
			free(selectMatrix);
		}
	}
	free(all);
	return (n+n/2);
}


void setMatrixValue(int *matrix, 
							   int rows, 
							   int cols,
							   int r1,
							   int r2,
							   int c1,
							   int c2,
							   int value){
	int i, j;
	if(r2 < r1 || c2 < c1|| r2 > rows - 1|| c2 > cols -1 || r1 < 0 || c1 < 0){
		printf("Error in setMatrixValue()!\n");
		exit(1);
	}
	for (i = r1; i <= r2; i++)
		for(j = c1; j <= c2; j++)
			matrix[i*cols + j] = value;
}

void assignMatrixToMatrix(int *matrix,
										 int rows,
										 int cols,
										 int r1,
										 int r2,
										 int c1,
										 int c2,
										 int * matrix2,
										 int r3,
										 int r4,
										 int c3,
										 int c4){
	int i,j,i2,j2;
	if(r2 < r1 || c2 < c1|| r2 > rows - 1|| c2 > cols -1 || r1 < 0 || c1 < 0 || (r2-r1 != r4-r3) || (c2-c1 != c4-c3)){
		printf("Error in assignMatrixToMatrix()!\n");
		exit(1);
	}
	for (i = r3, i2 = r1; i <= r4, i2 <= r2; i++, i2++)
		for(j = c3, j2 = c1; j <= c4, j2 <= c2; j++, j2++)
			matrix2[i*(c4-c3+1) + j] =  matrix[i2*cols+j2];
}

void printMatrix(int *matrix, 
						  int rows,
						  int cols,
						  int r1, 
						  int r2,
						  int c1,
						  int c2){
	int i, j;
	if(r2 < r1 || c2 < c1|| r2 > rows - 1|| c2 > cols -1 || r1 < 0 || c1 < 0){
		printf("Error in printMatrix()!\n");
		exit(1);
	}
	for (i = r1; i <= r2; i++){
		for(j = c1; j <= c2; j++)
			printf("%d ", matrix[i*cols+j]);
		printf("\n");
	}
}

void test_rc_msr_repair_proc_one_time(){
	int i, j;
	int * eMatrix;			//Encoding Matrix
	int t;						//Failed Node ID
	int * all;					//all nodes ID and num
	int * failed;			//failed nodes ID and num
	int * par;				//participant nodes ID and num
	int ** selectedMatrix;	//encoding Matrix in participant nodes
	int ** parMatrix;			//encoidng Matrix generated by participant nodes
	int * dmatrix;				//decoding matrix genereted by new comer
	int * temp;				//temp matrix in Decoding process of New comer
	int * sendMatrix;       //participant nodes send to new comer
	int * newComerMatrix; //new Comer Matrix
	int * stest;						//test if independent node IDs
	int * stestMatrix;			//test if independnet node Matrix
	int testTimes;	//test times if independnet node Matrix
	/*Generate Random Encoding Matrix */
	eMatrix = generateRandomMatrix(N*(D-K+1), K*(D-K+1), W);
	printf("\nGenerate Random Encoding Matrix:\n"); jerasure_print_matrix(eMatrix, N*(D-K+1), K*(D-K+1), W);

	/*One Node Failed*/
	t = rand()%N;
	setMatrixValue(eMatrix, 
							N*(D-K+1),
							K*(D-K+1),
							t*(D-K+1),
							(t+1)*(D-K+1)-1,
							0,
							K*(D-K+1)-1, 0);
	printf("\n Failed Node ID := %d\n", t);
	printf("\nEncoding Matrix After Node Failed (all zeros):\n"); jerasure_print_matrix(eMatrix, N*(D-K+1), K*(D-K+1), W);

	/*Random Select Participant Nodes*/
	all = Malloc(int, N);
	failed = Malloc(int, 1);
	par = Malloc(int, D);
	for(i = 0; i < N; i++) all[i] = i;
	failed[0] = t;
	random_select(all, failed, par, N, 1, D);
	selectedMatrix = Malloc(int *, D);
	for (i = 0; i < D; i++){
		selectedMatrix[i] = Malloc(int, (D-K+1)*K*(D-K+1));
		assignMatrixToMatrix(eMatrix,
										  N*(D-K+1),
										  K*(D-K+1),
										  par[i]*(D-K+1), 
										  (par[i]+1)*(D-K+1)-1,
										  0,
										  K*(D-K+1)-1,
										  selectedMatrix[i],
										  0,
										  (D-K+1)-1,
										  0,
										  K*(D-K+1)-1);
	}
	printf("\nRandom Select Participant Nodes's Encoding Matrix:\n");
	for(i = 0; i < D; i++)
		printMatrix(selectedMatrix[i],
						 (D-K+1),
						 K*(D-K+1),
						 0, 
						 (D-K+1)-1,
						 0,
						 K*(D-K+1)-1);
	free(all);free(failed);free(par);
	/*Random Generate Participant Nodes' Encoding Vector*/
	parMatrix = Malloc(int *, D);
	for(i = 0; i < D; i++)
		parMatrix[i] = generateRandomMatrix(1, (D-K+1), W);
	printf("\nRandom Generate Participant Nodes' Encoding Vector:\n"); 
	for(i = 0; i < D; i++)
		jerasure_print_matrix(parMatrix[i], 1, (D-K+1), W);
	
	/*Random Generate New Comer's Decoding Matrix:*/
	dmatrix = generateRandomMatrix((D-K+1), D, W);
	printf("\nRandom Generate New Comer's Decoding Matrix:\n"); jerasure_print_matrix(dmatrix,(D-K+1), D, W);

	/*Decoding Failed Node's Encoding Matrix*/
	sendMatrix = Malloc(int, D*K*(D-K+1));
	for(i = 0; i < D; i++){
		temp = jerasure_matrix_multiply(parMatrix[i], selectedMatrix[i], 1, (D-K+1), (D-K+1), K*(D-K+1),W);
		assignMatrixToMatrix(temp,
										  1,
										  K*(D-K+1),
										  0,
										  0,
										  0,
										  K*(D-K+1)-1,
										  sendMatrix,
										  i,
										  i,
										  0,
										  K*(D-K+1)-1);	
	}free(temp);
	
	//printf("\nSend Matrix:\n");jerasure_print_matrix(sendMatrix, D, K*(D-K+1), W);
	newComerMatrix = jerasure_matrix_multiply(dmatrix, sendMatrix, (D-K+1), D, D, K*(D-K+1),W);
	printf("\nNew Comer Matrix:\n"); jerasure_print_matrix(newComerMatrix, (D-K+1), K*(D-K+1),W);
	for (i = 0; i < D; i++){
		free(selectedMatrix[i]);
		free(parMatrix[i]);
	}
	free(selectedMatrix);
	free(parMatrix);
	free(dmatrix);
	free(sendMatrix);
	/*Add the New Comer Metrix to original System's Encoding Matrix*/
	assignMatrixToMatrix(newComerMatrix, 
									  (D-K+1),
									  K*(D-K+1),
									  0,
									  (D-K+1)-1,
									  0,
									  K*(D-K+1)-1,
									  eMatrix,
									  t*(D-K+1),
									  (t+1)*(D-K+1)-1,
									  0,
									  K*(D-K+1)-1);
	printf("\nAdd the New Comer Metrix to original System's Encoding Matrix:\n"); jerasure_print_matrix(eMatrix, N*(D-K+1), K*(D-K+1), W);
	
	/*test if Linear Independent */
	testTimes = 100;
	printf("\nTest Repair Successful?:\n");
	while(testTimes--){
	stest = Malloc(int, K-1);
	all = Malloc(int, N);
	failed = Malloc(int, 1);
	for(i = 0; i < N; i++) all[i] = i;
	failed[0] = t;
	random_select(all, failed, stest, N, 1, K-1);
	free(all);free(failed);
	stestMatrix = Malloc(int, K*(D-K+1)*K*(D-K+1));
	for(i = 0; i < K; i++){
		if(i < K-1)
			assignMatrixToMatrix(eMatrix, 
											  N*(D-K+1), 
											  K*(D-K+1), 
											  stest[i]*(D-K+1), 
											  (stest[i]+1)*(D-K+1)-1, 
											  0,
											  K*(D-K+1)-1, 
											  stestMatrix, 
											  i*(D-K+1),
											  (i+1)*(D-K+1)-1,
											  0, 
											  K*(D-K+1)-1);
		else
			assignMatrixToMatrix(newComerMatrix,
											  (D-K+1),
											  K*(D-K+1),
											  0,
											  (D-K+1) -1,
											  0,
											  K*(D-K+1)-1,
											  stestMatrix,
											  (K-1)*(D-K+1),
											  K*(D-K+1)-1,
											  0,
											  K*(D-K+1)-1);
	}
	//printf("\nSelected test if Independent Matrix:\n"); jerasure_print_matrix(stestMatrix, K*(D-K+1), K*(D-K+1), W);
	
	i = jerasure_invertible_matrix(stestMatrix, K*(D-K+1), W);
	if(i != 1){
		printf("failed repair, the generated matrix is linear dependent with others\n");
		exit(1);
	}
	else
		printf("Repair Successful!\n");
	}
	free(newComerMatrix);
	free(stest);
	free(stestMatrix);
	free(eMatrix);
}

void test_rc_msr_repiar_proc_multi_times(){
	int i, j;
	int * eMatrix;			//Encoding Matrix
	int t;						//Failed Node ID
	int * all;					//all nodes ID and num
	int * failed;			//failed nodes ID and num
	int * par;				//participant nodes ID and num
	int ** selectedMatrix;	//encoding Matrix in participant nodes
	int ** parMatrix;			//encoidng Matrix generated by participant nodes
	int * dmatrix;				//decoding matrix genereted by new comer
	int * temp;				//temp matrix in Decoding process of New comer
	int * sendMatrix;       //participant nodes send to new comer
	int * newComerMatrix; //new Comer Matrix
	int * stest;						//test if independent node IDs
	int * stestMatrix;			//test if independnet node Matrix
	int testTimes;	//test times if independnet node Matrix
	int totalSucRepTimes = 0; 

	/*Generate Random Encoding Matrix */
	eMatrix = generateRandomMatrix(N*(D-K+1), K*(D-K+1), W);
	printf("\nGenerate Random Encoding Matrix:\n"); jerasure_print_matrix(eMatrix, N*(D-K+1), K*(D-K+1), W);

	/*One Node Failed*/
	while(1){
		t = rand()%N;
		setMatrixValue(eMatrix, 
			N*(D-K+1),
			K*(D-K+1),
			t*(D-K+1),
			(t+1)*(D-K+1)-1,
			0,
			K*(D-K+1)-1, 0);
//		printf("\n Failed Node ID := %d\n", t);
//		printf("\nEncoding Matrix After Node Failed (all zeros):\n"); jerasure_print_matrix(eMatrix, N*(D-K+1), K*(D-K+1), W);

		/*Random Select Participant Nodes*/
		all = Malloc(int, N);
		failed = Malloc(int, 1);
		par = Malloc(int, D);
		for(i = 0; i < N; i++) all[i] = i;
		failed[0] = t;
		random_select(all, failed, par, N, 1, D);
		selectedMatrix = Malloc(int *, D);
		for (i = 0; i < D; i++){
			selectedMatrix[i] = Malloc(int, (D-K+1)*K*(D-K+1));
			assignMatrixToMatrix(eMatrix,
				N*(D-K+1),
				K*(D-K+1),
				par[i]*(D-K+1), 
				(par[i]+1)*(D-K+1)-1,
				0,
				K*(D-K+1)-1,
				selectedMatrix[i],
				0,
				(D-K+1)-1,
				0,
				K*(D-K+1)-1);
		}
//		printf("\nRandom Select Participant Nodes's Encoding Matrix:\n");
//		for(i = 0; i < D; i++)
//			printMatrix(selectedMatrix[i],
//			(D-K+1),
//			K*(D-K+1),
//			0, 
//			(D-K+1)-1,
//			0,
//			K*(D-K+1)-1);
		free(all);free(failed);free(par);
		/*Random Generate Participant Nodes' Encoding Vector*/
		parMatrix = Malloc(int *, D);
		for(i = 0; i < D; i++)
			parMatrix[i] = generateRandomMatrix(1, (D-K+1), W);
//		printf("\nRandom Generate Participant Nodes' Encoding Vector:\n"); 
//		for(i = 0; i < D; i++)
//			jerasure_print_matrix(parMatrix[i], 1, (D-K+1), W);

		/*Random Generate New Comer's Decoding Matrix:*/
		dmatrix = generateRandomMatrix((D-K+1), D, W);
//		printf("\nRandom Generate New Comer's Decoding Matrix:\n"); jerasure_print_matrix(dmatrix,(D-K+1), D, W);

		/*Decoding Failed Node's Encoding Matrix*/
		sendMatrix = Malloc(int, D*K*(D-K+1));
		for(i = 0; i < D; i++){
			temp = jerasure_matrix_multiply(parMatrix[i], selectedMatrix[i], 1, (D-K+1), (D-K+1), K*(D-K+1),W);
			assignMatrixToMatrix(temp,
				1,
				K*(D-K+1),
				0,
				0,
				0,
				K*(D-K+1)-1,
				sendMatrix,
				i,
				i,
				0,
				K*(D-K+1)-1);	
		}free(temp);

		//printf("\nSend Matrix:\n");jerasure_print_matrix(sendMatrix, D, K*(D-K+1), W);
		newComerMatrix = jerasure_matrix_multiply(dmatrix, sendMatrix, (D-K+1), D, D, K*(D-K+1),W);
//		printf("\nNew Comer Matrix:\n"); jerasure_print_matrix(newComerMatrix, (D-K+1), K*(D-K+1),W);
		for (i = 0; i < D; i++){
			free(selectedMatrix[i]);
			free(parMatrix[i]);
		}
		free(selectedMatrix);
		free(parMatrix);
		free(dmatrix);
		free(sendMatrix);
		/*Add the New Comer Metrix to original System's Encoding Matrix*/
		assignMatrixToMatrix(newComerMatrix, 
			(D-K+1),
			K*(D-K+1),
			0,
			(D-K+1)-1,
			0,
			K*(D-K+1)-1,
			eMatrix,
			t*(D-K+1),
			(t+1)*(D-K+1)-1,
			0,
			K*(D-K+1)-1);
//		printf("\nAdd the New Comer Metrix to original System's Encoding Matrix:\n"); jerasure_print_matrix(eMatrix, N*(D-K+1), K*(D-K+1), W);

		/*test if Linear Independent */
		testTimes = 100;
//		printf("\nTest Repair Successful?:\n");
		while(testTimes--){
			stest = Malloc(int, K-1);
			all = Malloc(int, N);
			failed = Malloc(int, 1);
			for(i = 0; i < N; i++) all[i] = i;
			failed[0] = t;
			random_select(all, failed, stest, N, 1, K-1);
			free(all);free(failed);
			stestMatrix = Malloc(int, K*(D-K+1)*K*(D-K+1));
			for(i = 0; i < K; i++){
				if(i < K-1)
					assignMatrixToMatrix(eMatrix, 
					N*(D-K+1), 
					K*(D-K+1), 
					stest[i]*(D-K+1), 
					(stest[i]+1)*(D-K+1)-1, 
					0,
					K*(D-K+1)-1, 
					stestMatrix, 
					i*(D-K+1),
					(i+1)*(D-K+1)-1,
					0, 
					K*(D-K+1)-1);
				else
					assignMatrixToMatrix(newComerMatrix,
					(D-K+1),
					K*(D-K+1),
					0,
					(D-K+1) -1,
					0,
					K*(D-K+1)-1,
					stestMatrix,
					(K-1)*(D-K+1),
					K*(D-K+1)-1,
					0,
					K*(D-K+1)-1);
			}
			//printf("\nSelected test if Independent Matrix:\n"); jerasure_print_matrix(stestMatrix, K*(D-K+1), K*(D-K+1), W);

			i = jerasure_invertible_matrix(stestMatrix, K*(D-K+1), W);
			if(i != 1){
				printf("failed repair, the generated matrix is linear dependent with others\n");
				exit(1);
			}
//			else
//				printf("Repair Successful!\n");
		}
		free(newComerMatrix);
		free(stest);
		free(stestMatrix);
		totalSucRepTimes++;
		printf(" No. %d is repair successful!\n", totalSucRepTimes);
	}
	free(eMatrix);
}

void test_store_coding_forward(){
	int *buf;
	int rev_n = 5;
	int r1 = 1; 
	int c1 = K*(D-K+1);
	int *eVector;
	int r2 = (D-K+1);
	int c2 = K*(D-K+1);
	int *p;
	int r3 = 1;
	int c3 = (D-K+1);
	int sed_n_max = 6;
	int * send_matrix;
	int sed_rows;

	eVector = generateRandomMatrix(r2, c2, W);
	printf("\neVector:\n"); jerasure_print_matrix(eVector, r2, c2, W);
	
	printf("\nrev_n = %d\nsed_n_max = %d\n", rev_n, sed_n_max);

	buf = generateRandomMatrix(rev_n*r1, c1, W);
	printf("\nbuf:\n"); jerasure_print_matrix(buf, rev_n*r1, c1, W);

	p = generateRandomMatrix(r3, c3, W);
	printf("\np:\n"); jerasure_print_matrix(p, r3, c3, W);

	//sed_rows = (rev_n + 1 > sed_n_max) ? sed_n_max : rev_n + 1;
	send_matrix = store_coding_forward(buf, rev_n, r1, c1, eVector, r2, c2, p, r3, c3, sed_n_max, &sed_rows, W);
	printf("\nsend_matrix:\n"); jerasure_print_matrix(send_matrix, sed_rows * r1, c1, W);

	//printf("send_matrix:%d\n",sizeof(send_matrix));
	free(send_matrix);
	free(eVector);
	free(p);
	free(buf);
}

void test_store_coding_forward2(){
	int *buf = NULL;
	int rev_n = 5;
	int r1 = 1; 
	int c1 = K*(D-K+1);
	int *eVector;
	int r2 = (D-K+1);
	int c2 = K*(D-K+1);
	int *p;
	int r3 = 1;
	int c3 = (D-K+1);
	int sed_n_max = 6;
	int * send_matrix;
	int sed_rows;

	eVector = generateRandomMatrix(r2, c2, W);
	printf("\neVector:\n"); jerasure_print_matrix(eVector, r2, c2, W);

	p = generateRandomMatrix(r3, c3, W);
	printf("\np:\n"); jerasure_print_matrix(p, r3, c3, W);

	send_matrix = store_coding_forward(buf, rev_n, r1, c1, eVector, r2, c2, p, r3, c3, sed_n_max, &sed_rows, W);
	printf("\nsend_matrix:\n"); jerasure_print_matrix(send_matrix, 1 * r1, c1, W);

	free(send_matrix);
	free(eVector);
	free(p);
}

void test_topo_store_coding_forward(){
	int **linkGraph;
	int * rootBuffer;
	int * rootRandomMatrix;
	int * newComerVector;
	int *eMatrix;
	int buf_rows;
	int sed_rows;
	int nodeNum = 5;
	int r = 0;
	int sed_n_max = D-K+1;
	//int sed_n_max = 1;
	int i, j, t;
	linkGraph = Malloc(int *, nodeNum);
	for(i = 0; i < nodeNum; i++){
		linkGraph[i] = Malloc(int, nodeNum);
		memset(linkGraph[i], 0, sizeof(int)*nodeNum);
	}
	rootBuffer = Malloc(int, nodeNum*K*(D-K+1));
	/*
	linkGraph[0][1] = 1;
	linkGraph[1][3] = 1;
	linkGraph[1][4] = 1;
	linkGraph[4][2] = 1;
	*/
	linkGraph[0][1] = 1;
	linkGraph[0][4] = 1;
	linkGraph[1][3] = 1;
	linkGraph[4][2] = 1;
	printf("\nlinkGraph:\n");
	for(i = 0; i < nodeNum; i++){
		for(j = 0; j < nodeNum; j++)
			printf("%d ", linkGraph[i][j]);
		printf("\n");
	}
	eMatrix = Malloc(int, N*(D-K+1)*K*(D-K+1));
	topo_store_coding_forward(r, linkGraph, nodeNum, (D-K+1), K*(D-K+1), 1, (D-K+1), sed_n_max, eMatrix, rootBuffer, &buf_rows, &sed_rows, W);
	printf("\nbuf_rows : = %d\n", buf_rows);
	printf("\nrootBuffer:\n"); jerasure_print_matrix(rootBuffer, buf_rows, K*(D-K+1), W);
	
	rootRandomMatrix = generateRandomMatrix((D-K+1), buf_rows, W);
	newComerVector = jerasure_matrix_multiply(rootRandomMatrix, rootBuffer, (D-K+1), buf_rows, buf_rows, K*(D-K+1),W);

	printf("\neMatrix:\n");jerasure_print_matrix(eMatrix, N*(D-K+1), K*(D-K+1), W);
	printf("\nnewComerVector:\n"); jerasure_print_matrix(newComerVector, (D-K+1), K*(D-K+1), W);
	i = judge_linear_independent(eMatrix, r, newComerVector, 100, N, D, K, W);
	if(i == 0) printf("\nrepair successful!\n");
	else printf("\nrepair failed!\n");

	for(i = 0; i < nodeNum; i++)
		free(linkGraph[i]);
	free(eMatrix);
	free(newComerVector);
	free(rootBuffer);
	free(rootRandomMatrix);
}

void test_topo_store_coding_forward_Jun(){
	int **linkGraph;
	int * rootBuffer;
	int * rootRandomMatrix;
	int * newComerVector;
	int *eMatrix;
	int buf_rows;
	int sed_rows;
	int nodeNum = 5;
	int r = 0;
	//int sed_n_max = D-K+1;
	int sed_n_max = 1;
	int i, j, t;
	linkGraph = Malloc(int *, nodeNum);
	for(i = 0; i < nodeNum; i++){
		linkGraph[i] = Malloc(int, nodeNum);
		memset(linkGraph[i], 0, sizeof(int)*nodeNum);
	}
	rootBuffer = Malloc(int, nodeNum*K*(D-K+1));
	
	linkGraph[0][1] = 1;
	linkGraph[1][3] = 1;
	linkGraph[1][4] = 1;
	linkGraph[4][2] = 1;
	/*
	linkGraph[0][1] = 1;
	linkGraph[0][4] = 1;
	linkGraph[1][3] = 1;
	linkGraph[4][2] = 1;
	*/
	printf("\nlinkGraph:\n");
	for(i = 0; i < nodeNum; i++){
		for(j = 0; j < nodeNum; j++)
			printf("%d ", linkGraph[i][j]);
		printf("\n");
	}
	eMatrix = Malloc(int, N*(D-K+1)*K*(D-K+1));
	topo_store_coding_forward(r, linkGraph, nodeNum, (D-K+1), K*(D-K+1), 1, (D-K+1), sed_n_max, eMatrix, rootBuffer, &buf_rows, &sed_rows, W);
	printf("\nbuf_rows : = %d\n", buf_rows);
	printf("\nrootBuffer:\n"); jerasure_print_matrix(rootBuffer, buf_rows, K*(D-K+1), W);
	
	rootRandomMatrix = generateRandomMatrix((D-K+1), buf_rows, W);
	newComerVector = jerasure_matrix_multiply(rootRandomMatrix, rootBuffer, (D-K+1), buf_rows, buf_rows, K*(D-K+1),W);

	printf("\neMatrix:\n");jerasure_print_matrix(eMatrix, N*(D-K+1), K*(D-K+1), W);
	printf("\nnewComerVector:\n"); jerasure_print_matrix(newComerVector, (D-K+1), K*(D-K+1), W);
	i = judge_linear_independent(eMatrix, r, newComerVector, 100, N, D, K, W);
	if(i == 0) printf("\nrepair successful!\n");
	else printf("\nrepair failed!\n");

	for(i = 0; i < nodeNum; i++)
		free(linkGraph[i]);
	free(eMatrix);
	free(newComerVector);
	free(rootBuffer);
	free(rootRandomMatrix);
}

void test_topo_store_coding_forward_star(){
		int **linkGraph;
	int * rootBuffer;
	int * rootRandomMatrix;
	int * newComerVector;
	int *eMatrix;
	int buf_rows;
	int sed_rows;
	int nodeNum = 5;
	int r = 0;
	int sed_n_max = D-K+1;
	//int sed_n_max = 1;
	int i, j, t;
	linkGraph = Malloc(int *, nodeNum);
	for(i = 0; i < nodeNum; i++){
		linkGraph[i] = Malloc(int, nodeNum);
		memset(linkGraph[i], 0, sizeof(int)*nodeNum);
	}
	rootBuffer = Malloc(int, nodeNum*K*(D-K+1));
	
	linkGraph[0][1] = 1;
	linkGraph[0][2] = 1;
	linkGraph[0][3] = 1;
	linkGraph[0][4] = 1;
	/*
	linkGraph[0][1] = 1;
	linkGraph[0][4] = 1;
	linkGraph[1][3] = 1;
	linkGraph[4][2] = 1;
	*/
	printf("\nlinkGraph:\n");
	for(i = 0; i < nodeNum; i++){
		for(j = 0; j < nodeNum; j++)
			printf("%d ", linkGraph[i][j]);
		printf("\n");
	}
	eMatrix = Malloc(int, N*(D-K+1)*K*(D-K+1));
	topo_store_coding_forward(r, linkGraph, nodeNum, (D-K+1), K*(D-K+1), 1, (D-K+1), sed_n_max, eMatrix, rootBuffer, &buf_rows, &sed_rows, W);
	printf("\nbuf_rows : = %d\n", buf_rows);
	printf("\nrootBuffer:\n"); jerasure_print_matrix(rootBuffer, buf_rows, K*(D-K+1), W);
	
	rootRandomMatrix = generateRandomMatrix((D-K+1), buf_rows, W);
	newComerVector = jerasure_matrix_multiply(rootRandomMatrix, rootBuffer, (D-K+1), buf_rows, buf_rows, K*(D-K+1),W);

	printf("\neMatrix:\n");jerasure_print_matrix(eMatrix, N*(D-K+1), K*(D-K+1), W);
	printf("\nnewComerVector:\n"); jerasure_print_matrix(newComerVector, (D-K+1), K*(D-K+1), W);
	i = judge_linear_independent(eMatrix, r, newComerVector, 100, N, D, K, W);
	if(i == 0) printf("\nrepair successful!\n");
	else printf("\nrepair failed!\n");

	for(i = 0; i < nodeNum; i++)
		free(linkGraph[i]);
	free(eMatrix);
	free(newComerVector);
	free(rootBuffer);
	free(rootRandomMatrix);
}

int judge_MDS_property(int r,
									  int **linkGraph,
									  int nodeNum,
									  int n,
									  int d,
									  int k,
									  int sed_n_max,
									  int w)
{
	int *rootBuffer;
	int *eMatrix;
	int *rootRandomMatrix;
	int *newComerVector;
	int buf_rows;
	int sed_rows;
	int i,j;


	srand((unsigned)time(NULL));
	rootBuffer = Malloc(int, nodeNum*k*(d-k+1));
	eMatrix = Malloc(int, n*(d-k+1)*k*(d-k+1));
	topo_store_coding_forward(r, linkGraph, nodeNum, (d-k+1), k*(d-k+1), 1, (d-k+1), sed_n_max, eMatrix, rootBuffer, &buf_rows, &sed_rows, w);
	rootRandomMatrix = generateRandomMatrix((d-k+1), buf_rows, w);
	newComerVector = jerasure_matrix_multiply(rootRandomMatrix, rootBuffer, (d-k+1), buf_rows, buf_rows, k*(d-k+1),w);
	i = judge_linear_independent(eMatrix, r, newComerVector, 100, n, d, k, w);

	free(rootBuffer);
	free(eMatrix);
	free(rootRandomMatrix);
	free(newComerVector);

	return i;
//if(i == 0) printf("\nrepair successful!\n");
//else printf("\nrepair failed!\n");
}

int judge_MDS_property2(int r,
										int *linkGraph2,
										int nodeNum,
										int n,
										int d,
										int k,
										int sed_n_max,
										int w)
{
	int *rootBuffer;
	int *eMatrix;
	int *rootRandomMatrix;
	int *newComerVector;
	int **linkGraph;
	int buf_rows;
	int sed_rows;
	int i, j;

	linkGraph = Malloc(int *, nodeNum);
	for(i = 0; i < nodeNum; i++)
		linkGraph[i] = Malloc(int, nodeNum);
	for(i = 0; i < nodeNum; i++)
		for(j = 0; j < nodeNum; j++)
			linkGraph[i][j] = linkGraph2[i*nodeNum +j];
	/*
	for(i = 0; i < nodeNum; i++)
		for(j = 0; j < nodeNum; j++)
			printf("%d ", linkGraph2[i*nodeNum+j]);
	printf("\n\n\n");

	for(i = 0; i < nodeNum; i++)
		for(j = 0; j < nodeNum; j++)
			printf("%d ", linkGraph[i][j]);
	printf("\n");
	*/

	srand((unsigned)time(NULL));
	rootBuffer = Malloc(int, nodeNum*k*(d-k+1));
	eMatrix = Malloc(int, n*(d-k+1)*k*(d-k+1));
	topo_store_coding_forward(r, linkGraph, nodeNum, (d-k+1), k*(d-k+1), 1, (d-k+1), sed_n_max, eMatrix, rootBuffer, &buf_rows, &sed_rows, w);
	rootRandomMatrix = generateRandomMatrix((d-k+1), buf_rows, w);
	newComerVector = jerasure_matrix_multiply(rootRandomMatrix, rootBuffer, (d-k+1), buf_rows, buf_rows, k*(d-k+1),w);
	i = judge_linear_independent(eMatrix, r, newComerVector, 100, n, d, k, w);

	free(rootBuffer);
	free(eMatrix);
	free(rootRandomMatrix);
	free(newComerVector);
	for(j = 0; j < nodeNum; j++)
		free(linkGraph[j]);
	free(linkGraph);
	return i;
	//if(i == 0) printf("\nrepair successful!\n");
	//else printf("\nrepair failed!\n");
}

void test_file_retrievable( FILE *fp,
										int r2,
										int **linkGraph2,
										int nodeNum,
										int n,
										int d,
										int k,
										int sed_n_max,
										int w,
										int testNum,
										int eachtestTimes)
{
	int *rootBuffer;
	int *eMatrix;
	int *rootRandomMatrix;
	int *newComerVector;
	int **linkGraph;
	int buf_rows;
	int sed_rows;
	int i, j;
	int r;
	int sucRevNum;

	srand((unsigned)time(NULL));
	//eMatrix = Malloc(int, n*(d-k+1)*k*(d-k+1));
	eMatrix = generateRandomMatrix(n*(d-k+1), k*(d-k+1), w);
	
	sucRevNum = judge_retrievable_original_file(eMatrix, -1, NULL,eachtestTimes, n, d, k, w);
	printf("sucRevNum = %d\n", sucRevNum);
	fprintf(fp, "%d\t", sucRevNum);

	while(testNum--){
		if(r2 == -1)
			r = rand()%nodeNum;

		linkGraph = Malloc(int *, nodeNum);
		for(i = 0; i < nodeNum; i++)
			linkGraph[i] = Malloc(int, nodeNum);
		for(i = 0; i < nodeNum; i++)
			for(j = 0; j < nodeNum; j++)
				linkGraph[i][j] = linkGraph2[r][i*nodeNum +j];

//		jerasure_print_matrix(eMatrix, n*(d-k+1), k*(d-k+1),w);printf("\n");
		sucRevNum = 0;
		rootBuffer = Malloc(int, nodeNum*k*(d-k+1));
		topo_store_coding_forward2(r, linkGraph, nodeNum, (d-k+1), k*(d-k+1), 1, (d-k+1), sed_n_max, eMatrix, rootBuffer, &buf_rows, &sed_rows, w);
		
//		jerasure_print_matrix(eMatrix, n*(d-k+1), k*(d-k+1),w);printf("\n");
//		printf("\nrootbuffer\n"); jerasure_print_matrix(rootBuffer, buf_rows, k*(d-k+1), w);
		rootRandomMatrix = generateRandomMatrix((d-k+1), buf_rows, w);
//		printf("\nrootRandomMatrix:\n");  jerasure_print_matrix(rootRandomMatrix, (d-k+1), buf_rows, w);
 		newComerVector = jerasure_matrix_multiply(rootRandomMatrix, rootBuffer, (d-k+1), buf_rows, buf_rows, k*(d-k+1),w);
//		printf("\nnewComerVector\n"); jerasure_print_matrix(newComerVector, (d-k+1), k*(d-k+1),w);
		sucRevNum = judge_retrievable_original_file(eMatrix, r, newComerVector,eachtestTimes, n, d, k, w);
		printf("sucRevNum = %d\n", sucRevNum);
		fprintf(fp, "%d\t", sucRevNum);
		free(rootRandomMatrix);
		free(newComerVector);
		free(rootBuffer);

		for(j = 0; j < nodeNum; j++)
			free(linkGraph[j]);
		free(linkGraph);
	}
	free(eMatrix);

}

void test_file_retrievable_W( FILE *fp,
	int *eMatrix,
	int r,
	int *linkGraph2,
	int nodeNum, //is n
	int n,
	int d,
	int k,
	int sed_n_max,
	int w,
	int eachtestTimes)
{
	int *rootBuffer;
	
	int *rootRandomMatrix;
	int *newComerVector;
	int **linkGraph;
	int buf_rows;
	int sed_rows;
	int i, j;
	int sucRevNum;

	srand((unsigned)time(NULL));
	//eMatrix = Malloc(int, n*(d-k+1)*k*(d-k+1));

//	sucRevNum = judge_retrievable_original_file(eMatrix, -1, NULL,eachtestTimes, n, d, k, w);
//	printf("sucRevNum = %d\n", sucRevNum);
//	fprintf(fp, "%d\t", sucRevNum);


		linkGraph = Malloc(int *, n);
		for(i = 0; i < n; i++)
			linkGraph[i] = Malloc(int, n);
		for(i = 0; i < n; i++)
			for(j = 0; j < n; j++)
				linkGraph[i][j] = linkGraph2[i*n +j];

		//		jerasure_print_matrix(eMatrix, n*(d-k+1), k*(d-k+1),w);printf("\n");
		sucRevNum = 0;
		rootBuffer = Malloc(int, nodeNum*k*(d-k+1));
		topo_store_coding_forward2(r, linkGraph, nodeNum, (d-k+1), k*(d-k+1), 1, (d-k+1), sed_n_max, eMatrix, rootBuffer, &buf_rows, &sed_rows, w);

		//		jerasure_print_matrix(eMatrix, n*(d-k+1), k*(d-k+1),w);printf("\n");
		//		printf("\nrootbuffer\n"); jerasure_print_matrix(rootBuffer, buf_rows, k*(d-k+1), w);
		rootRandomMatrix = generateRandomMatrix((d-k+1), buf_rows, w);
		//		printf("\nrootRandomMatrix:\n");  jerasure_print_matrix(rootRandomMatrix, (d-k+1), buf_rows, w);
		newComerVector = jerasure_matrix_multiply(rootRandomMatrix, rootBuffer, (d-k+1), buf_rows, buf_rows, k*(d-k+1),w);
		//		printf("\nnewComerVector\n"); jerasure_print_matrix(newComerVector, (d-k+1), k*(d-k+1),w);
		sucRevNum = judge_retrievable_original_file(eMatrix, r, newComerVector,eachtestTimes, n, d, k, w);
		printf("sucRevNum = %d\n", sucRevNum);
		fprintf(fp, "%d\t", sucRevNum);
		free(rootRandomMatrix);
		free(newComerVector);
		free(rootBuffer);

		for(j = 0; j < n; j++)
			free(linkGraph[j]);
		free(linkGraph);
	//free(eMatrix);

}

void test_file_retrievable_W_block_num( FILE *fp,
	int *eMatrix,
	int r,
	int *linkGraph2,
	int nodeNum, //is n
	int n,
	int d,
	int k,
	int sed_n_max,
	int w,
	int intertestTime,
	int outerTestTimes)
{
	int *rootBuffer;

	int *rootRandomMatrix;
	int *newComerVector;
	int **linkGraph;
	int buf_rows;
	int sed_rows;
	int i, j;
	int avgNeedNum = 0;

	srand((unsigned)time(NULL));
	//eMatrix = Malloc(int, n*(d-k+1)*k*(d-k+1));

	//	sucRevNum = judge_retrievable_original_file(eMatrix, -1, NULL,eachtestTimes, n, d, k, w);
	//	printf("sucRevNum = %d\n", sucRevNum);
	//	fprintf(fp, "%d\t", sucRevNum);


	linkGraph = Malloc(int *, n);
	for(i = 0; i < n; i++)
		linkGraph[i] = Malloc(int, n);
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			linkGraph[i][j] = linkGraph2[i*n +j];

	//		jerasure_print_matrix(eMatrix, n*(d-k+1), k*(d-k+1),w);printf("\n");
	avgNeedNum = 0;
	rootBuffer = Malloc(int, nodeNum*k*(d-k+1));
	topo_store_coding_forward2(r, linkGraph, nodeNum, (d-k+1), k*(d-k+1), 1, (d-k+1), sed_n_max, eMatrix, rootBuffer, &buf_rows, &sed_rows, w);

	//		jerasure_print_matrix(eMatrix, n*(d-k+1), k*(d-k+1),w);printf("\n");
	//		printf("\nrootbuffer\n"); jerasure_print_matrix(rootBuffer, buf_rows, k*(d-k+1), w);
	rootRandomMatrix = generateRandomMatrix((d-k+1), buf_rows, w);
	//		printf("\nrootRandomMatrix:\n");  jerasure_print_matrix(rootRandomMatrix, (d-k+1), buf_rows, w);
	newComerVector = jerasure_matrix_multiply(rootRandomMatrix, rootBuffer, (d-k+1), buf_rows, buf_rows, k*(d-k+1),w);
	j = outerTestTimes;
	while(j--){
		avgNeedNum += linear_independent_num2(eMatrix, r, newComerVector, n, d, k, intertestTime, w);
	}
	avgNeedNum /= outerTestTimes;
	printf("avgNeedNum = %d\n", avgNeedNum);
	fprintf(fp, "%d\t", avgNeedNum);
	//		printf("\nnewComerVector\n"); jerasure_print_matrix(newComerVector, (d-k+1), k*(d-k+1),w);
	//sucRevNum = judge_retrievable_original_file(eMatrix, r, newComerVector,eachtestTimes, n, d, k, w);
//	printf("sucRevNum = %d\n", sucRevNum);
//	fprintf(fp, "%d\t", sucRevNum);
	free(rootRandomMatrix);
	free(newComerVector);
	free(rootBuffer);

	for(j = 0; j < n; j++)
		free(linkGraph[j]);
	free(linkGraph);
	//free(eMatrix);

}


void test1(){
		//srand((unsigned)time(NULL));
	
	//test_rc_msr_repair_proc_one_time();
	//test_rc_msr_repiar_proc_multi_times();
	//test_store_coding_forward();
	//test_store_coding_forward2();
	//test_topo_store_coding_forward();
	//test_topo_store_coding_forward_Jun();
	//test_topo_store_coding_forward_star();
	int r = 0;
	int **linkGraph;
	int nodeNum = 5;
	int sed_n_max = D-K+1;
	//int sed_n_max = 1;
	int i, j, k;
	
	linkGraph = Malloc(int *, nodeNum);
	for(i = 0; i < nodeNum; i++){
		linkGraph[i] = Malloc(int, nodeNum);
		memset(linkGraph[i], 0, sizeof(int)*nodeNum);
	}
	
	linkGraph[0][1] = 1;
	linkGraph[1][3] = 1;
	linkGraph[1][4] = 1;
	linkGraph[4][2] = 1;
	/*
	linkGraph[0][1] = 1;
	linkGraph[0][4] = 1;
	linkGraph[1][3] = 1;
	linkGraph[4][2] = 1;
	*/
	i = judge_MDS_property(r, linkGraph, nodeNum, N, D, K, sed_n_max, W);
	if(i == 0) printf("Yes!\n");
	else printf("No!\n");
	for(i = 0; i < nodeNum; i++)
		free(linkGraph[i]);
	free(linkGraph);
}

void test2(){
		//srand((unsigned)time(NULL));
	
	//test_rc_msr_repair_proc_one_time();
	//test_rc_msr_repiar_proc_multi_times();
	//test_store_coding_forward();
	//test_store_coding_forward2();
	//test_topo_store_coding_forward();
	//test_topo_store_coding_forward_Jun();
	//test_topo_store_coding_forward_star();
	int r = 0;
	int *linkGraph;
	int nodeNum = 5;
	int sed_n_max = D-K+1;
	//int sed_n_max = 1;
	int i, j, k;
	
	linkGraph = Malloc(int, nodeNum*nodeNum);
	memset(linkGraph, 0, sizeof(int)*nodeNum*nodeNum);
	
	linkGraph[0*nodeNum+1] = 1;
	linkGraph[1*nodeNum+3] = 1;
	linkGraph[1*nodeNum+4] = 1;
	linkGraph[4*nodeNum+2] = 1;
	//jerasure_print_matrix(linkGraph, nodeNum, nodeNum, W);
	/*
	linkGraph[0][1] = 1;
	linkGraph[0][4] = 1;
	linkGraph[1][3] = 1;
	linkGraph[4][2] = 1;
	*/
	i = judge_MDS_property2(r, linkGraph, nodeNum, N, D, K, sed_n_max, W);
	if(i == 0) printf("Yes!\n");
	else printf("No!\n");
	free(linkGraph);
}

void test_one_repiar_mds_from_file(){
	char filename[] = "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\JunSingleMaxSpanningTreeLinkGraph.txt";
//	char filename[] = "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\maxSpannTreeLinkGraph.txt";
//	char filename[] = "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\starLinkGraph.txt";
	FILE* fp;
	int i;
	int j, j2, j3;
	int n = 50;
	int d = n-1;
	int k = n-2;
	//int sed_n_max = d-k+1;
	int sed_n_max = 1;
	int w = 16;
	int r = 0;
	
	int nodeNum;
	int *linkGraph;
	char StrLine[MAX_LENGTH];
	char c[100];
	int totalSucNum = 0;
	if((fp = fopen(filename, "r")) == NULL){
		printf("error!");
		return;
	}
	nodeNum = d+1;
	linkGraph = Malloc(int, nodeNum*nodeNum);
	//printf("nodeNum = %d\n", nodeNum);
	while(!feof(fp)){
		fgets(StrLine, MAX_LENGTH, fp);		
		j2 = 0;
		j3 = 0;
	//	printf("%s\n\n\n", StrLine);
		for(j = 0; StrLine[j] != '\0'; j++){
			if(StrLine[j] != ' ')
				c[j2++] = StrLine[j];
			else{
				c[j2++] = '\0';
				linkGraph[j3++] = atoi(c);
				
				j2 = 0;
			}
		}
//		for(j = 0; j < nodeNum*nodeNum; j++)
//			printf("%d ", linkGraph[j]);
//		printf("\n"); 
		/*
		for(i = 0; i < nodeNum; i++)
			for(j = 0; j < nodeNum; j++)
				printf("%d ", linkGraph[i*nodeNum+j]);
		printf("\n\n\n");
		*/
		i = judge_MDS_property2(r, linkGraph, nodeNum, n, d, k, sed_n_max, w);
		if(i == 0) {printf("Yes!\n"); totalSucNum++;}
		else printf("No!\n");
	}
	free(linkGraph);
	printf("totalSucNum = %d\n", totalSucNum);
	fclose(fp);
}


void test_one_tree_file_retrievable(){
	char filename[] = "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\JunEachNodeAsRootLinkGraph.txt";
//	char filename[] = "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\maxSpannTreeLinkGraph.txt";
//	char filename[] = "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\starLinkGraph.txt";
	char outputfilename[] =  "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\JunEachNodeAsRootLinkGraphoutput.txt";
	FILE* fp;
	int i;
	int j, j2, j3;
	int n = 5;
	int d = n-1;
	int k = n-2;
	//int sed_n_max = d-k+1;
	int sed_n_max = 1;
	int w = 16;
	int r = -1;  //random node failed
	
	int nodeNum;
	int **linkGraph;
	char StrLine[MAX_LENGTH];
	char c[100];
	int totalSucNum = 0;
	if((fp = fopen(filename, "r")) == NULL){
		printf("error!");
		return;
	}
	nodeNum = d+1;
	linkGraph = Malloc(int *, nodeNum+1);
	for(i = 0; i < nodeNum+1; i++)
		linkGraph[i] = Malloc(int, nodeNum*nodeNum);
	
	i = 0;
	while(!feof(fp))
	{
		fgets(StrLine, MAX_LENGTH, fp);		
		j2 = 0;
		j3 = 0;
		for(j = 0; StrLine[j] != '\0'; j++){
			if(StrLine[j] != ' ')
				c[j2++] = StrLine[j];
			else{
				c[j2++] = '\0';
				linkGraph[i][j3++] = atoi(c);
				
				j2 = 0;
			}
		}
		i++;
	}
	fclose(fp);
	/*
	for(i = 0; i < nodeNum+1; i++){
		for(j = 0; j < nodeNum*nodeNum; j++)
			printf("%d ", linkGraph[i][j]);
		printf("\n");
	}
	*/

	if((fp = fopen(outputfilename, "w")) == NULL){
		printf("error!");
		return;
	}

	test_file_retrievable(fp, r,linkGraph,nodeNum, n, d,k, sed_n_max, w, 100, 10000);
	fclose(fp);
	for(i = 0; i < nodeNum+1; i++)
		free(linkGraph[i]);
	free(linkGraph);
	
}

void test_avg_trees_file_retrievable(){
	char filename[MAX_LENGTH2];
//	char filename[] = "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\JunEachNodeAsRootLinkGraph.txt";
//	char filename[] = "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\maxSpannTreeLinkGraph.txt";
//	char filename[] = "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\starLinkGraph.txt";
	char outputfilename[] =  "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\JunEachNodeAsRootLinkGraphoutput.txt";
	FILE* fp, *fp2;
	int i, i2;
	int j, j2, j3;
	int n = 10;
	int d = n-1;
	int k = n-2;
	//int sed_n_max = d-k+1;
	int sed_n_max = 1;
	int w = 16;
	int r = -1;  //random node failed
	int avgNum = 100;
	int nodeNum;
	int **linkGraph;
	char StrLine[MAX_LENGTH];
	char c[100];
	int totalSucNum = 0;
	
	if((fp2 = fopen(outputfilename, "w")) == NULL){
		printf("error!");
		return;
	}

	for(i2 = 0; i2 < avgNum; i2++){
		sprintf(filename, "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\avg\\JunEachNodeAsRootLinkGraph_%d.txt", i2);
		printf("%s\n", filename);
		if((fp = fopen(filename, "r")) == NULL){
			printf("error!\n");
			return;
		}
		nodeNum = d+1;
		linkGraph = Malloc(int *, nodeNum+1);
		for(i = 0; i < nodeNum+1; i++)
			linkGraph[i] = Malloc(int, nodeNum*nodeNum);
	
		i = 0;
		while(!feof(fp))
		{
			fgets(StrLine, MAX_LENGTH, fp);		
			j2 = 0;
			j3 = 0;
			for(j = 0; StrLine[j] != '\0'; j++){
				if(StrLine[j] != ' ')
					c[j2++] = StrLine[j];
				else{
					c[j2++] = '\0';
					linkGraph[i][j3++] = atoi(c);
				
					j2 = 0;
				}
			}
			i++;
		}
		fclose(fp);

		test_file_retrievable(fp2, r,linkGraph,nodeNum, n, d,k, sed_n_max, w, 50, 1000);
		fprintf(fp2, "\n");
		for(i = 0; i < nodeNum+1; i++)
			free(linkGraph[i]);
		free(linkGraph);
	}
	fclose(fp2);
}

void test_avg_trees_file_retrievable_random_d_W(){
		char filename[MAX_LENGTH2];
		char filename2[MAX_LENGTH2];
		char outputfilename[] =  "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\JunEachNodeAsRootLinkGraphRandom_d_output.txt";
		FILE* fp, *fp1, *fp2;
		int i, i2, ii;
		int j, j2, j3, j4;
		int n = 11;
		int d = 8;
		int k = 6;
		//int sed_n_max = d-k+1;
		int sed_n_max = 1;
		int w = 16;
		int r;
		int totalTopoNum = 100;//
		int repairRounds = 100;
		int eachTestTimes = 100;
		int nodeNum;
		int *eMatrix;
		int *linkGraph;
		char StrLine[MAX_LENGTH];
		char StrLine2[MAX_LENGTH];
		char c[100];
		char c2[100];
		int totalSucNum = 0;

		if((fp2 = fopen(outputfilename, "w")) == NULL){
			printf("error!");
			return;
		}
		for(ii = 0; ii < totalTopoNum; ii++){
			sprintf(filename, "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\linkGraph\\JunEachNodeAsRootLinkGraphRandom_d_%d.txt", ii);
			printf("%s\n", filename);
			if((fp = fopen(filename, "r")) == NULL){
				printf("error!\n");
				return;
			}

			sprintf(filename2, "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\root\\JunEachNodeAsRootLinkGraphRandom_d_root_%d.txt", ii);
			printf("%s\n", filename2);
			if((fp1 = fopen(filename2, "r")) == NULL){
				printf("error!\n");
				return;
			}

			nodeNum = d+1;
			linkGraph = Malloc(int, n*n);

			eMatrix = generateRandomMatrix(n*(d-k+1), k*(d-k+1), w);

			while(!feof(fp) && !feof(fp1))
			{
				fgets(StrLine, MAX_LENGTH, fp);
				fgets(StrLine2, 10, fp1);

	//			printf("%s\n", StrLine);
	//			printf("%s\n", StrLine2);

				j2 = 0;
				j3 = 0;
				for(j = 0; StrLine[j] != '\0'; j++){
					if(StrLine[j] != ' ' )
					{
						c[j2] = StrLine[j];
						j2++;
					}
					else{
						c[j2] = '\0';
						linkGraph[j3++] = atoi(c);
						j2 = 0;
					}
				}


				r = atoi(StrLine2);
				/*
				j4 = 0;
			
				for(j = 0; StrLine2[j] != '\0'; j++){
					if(StrLine2[j] != ' ' )
					{
						c2[j4] = StrLine2[j];
						j4++;
					}
					else{
						c2[j4] = '\0';
						r = atoi(c2);
						j4 = 0;
					}
				}
				*/
	//			jerasure_print_matrix(linkGraph,n,n,w);
	//			printf("r = %d\n", r);
				test_file_retrievable_W(fp2, eMatrix, r,linkGraph,n, n, d, k, sed_n_max, w,eachTestTimes);
			}
			fprintf(fp2, "\n");
			free(linkGraph);
			free(eMatrix);
			fclose(fp1);
			fclose(fp);
		}
	fclose(fp2);

}

void test_avg_trees_file_retrievable_random_d_W_block_num(){
		char filename[MAX_LENGTH2];
		char filename2[MAX_LENGTH2];
		char outputfilename[] =  "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\JunEachNodeAsRootLinkGraphRandom_d_output.txt";
		FILE* fp, *fp1, *fp2;
		int i, i2, ii;
		int j, j2, j3, j4;
		int n = 10;
		int d = 8;
		int k = 6;
		//int sed_n_max = d-k+1;
		int sed_n_max = 1;
		int w = 16;
		int r;
		int totalTopoNum = 10;//
		//int repairRounds = 30;

		int interTestTime = 5;//内圈循环，由于测试能够恢复源文件，不宜过大。
		int outerTestTime = 100; //外圈循环，用来平均的
		int nodeNum;
		int *eMatrix;
		int *linkGraph;
		char StrLine[MAX_LENGTH];
		char StrLine2[MAX_LENGTH];
		char c[100];
		char c2[100];
		int totalSucNum = 0;

		if((fp2 = fopen(outputfilename, "w")) == NULL){
			printf("error!");
			return;
		}
		for(ii = 0; ii < totalTopoNum; ii++){
			sprintf(filename, "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\linkGraph\\JunEachNodeAsRootLinkGraphRandom_d_%d.txt", ii);
			printf("%s\n", filename);
			if((fp = fopen(filename, "r")) == NULL){
				printf("error!\n");
				return;
			}

			sprintf(filename2, "D:\\Dropbox\\wds\\yanwang\\infocom14_wy\\result\\root\\JunEachNodeAsRootLinkGraphRandom_d_root_%d.txt", ii);
			printf("%s\n", filename2);
			if((fp1 = fopen(filename2, "r")) == NULL){
				printf("error!\n");
				return;
			}

			nodeNum = d+1;
			linkGraph = Malloc(int, n*n);

			eMatrix = generateRandomMatrix(n*(d-k+1), k*(d-k+1), w);

			while(!feof(fp) && !feof(fp1))
			{
				fgets(StrLine, MAX_LENGTH, fp);
				fgets(StrLine2, 10, fp1);

	//			printf("%s\n", StrLine);
	//			printf("%s\n", StrLine2);

				j2 = 0;
				j3 = 0;
				for(j = 0; StrLine[j] != '\0'; j++){
					if(StrLine[j] != ' ' )
					{
						c[j2] = StrLine[j];
						j2++;
					}
					else{
						c[j2] = '\0';
						linkGraph[j3++] = atoi(c);
						j2 = 0;
					}
				}


				r = atoi(StrLine2);
				/*
				j4 = 0;
			
				for(j = 0; StrLine2[j] != '\0'; j++){
					if(StrLine2[j] != ' ' )
					{
						c2[j4] = StrLine2[j];
						j4++;
					}
					else{
						c2[j4] = '\0';
						r = atoi(c2);
						j4 = 0;
					}
				}
				*/
	//			jerasure_print_matrix(linkGraph,n,n,w);
	//			printf("r = %d\n", r);
				//test_file_retrievable_W(fp2, eMatrix, r,linkGraph,n, n, d, k, sed_n_max, w,eachTestTimes);
				test_file_retrievable_W_block_num(fp2, eMatrix, r, linkGraph, n, n, d, k, sed_n_max, w,interTestTime, outerTestTime);
			}
			fprintf(fp2, "\n");
			free(linkGraph);
			free(eMatrix);
			fclose(fp1);
			fclose(fp);
		}
	fclose(fp2);

}