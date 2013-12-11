#ifndef _RC_MSR_TOOLS_H
#define _RC_MSR_TOOLS_H

#include <stdio.h>
#define Malloc(type, n) (type *)malloc((n)*sizeof(type))

int * generateRandomMatrix(int rows, int cols, int w);
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
	int w);

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
	int w);

int * topo_store_coding_forward(int r, 
	int ** linkGraph,
	int nodeNum,
	int r1,
	int c1,
	int r2,
	int c2,
	int sed_n_max,
	int *eMatrix,
	int *buf,
	int *buf_rows,
	int w
	);
int judge_linear_independent(int *eMatrix, 
	int t,
	int *newComerVector,
	int testTimes,
	int n,
	int d,
	int k,
	int w);
int judge_retrievable_original_file(int *eMatrix,
	int t,
	int *newComerVector,
	int testTimes,
	int n,
	int d,
	int k,
	int w);
int linear_independent_num(int *eMatrix,
	int n,
	int d,
	int k,
	int avgTimes,
	//											int *select,
	//											int select_num,
	int w
	);
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
	);
int random_select(int *a,  int *b, int *c, int n_a , int n_b, int n_c );
void setMatrixValue(int *matrix, int rows, int cols, int r1, int r2, int c1, int c2, int value);
void printMatrix(int *matrix, int rows, int cols, int r1, int r2, int c1, int c2);
void assignMatrixToMatrix(int *matrix, int rows, int cols, int r1, int r2, int c1, int c2, int * matrix2, int r3, int r4, int c3, int c4);
void test_rc_msr_repair_proc_one_time();
void test_rc_msr_repiar_proc_multi_times();
void test_store_coding_forward();
void test_store_coding_forward2();
void test_topo_store_coding_forward();
void test_topo_store_coding_forward_Jun();
void test_topo_store_coding_forward_star();
int judge_MDS_property(int r,
	int **linkGraph,
	int nodeNum,
	int n,
	int d,
	int k,
	int sed_n_max,
	int w);

int judge_MDS_property2(int r,
	int *linkGraph2,
	int nodeNum,
	int n,
	int d,
	int k,
	int sed_n_max,
	int w);
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
	int eachtestTimes);
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
	int eachtestTimes);
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
	int outerTestTimes);
void test1();
void test2();
void test_one_repiar_mds_from_file();
void test_one_repiar_mds_from_file_debug();

void test_one_tree_file_retrievable();
void test_avg_trees_file_retrievable();
void test_avg_trees_file_retrievable_random_d_W(); //用来平均求出Jun的工作里随着修复次数的增多，用户无法恢复文件的概率
void test_avg_trees_file_retrievable_random_d_W_block_num();  //用来平均求出Jun的工作里随着修复次数的增多，用户如果想恢复文件需要多少的数据块（如果大于n，则表示不能恢复）
#endif