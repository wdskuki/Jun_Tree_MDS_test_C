/* jerasure.h - header of kernel procedures
 * James S. Plank

JERASURE - Library for Erasure Coding
Copright (C) 2007 James S. Plank

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

James S. Plank
Department of Electrical Engineering and Computer Science
University of Tennessee
Knoxville, TN 37996
plank@cs.utk.edu
*/

/*
 * $Revision: 1.2 $
 * $Date: 2008/08/19 17:40:58 $
 */

#ifndef _JERASURE_H
#define _JERASURE_H

/* This uses procedures from the Galois Field arithmetic library */

#include "galois.h"


/*
int *jerasure_matrix_to_bitmatrix(int k, int m, int w, int *matrix);
int **jerasure_dumb_bitmatrix_to_schedule(int k, int m, int w, int *bitmatrix);
int **jerasure_smart_bitmatrix_to_schedule(int k, int m, int w, int *bitmatrix);
int ***jerasure_generate_schedule_cache(int k, int m, int w, int *bitmatrix, int smart);

void jerasure_free_schedule(int **schedule);
void jerasure_free_schedule_cache(int k, int m, int ***cache);




void jerasure_do_parity(int k, char **data_ptrs, char *parity_ptr, int size);

void jerasure_matrix_encode(int k, int m, int w, int *matrix,
                          char **data_ptrs, char **coding_ptrs, int size);

void jerasure_bitmatrix_encode(int k, int m, int w, int *bitmatrix,
                            char **data_ptrs, char **coding_ptrs, int size, int packetsize);

void jerasure_schedule_encode(int k, int m, int w, int **schedule,
                                  char **data_ptrs, char **coding_ptrs, int size, int packetsize);



int jerasure_matrix_decode(int k, int m, int w, 
                          int *matrix, int row_k_ones, int *erasures,
                          char **data_ptrs, char **coding_ptrs, int size);
                          
int jerasure_bitmatrix_decode(int k, int m, int w, 
                            int *bitmatrix, int row_k_ones, int *erasures,
                            char **data_ptrs, char **coding_ptrs, int size, int packetsize);

int jerasure_schedule_decode_lazy(int k, int m, int w, int *bitmatrix, int *erasures,
                            char **data_ptrs, char **coding_ptrs, int size, int packetsize,
                            int smart);

int jerasure_schedule_decode_cache(int k, int m, int w, int ***scache, int *erasures,
                            char **data_ptrs, char **coding_ptrs, int size, int packetsize);

int jerasure_make_decoding_matrix(int k, int m, int w, int *matrix, int *erased, 
                                  int *decoding_matrix, int *dm_ids);

int jerasure_make_decoding_bitmatrix(int k, int m, int w, int *matrix, int *erased, 
                                  int *decoding_matrix, int *dm_ids);

int *jerasure_erasures_to_erased(int k, int m, int *erasures);

 
void jerasure_matrix_dotprod(int k, int w, int *matrix_row,
                          int *src_ids, int dest_id,
                          char **data_ptrs, char **coding_ptrs, int size);

void jerasure_bitmatrix_dotprod(int k, int w, int *bitmatrix_row,
                             int *src_ids, int dest_id,
                             char **data_ptrs, char **coding_ptrs, int size, int packetsize);

void jerasure_do_scheduled_operations(char **ptrs, int **schedule, int packetsize);


int jerasure_invert_matrix(int *mat, int *inv, int rows, int w);
int jerasure_invert_bitmatrix(int *mat, int *inv, int rows);
int jerasure_invertible_matrix(int *mat, int rows, int w);
int jerasure_invertible_bitmatrix(int *mat, int rows);


void jerasure_print_matrix(int *matrix, int rows, int cols, int w);
void jerasure_print_bitmatrix(int *matrix, int rows, int cols, int w);


int *jerasure_matrix_multiply(int *m1, int *m2, int r1, int c1, int r2, int c2, int w);


void jerasure_get_stats(double *fill_in);

*/

void jerasure_print_matrix(int *m, int rows, int cols, int w);
void jerasure_print_bitmatrix(int *m, int rows, int cols, int w);
int jerasure_make_decoding_matrix(int k, int m, int w, int *matrix, int *erased, int *decoding_matrix, int *dm_ids);
int jerasure_make_decoding_bitmatrix(int k, int m, int w, int *matrix, int *erased, int *decoding_matrix, int *dm_ids);
int jerasure_matrix_decode(int k, int m, int w, int *matrix, int row_k_ones, int *erasures,
	char **data_ptrs, char **coding_ptrs, int size);
int *jerasure_matrix_to_bitmatrix(int k, int m, int w, int *matrix);
void jerasure_matrix_encode(int k, int m, int w, int *matrix,char **data_ptrs, char **coding_ptrs, int size);
void jerasure_bitmatrix_dotprod(int k, int w, int *bitmatrix_row,
	int *src_ids, int dest_id,
	char **data_ptrs, char **coding_ptrs, int size, int packetsize);
void jerasure_do_parity(int k, char **data_ptrs, char *parity_ptr, int size);
int jerasure_invert_matrix(int *mat, int *inv, int rows, int w);
int jerasure_invertible_matrix(int *mat, int rows, int w);
int *jerasure_erasures_to_erased(int k, int m, int *erasures);
void jerasure_free_schedule(int **schedule);
void jerasure_free_schedule_cache(int k, int m, int ***cache);
void jerasure_matrix_dotprod(int k, int w, int *matrix_row,
	int *src_ids, int dest_id,
	char **data_ptrs, char **coding_ptrs, int size);
int jerasure_bitmatrix_decode(int k, int m, int w, int *bitmatrix, int row_k_ones, int *erasures,
	char **data_ptrs, char **coding_ptrs, int size, int packetsize);
static char **set_up_ptrs_for_scheduled_decoding(int k, int m, int *erasures, char **data_ptrs, char **coding_ptrs);
static int set_up_ids_for_scheduled_decoding(int k, int m, int *erasures, int *row_ids, int *ind_to_row);
static int **jerasure_generate_decoding_schedule(int k, int m, int w, int *bitmatrix, int *erasures, int smart);
int jerasure_schedule_decode_lazy(int k, int m, int w, int *bitmatrix, int *erasures,
	char **data_ptrs, char **coding_ptrs, int size, int packetsize,int smart);
int jerasure_schedule_decode_cache(int k, int m, int w, int ***scache, int *erasures,char **data_ptrs, char **coding_ptrs, int size, int packetsize);
int ***jerasure_generate_schedule_cache(int k, int m, int w, int *bitmatrix, int smart);
int jerasure_invert_bitmatrix(int *mat, int *inv, int rows);
int jerasure_invertible_bitmatrix(int *mat, int rows);
int *jerasure_matrix_multiply(int *m1, int *m2, int r1, int c1, int r2, int c2, int w);
void jerasure_get_stats(double *fill_in);
void jerasure_get_stats(double *fill_in);
void jerasure_do_scheduled_operations(char **ptrs, int **operations, int packetsize);
void jerasure_schedule_encode(int k, int m, int w, int **schedule,
	char **data_ptrs, char **coding_ptrs, int size, int packetsize);
int **jerasure_dumb_bitmatrix_to_schedule(int k, int m, int w, int *bitmatrix);
int **jerasure_smart_bitmatrix_to_schedule(int k, int m, int w, int *bitmatrix);
void jerasure_bitmatrix_encode(int k, int m, int w, int *bitmatrix,char **data_ptrs, char **coding_ptrs, int size, int packetsize);

#endif
