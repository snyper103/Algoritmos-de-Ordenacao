#ifndef STDS_H_INCLUDED
#define STDS_H_INCLUDED

void selectionSort ( long* array, long size );
void bubbleSort ( long* array, long size );
void recursiveBubbleSort ( long* array, long size );
void insertionSort ( long* array, long size );
void recursiveInsertionSort ( long* array, long size );
void mergeSort ( long* array, long left_index, long right_index );
void quickSort ( long* array, long low_index, long high_index );
void heapSort ( long* array, long size );
void shellSort ( long* array, long size );
void countingSort ( long* array, long size );
void radixSort ( long* array, long size );
long exemplo (long** matriz, long coluna, long linha, long* count);

#endif
