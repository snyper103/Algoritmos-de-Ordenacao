#include "stds.h"
#include <stdlib.h>

/** just in the .c file. */
void SWAP ( long* a, long* b )
{
    long temp = *a;
    *a = *b;
    *b = temp;
}

void mergeArray ( long* array, long left_index, long middle_index, long right_index )
{
    register long i, j, l;
    long left_size, right_size,* left_array,* right_array;

    left_size = middle_index - left_index + 1; /** Calculates the size of the left array. */
    right_size = right_index - middle_index; /** Calculates the size of the right array. */

    /* Creates the left and right auxiliary arrays. */
    left_array = (long*)malloc(left_size * sizeof(long));   /** Makes the left array */
    right_array = (long*)malloc(right_size * sizeof(long)); /** Makes the right array */

    /* Copy data to the auxiliary arrays, left and right. */
    for ( i = 0; i < left_size; ++i )
        left_array[i] = array[left_index+i];    /** Make the copy of what is on the left of original array to the left auxiliary array. */
    for ( i = 0; i < right_size; ++i )
        right_array[i] = array[middle_index+1+i]; /** Make the copy of what is on the left of original array to the right auxiliary array. */

    /* Merge the auxiliary arrays back into the original array. */
    for ( i = 0, j = 0, l = left_index; (i < left_size) && (j < right_size); ++l ) {
        if ( left_array[i] <= right_array[j] )
        {
            array[l] = left_array[i];
            ++i;
        }

        else
        {
            array[l] = right_array[j];
            ++j;
        }
    }

    /* Copy the remaining elements of left_array, if there are any */
    for ( ; i < left_size; ++l, ++i )
        array[l] = left_array[i];

    /* Copy the remaining elements of right_array, if there are any */
    for ( ; j < right_size; ++l, ++j )
        array[l] = right_array[j];

    free(left_array);
    free(right_array);
}

/*  This function takes last element as pivot, places the pivot element
    at its correct position in sorted array, and places all smaller
    (smaller than pivot) to left of pivot and all greater elements to
    right of pivot  */
long partitionArray ( long* array, long low_index, long high_index )
{
    register long i, j;
    long pivot = array[high_index];

    for ( i = (low_index - 1), j = low_index; j < high_index; ++j ) {
        if ( array[j] <= pivot )
        {
            ++i;
            SWAP(&array[i], &array[j]);
        }
    }

    SWAP(&array[i+1], &array[high_index]);

    return i + 1;
}

/*  here is where we build the heap and the binary tree. */
void heapify ( long* array, long size, long index )
{
    /* Initialize largest as root, left = 2*i + 1, right = 2*i + 2 */
    long left_index = index * 2 + 1, right_index = index * 2 + 2, largest_index = index;

    /* If left child is larger than root */
    if ( (left_index < size) && (array[left_index] > array[largest_index]) )
        largest_index = left_index;

    /* If right child is larger than largest so far */
    if ( (right_index < size) && (array[right_index] > array[largest_index]) )
        largest_index = right_index;

    /* If largest is not root */
    if ( largest_index != index )
    {
        SWAP(&array[index], &array[largest_index]);

        /* Recursively heapify the affected sub-tree */
        heapify(array, size, largest_index);
    }
}

/*  A function to do counting sort of array according to
    the digit represented by exponent.  */
void radix_countSort ( long* array, long size, long exp)
{
    long* bucketArray, auxiliaryArray[10] = {0};
    register long i;

    bucketArray = (long*)malloc(size * sizeof(long));

    /* Store count of occurrences in auxiliary array */
    for ( i = 0; i < size; ++i )
        ++auxiliaryArray[(array[i]/exp)%10];

    /*  Change auxiliaryArray at the index i so that auxiliaryArray[i] now contains actual
        position of this digit in bucket array.  */
    for ( i = 1; i < 10; ++i )
        auxiliaryArray[i] += auxiliaryArray[i - 1];


    /*  Build the bucket array.  */
    for ( i = size - 1; i >= 0; --i ) {
        bucketArray[auxiliaryArray[ (array[i]/exp)%10 ] - 1] = array[i];

        --auxiliaryArray[(array[i]/exp)%10];
    }

    /*  Copy the bucket array to the original array, so that original array now
        contains sorted numbers according to current digit. */
    for ( i = 0; i < size; ++i )
        array[i] = bucketArray[i];

    free(bucketArray);
}

/** In the library .h. */
void selectionSort ( long* array, long size ) /** I think it don't need an explanation. */
{
    register long i, j;

    for ( i = 0; i < size - 1; ++i )
        for ( j = i + 1; j < size; ++j )
            if(array[i] > array[j])
                SWAP(&array[i], &array[j]);
}

void bubbleSort ( long* array, long size ) /** The same of the previously function. */
{
    register long i, j;

    for ( i = 0; i < size - 1; ++i ) {
        for ( j = 0; j < size - i - 1; ++j )
            if( array[j] > array[j+1] )
                SWAP(&array[j], &array[j+1]);
    }
}

void recursiveBubbleSort ( long* array, long size ) /** If you got it how the bubble sort works, this is the same, but with recursion */
{
    register long i;

    if ( size <= 1 )
        return ;

    for ( i = 0; i < size - 1; ++i )
        if ( array[i] > array[i+1] )
            SWAP(&array[i], &array[i+1]);

    recursiveBubbleSort( array, size - 1 );
}

void insertionSort ( long* array, long size )
{
    long temp; /** Is necessary create a temporary variable. */
    register long i, j;

    /* Sort the array. */
    for ( i = 1; i < size; ++i ) {
        temp = array[i]; /** initializes the temp variable in the index that the array is pointing to. */

        /* Compare what's in the array with what's ahead, if what's in the array is larger, then copy forward. */
        for ( j = i - 1; (j >= 0) && (array[j] > temp); --j )
            array[j+1] = array[j];

        array[j+1] = temp;  /** After that copy the temporary variable into the next index of the array. */
    }
}

void recursiveInsertionSort ( long* array, long size ) /** The same of above function, but with recursion. */
{
    long temp;
    register long i = size - 1, j;

    if ( size <= 1)
        return ;

    recursiveInsertionSort(array, size - 1);

    temp = array[i];
    for ( j = i - 1; (j >= 0) && (array[j] > temp); --j )
        array[j+1] = array[j];

    array[j+1] = temp;
}

void mergeSort ( long* array, long left_index, long right_index )
{
    long middle_index =  left_index + (right_index - left_index) / 2; /** Calculates the middle index, avoiding overflow by larger l. */

    if ( left_index < right_index )
    {
        /* Sort first and second halves */
        mergeSort(array, left_index, middle_index);
        mergeSort(array, middle_index+1, right_index);

        mergeArray(array, left_index, middle_index, right_index);
    }
}

void quickSort ( long* array, long low_index, long high_index ) /** I think its too simple, therefore there's no need of a explanation. */
{
    long partition_index;

    if ( low_index < high_index )
    {
        partition_index = partitionArray(array, low_index, high_index);

        quickSort(array, low_index, partition_index - 1);
        quickSort(array, partition_index + 1, high_index);
    }
}

void heapSort ( long* array, long size )
{
    register long i;

    /* Build heap (rearrange array) */
    for ( i = size / 2 - 1; i >= 0; --i )
        heapify(array, size, i);

    /* One by one extract an element from heap */
    for ( i = size - 1; i >= 0; --i ) {
        /* Move current root to end */
        SWAP(&array[0], &array[i]);

        /* call max heapify on the reduced heap */
        heapify(array, i, 0);
    }
}

void shellSort ( long* array, long size )
{
    register long gap, i, j;
    long temp;

    /* Start with a big gap, then reduce the gap */
    for ( gap = size / 2; gap > 0; gap /= 2 )
    /*  Do a gapped insertion sort for this gap size.
        The first gap elements array[0 ... gap-1] are
        already in gapped order keep adding one more
        element until the entire array is gap sorted  */
        for ( i = gap; i < size; ++i ) {
        /*  add array[index] to the elements that have been gap sorted
            save array[index] in temp and make a hole at position i  */
            temp = array[i];

            /*  shift earlier gap-sorted elements up until the correct
                location for array[index] is found   */
            for ( j = i; (j >= gap) && (array[j - gap] > temp); j -= gap )
                array[j] = array[j - gap];

            /* put temp (the original array[index]) in its correct location */
            array[j] = temp;
    }
}

void countingSort ( long* array, long size )
{
    long* bucketArray,* auxiliaryArray, max;
    register long i;

    bucketArray = (long*)malloc(size * sizeof(long));

    max = array[0];
    for ( i = 1; i < size; ++i ) /** Determines the maximum value */
        if ( array[i] > max )
            max = array[i];

    /* Creates and initializes the auxiliary array */
    auxiliaryArray = (long*)calloc(max + 1, sizeof(long));

    /* Counts the numbers of each number */
    for ( i = 0; i < size; ++i )
        ++auxiliaryArray[array[i]];

    /* Setting up the bucket array */
    for ( i = 1; i <= max; ++i )
        auxiliaryArray[i] += auxiliaryArray[i - 1];

  	for ( i = size - 1; i >= 0; --i ) {
        bucketArray[ auxiliaryArray[ array[i] ] - 1] = array[i];

        --auxiliaryArray[ array[i] ];
    }

    /* Copy the bucket array to the original array. */
    for ( i = 0; i < size; ++i )
        array[i] = bucketArray[i];

    free(auxiliaryArray);
    free(bucketArray);
}

void radixSort ( long* array, long size )
{
    long maxNum, exp;

    /*  Find the maximum number, to know the number of digits.   */
    for ( maxNum = array[0], exp = 1; exp < size; ++exp )
        if ( array[exp] > maxNum )
            maxNum = array[exp];

    /*  Do counting sort for every digit. Note that instead
        of passing digit number, exp is passed. exp is 10^i
        where i is current digit number. */
    for ( exp = 1; maxNum / exp > 0; exp *= 10 )
        radix_countSort(array, size, exp);
}
