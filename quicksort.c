// NOTE: Only used to sort vertex costs in nealder mead method now,
// Previously required if sorting concatenated basis positions device side,
// now expecting this to be passed from host.

#ifndef EXTERNAL_CONCAT
#include "shared_defines.h"
#endif

void swapf(FPN *arr, int i, int j) {
  FPN temp = arr[i];
  arr[i] = arr[j];
  arr[j] = temp;
}

void swapi(int *arr, int i, int j) {
  int temp = arr[i];
  arr[i] = arr[j];
  arr[j] = temp;
}

int partition(FPN *arr, int *ids, int l, int h) {
  // index of smaller element
  int i = l - 1;

  FPN pivot = arr[h];

  for (int j = l; j < h; j++) {
    // If current element is smaller than or equal to pivot
    if (arr[j] <= pivot) {
      // increment index of smaller element
      i += 1;
      swapf(arr, i, j);
      swapi(ids, i, j);
    }
  }

  swapf(arr, i + 1, h);
  swapi(ids, i + 1, h);

  return i + 1;
}

void quickSortIterative(const FPN *_arr, int *ids, int *stack, FPN *arr, int Na)
// Source: originally from
// https://github.com/numba/numba/issues/4283#issuecomment-908511793
{

  for (int i = 0; i < Na; i++) {
    arr[i] = _arr[i];
  }

  // low and high indices.
  int l = 0;
  int h = Na - 1;
  int p = 0;

  // stack = [0]  *size
  // ids = list(range(len(arr)))

  // initialize top of stack
  int top = -1;

  // fill ids with range(len(arr))
  for (int i = 0; i < Na; i++) {
    ids[i] = i;
    stack[i] = 0;
  }

  // push initial values of l and h to stack
  top = top + 1;
  stack[top] = l;
  top = top + 1;
  stack[top] = h;

  // Keep popping from stack while is not empty
  while (top >= 0) {

    // Pop h and l
    h = stack[top];
    top = top - 1;
    l = stack[top];
    top = top - 1;

    // Set pivot element at its correct position in
    // sorted array
    p = partition(arr, ids, l, h);

    // If there are elements on left side of pivot,
    // then push left side to stack
    if (p - 1 > l) {
      top = top + 1;
      stack[top] = l;
      top = top + 1;
      stack[top] = p - 1;
    }

    // If there are elements on right side of pivot,
    // then push right side to stack
    if (p + 1 < h) {
      top = top + 1;
      stack[top] = p + 1;
      top = top + 1;
      stack[top] = h;
    }
  }
}
