#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <stdlib.h>

//Implementation

int compare (const void *a,const void *b)
{
  return ( *(int*)a - *(int*)b);
}

void sort(int *a,int l) {
  qsort(a,l,sizeof(int),compare);
}

int is_sorted(int *a,int l) {
  int k=INT_MIN;
  for(int i=0;i<l;i++) {
    if (a[i] < k) 
      return 0;
    else
      k = a[i];
    printf("%i, %i\n",a[i],k);
  }
  return 1;
}

//Tests
void test_empty() {
  int a[] = {};
  sort(a,0);
  assert(is_sorted(a,0));
}

void test_pre_sorted() {
  int a[] = {1,2,3};
  sort(a,3);
  assert(is_sorted(a,3));
}

void test_reverse_sorted() {
  int a[] = {3,2,1};
  sort(a,3);
  assert(is_sorted(a,3));
}

int main(int argc, char const *argv[]) {
  test_empty();
  test_pre_sorted();
  test_reverse_sorted();
  printf("Test succes");
}
