#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <stdlib.h>

//Implementation
typedef struct
{
    int policy;
    int age;
} customer;

int compare(const void *l,const void *r)
{
  customer* lv = (customer*)l;
  customer* rv = (customer*)r;
  int value = lv->policy - rv->policy;
  if (value == 0)
    value = lv->age - rv->age;
  return value;
}

void sort(customer *c,int l) {
  qsort(c,l,sizeof(customer),compare);
}

// Test helper
int is_sorted(customer *c,int l) {
  customer pc;
  for(int i=0;i<l;i++) {
    if (i>0) {
      if (c[i].policy < pc.policy) {
        return 0;
      }else{
        if (c[i].policy == pc.policy) {
          if (c[i].age < pc.age) {
            return 0;
          }
        }
      }
    }
    pc = c[i];
  }
  return 1;
}

//Tests
void test_empty() {
  customer c[0];
  sort(c,0);
  assert(is_sorted(c,0));
}

void test_pre_sorted() {
  customer c[3];
  c[0].policy = 1;
  c[0].age = 0;
  c[1].policy = 2;
  c[1].age = 0;
  c[2].policy = 3;
  c[2].age = 0;
  sort(c,3);
  assert(is_sorted(c,3));
}

void test_reverse_sorted() {
  customer c[3];
  c[0].policy = 3;
  c[0].age = 0;
  c[1].policy = 2;
  c[1].age = 0;
  c[2].policy = 1;
  c[2].age = 0;
  sort(c,3);
  assert(is_sorted(c,3));
}

void test_age_sorted() {
  customer c[3];
  c[0].policy = 0;
  c[0].age = 3;
  c[1].policy = 0;
  c[1].age = 1;
  c[2].policy = 0;
  c[2].age = 2;
  sort(c,3);
  assert(is_sorted(c,3));
}

void test_long_sorted() {
  customer c[6];
  c[0].policy = 8;
  c[0].age = 3;
  c[1].policy = 3;
  c[1].age = 2;
  c[2].policy = 1;
  c[2].age = 1;
  c[3].policy = 1;
  c[3].age = 7;
  c[4].policy = 2;
  c[4].age = 4;
  c[5].policy = 1;
  c[5].age = 9;
  sort(c,6);
  assert(is_sorted(c,6));


  int policy[6];
  int age[6];
  for(int i=0;i<6;i++) {
    policy[i] = c[i].policy;
    age[i] = c[i].age;
  }

  for(int i=0;i<6;i++)
    printf("%i(%i),",policy[i],age[i]);
}

int main(int argc, char const *argv[]) {
  test_empty();
  test_pre_sorted();
  test_age_sorted();
  test_reverse_sorted();
  test_long_sorted();
  printf("Test succes");
}
