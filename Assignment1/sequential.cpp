#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

const int N = 20;
const int MAX_NUM = 1000;

int main()
{
    clock_t start, end;
    int i, n = N;
    int num[N];
    bool tag = true;

	start = clock();
	
    srand(clock());
    for (i = 0; i < n; i++)
        num[i] = rand() % MAX_NUM;

    printf("Name: Yuncong Cui\nStudent ID: 118010045\nHomework 1, Odd Even Sort, Sequential Implementation\n\n---------------Before Odd Even Sort---------------\nArray: ");
    for (i = 0; i < n; i++)
        printf("%d ", num[i]);

    while (tag) {
        tag = false;

        for (i = 1; i < n - 1; i += 2) {
            if (num[i] > num[i + 1]) {
                swap(num[i], num[i + 1]);
                tag = true;
            }
        }
        for (i = 0; i < n - 1; i += 2) {
            if (num[i] > num[i + 1]) {
                swap(num[i], num[i + 1]);
                tag = true;
            }
        }
    }

    end = clock();

    printf("\n\n---------------After Odd Even Sort----------------\nArray: ");
    for (i = 0; i < n; i++)
        printf("%d ", num[i]);
    printf("\n\n------------------Time Analysis-------------------\nThe number of elements: %d\nRuntime: %lf seconds\n", n, (double)(end - start) / CLOCKS_PER_SEC);

    return(0);
}
