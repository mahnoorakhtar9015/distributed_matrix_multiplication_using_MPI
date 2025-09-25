#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// FUNCTION TO WRITE MATRIX A INTO FILE
void WriteMatrixA(int MatrixSizeA)
{
    FILE *fptr = fopen("./MatrixA.txt" ,"w");

    if(!fptr)
    {
        printf("Error Opening file\n");
        exit(0);
    }

    int Rows = MatrixSizeA;
    int Columns = MatrixSizeA;

    for(int i = 0; i < Rows; i++)
    {
        for(int j = 0; j < Columns; j++)
        {
            fprintf(fptr, "%d,%d,%d,%d",0, i, j, rand() % 100);
            fprintf(fptr, "\n");
        }
    }

    fclose(fptr);

}

// FUNCTION TO WRITE MATRIX B INTO FILE
void WriteMatrixB(int MatrixSizeB)
{

    FILE *fptr = fopen("./MatrixB.txt" ,"w");

    if(!fptr)
    {
        printf("Error Opening file\n");
        exit(0);
    }

    int Rows = MatrixSizeB;
    int Columns = MatrixSizeB;

    for(int i = 0; i < Rows; i++)
    {
        for(int j = 0; j < Columns; j++)
        {
            fprintf(fptr, "%d,%d,%d,%d",1, i, j, rand() % 100);
            fprintf(fptr, "\n");
        }
    }

    fclose(fptr);

}



int main()
{
    // GENERATING RANDOM SIZE OF MATRIX A AND B BETWEEN 2^4 AND 2^8
    srand(time(NULL));  

    int Power = rand() % 5 + 4;
    int MatrixSizeA = 2, MatrixSizeB = 2;  

    for(int i = 0; i < Power - 1; i++)
    {
        MatrixSizeA *=2;
        MatrixSizeB *=2;
    }


    // WRITING MATRIX DATA INTO A FILE
    WriteMatrixA(MatrixSizeA);
    WriteMatrixB(MatrixSizeB);

}