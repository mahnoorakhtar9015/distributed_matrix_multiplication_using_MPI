4#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stddef.h>

typedef struct key_value
{
	int i;
	int k;
	int matrix;
	int j;
	int value_at_index;
};

typedef struct reducer_key_value
{
	int i;
	int k;
	int value_at_index;
};

struct reducer_key_value * reducer(struct key_value*array, int blocksize,int row_size)
{
	struct reducer_key_value * RKV=(struct reducer_key_value*)(malloc(sizeof(struct reducer_key_value)*(blocksize/(row_size+row_size))));
	int count=0;
	int mul;
	int value=0;
	int index_of_reducer=0;

	for (int i = 0 ; i<blocksize ;i++)     
	{ 
		if(count<(row_size+row_size))
		{   
				
			if(i%2==0)
			{
				mul=(array[i].value_at_index)*(array[i+1].value_at_index);

				value+=mul;
				mul=0;

				count+=2;   
			}
						
			if(count==(row_size+row_size))
			{
				count=0;
				RKV[index_of_reducer].i=array[i].i;
				RKV[index_of_reducer].k=array[i].k;
				RKV[index_of_reducer].value_at_index=value;

				value=0;
				index_of_reducer+=1;
			}
		}
	}

	return RKV;
}

struct key_value * mapper(int key,int array[4])
{
	struct key_value * KV = (struct key_value *)(malloc(sizeof(struct key_value)*key));
	for(int i=0 ; i<key ; i++)
	{
		if(array[0]==0)
		{
			KV[i].i=array[1];
			KV[i].k=i;
			KV[i].matrix=array[0];
			KV[i].j=array[2];
			KV[i].value_at_index=array[3];
		}
		else
		{
			KV[i].i=i;
			KV[i].k=array[2];
			KV[i].matrix=array[0];
			KV[i].j=array[1];
			KV[i].value_at_index=array[3];
		}
	}
   	return KV;
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int size, rank;

	//const int root_rank = 0;
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int * sendCount = NULL;
	int blocksize;
	int * myArray = NULL;
	int * matrix = NULL;

	struct key_value * final_keys = NULL;
	struct key_value * KV_pair = NULL;

	struct key_value * input_keys_red = NULL;


    FILE *fptr;
    fptr = fopen("./MatrixA.txt", "r");

    if (fptr == NULL) {
        printf("Error: Could not open file.\n");
        return 1;
    }

    int NumberLines = 0;
    char c;
    while ((c = fgetc(fptr)) != EOF) 
    { 
        if (c == '\n') 
        {
            NumberLines++;
        }
    }
    fclose(fptr);

	const int MATRIX_SIZE = NumberLines;					// The matrix size will be based on the file.
	const int row_size = (int)(sqrt(MATRIX_SIZE));			// Square root of MATRIX_SIZE

    int len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(name, &len);

	if (rank == 0) {
        printf("\nMaster with process_id %d running on %s\n", rank, name);
		
		// CREATE 1-D ARRAY OF SIZE MATRIX SIZE(NUMBER OF LINES) * 4 * 2
		matrix = (int *) (malloc(sizeof(int) * MATRIX_SIZE * 4 * 2));

		// READING DATA FROM MATRIX A INTO 1-D ARRAY
		fptr = fopen("./MatrixA.txt", "r");

		if (fptr == NULL) 
		{
			printf("Error: Could not open MatrixA.txt.\n");
			return 1;
		}

		int i = 0, value;
		while (fscanf(fptr, "%d,", &value) != EOF) 
		{
			matrix[i++] = value;
		}

		fclose(fptr);


		// READING DATA FROM MATRIX B INTO 1-D ARRAY
		fptr = fopen("./MatrixB.txt", "r");

		if (fptr == NULL) 
		{
			printf("Error: Could not open MatrixB.txt.\n");
			return 1;
		}

		while (fscanf(fptr, "%d,", &value) != EOF) 
		{
			matrix[i++] = value;
		}

		fclose(fptr);

		// for (int i=0; i < MATRIX_SIZE*2; i++) {
		// 	for (int j = 0; j < 4; j++) {
		// 		printf("matrix[%d][%d]: %d\n", i, j, matrix[(i*4)+j]);
		// 	}
		// }

		sendCount = (int *) (malloc(sizeof(int) * size));
        for (int i=0; i < size; i++) {
			sendCount[i] = 0;
        }

		int nproc = size-1;	// as one process will act as node manager
		int mappers = 0;
		for (int i=nproc; i > 0; i--) {
			if ((MATRIX_SIZE*2) % i == 0) {
				mappers = i;
				break;
			}
		}

		for (int i=1; i <= mappers; i++) {
			sendCount[i] = (MATRIX_SIZE*2)/mappers;
		}

        for (int i=0; i<size; i++)
        {
            if (sendCount[i] != 0)
            {
                printf("Task Map assigned to process %d\n", i);
            }
        }
	}

	// Informs every process the bllocksize for which it will be making keys
	MPI_Scatter(sendCount, 1, MPI_INT, &blocksize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//printf("Process: %d has block size of %d for mapper\n", rank, blocksize);

	MPI_Barrier(MPI_COMM_WORLD);

	// Scattering the matrix to available processes
	switch(rank)
	{
		case 0:
		{
			int * displacement = (int *) (malloc(sizeof(int) * size));
			displacement[0] = 0;

			for (int i=1; i < size; i++) {
				sendCount[i] *= 4;
				displacement[i] = displacement[i-1] + sendCount[i-1];
			}

			myArray = (int *) (malloc(sizeof(int) * 4 * blocksize));
			MPI_Scatterv(matrix, sendCount, displacement, MPI_INT, myArray, blocksize*4, MPI_INT, 0, MPI_COMM_WORLD);
			break;
		}
		default:
		{
			myArray = (int *) (malloc(sizeof(int) * 4 * blocksize));
			MPI_Scatterv(NULL, NULL, NULL, MPI_INT, myArray, blocksize*4, MPI_INT, 0, MPI_COMM_WORLD);
			
            if (blocksize != 0)
            {
                printf("Process %d received task Map on %s\n", rank, name);
            }
            
            break;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// for (int x=0; x < blocksize; x++) {
	// 	for (int y=0; y < 4; y++) {
	// 		printf("Process: %d val: %d\n", rank, myArray[(x*4)+y]);
	// 	}
	// 	printf("\n");
	// }

	// Mapping Operation
	KV_pair = (struct key_value *)(malloc(sizeof(struct key_value) * row_size * blocksize));
    struct key_value * temp = NULL;

	for(int i=0 ; i<blocksize ; i++)
	{
		int arr[4]= {myArray[(i*4)+0],myArray[(i*4)+1],myArray[(i*4)+2],myArray[(i*4)+3]};
		temp = mapper(row_size,arr);
		for(int j=0;j<row_size;j++) {
			KV_pair[(i*row_size) + j].i = temp[j].i;
			KV_pair[(i*row_size) + j].j = temp[j].j;
			KV_pair[(i*row_size) + j].k = temp[j].k;
			KV_pair[(i*row_size) + j].matrix = temp[j].matrix;
			KV_pair[(i*row_size) + j].value_at_index = temp[j].value_at_index;
			//printf("P:%d, i=>%d, k=>%d, m=>%d, j=>%d, v=>%d\n", rank, temp[j].i, temp[j].k, temp[j].matrix, temp[j].j, temp[j].value_at_index);
		}

		free(temp);
		temp = NULL;
	}

	// de-allocating memory from myArray
	free(myArray);
	myArray = NULL;
	MPI_Barrier(MPI_COMM_WORLD);

	int bl[5] = {1, 1, 1, 1, 1};
	MPI_Datatype types[5] = {MPI_INT, MPI_INT,MPI_INT,MPI_INT,MPI_INT};
	MPI_Aint offsets[5] = {
		offsetof(struct key_value, i),
		offsetof(struct key_value, k),
		offsetof(struct key_value, matrix),
		offsetof(struct key_value, j),
		offsetof(struct key_value, value_at_index)
	};
	MPI_Datatype myStruct;
	MPI_Type_create_struct(5, bl, offsets,types, &myStruct);
	MPI_Type_commit(&myStruct);

	if (rank == 0) {
		final_keys = (struct key_value *)(malloc(sizeof(struct key_value) * row_size * row_size * row_size * 2));
		
		int * displacement = (int *) (malloc(sizeof(int) * size));

		sendCount[0] = 0;
		displacement[0] = 0;

		for (int i=1; i < size; i++) {
			sendCount[i] /= 4;
			sendCount[i] *= row_size;	// as each value has now mapped based on row_size
			displacement[i] = displacement[i-1] + sendCount[i-1];
		}

		MPI_Gatherv(KV_pair, row_size * blocksize, myStruct, final_keys, sendCount, displacement, myStruct, 0, MPI_COMM_WORLD);
		
		free(KV_pair);
		KV_pair = NULL;
		MPI_Barrier(MPI_COMM_WORLD);

		// for(int j=0;j< row_size * row_size * row_size * 2 ;j++) {
		// 	printf("P:%d, i=>%d, k=>%d, m=>%d, j=>%d, v=>%d\n", rank, final_keys[j].i, final_keys[j].k, final_keys[j].matrix, final_keys[j].j, final_keys[j].value_at_index);
		// }

		// At this point the root process has an array containing all the key-value
		// from the mapper named final_keys.
		// -> Shuffler will be implemented here.
		// -> Then different key-value blocks will be send to each process
		//    using MPI_Scatterv
	}
	else
	{
		MPI_Gatherv(KV_pair, row_size * blocksize, myStruct, NULL, NULL, NULL, myStruct, 0, MPI_COMM_WORLD);
		
        if (blocksize != 0)
        {
            printf("Process %d has completed task Map\n", rank);
        }

		free(KV_pair);
		KV_pair = NULL;
		MPI_Barrier(MPI_COMM_WORLD);

		// All the other processes except root will get the blocks of key-value
		// for reduction here using MPI_Scatterv
	}

	// SHUFFLER IMPLEMENTATION

	if (rank == 0)
	{
		// re-arranging of the final keys from the mapper
		struct key_value * temp_keys = (struct key_value *)(malloc(sizeof(struct key_value) * row_size * row_size * row_size * 2));

		int index = 0;
		for (int i=0; i < (row_size * row_size * row_size * 2); i++)
		{
			index = (final_keys[i].i * row_size * row_size * 2) +
					(final_keys[i].k * row_size * 2) +
					(final_keys[i].j * 2) +
					final_keys[i].matrix;

			temp_keys[index].i = final_keys[i].i;
			temp_keys[index].j = final_keys[i].j;
			temp_keys[index].k = final_keys[i].k;
			temp_keys[index].matrix = final_keys[i].matrix;
			temp_keys[index].value_at_index = final_keys[i].value_at_index;
		}

		for (int i=0; i < (row_size * row_size * row_size * 2); i++)
		{
			final_keys[i].i = temp_keys[i].i;
			final_keys[i].j = temp_keys[i].j;
			final_keys[i].k = temp_keys[i].k;
			final_keys[i].matrix = temp_keys[i].matrix;
			final_keys[i].value_at_index = temp_keys[i].value_at_index;
		}

		free(temp_keys);
		temp_keys = NULL;

		// compute send count array and number of processes required

		free(sendCount);
		sendCount = (int *) (malloc(sizeof(int) * size));
        for (int i=0; i < size; i++) {
			sendCount[i] = 0;
        }

		int nproc = size-1;	// as one process will act as node manager
		int reducers = 0;
		for (int i=nproc; i > 0; i--) {
			if ((row_size * row_size) % i == 0) {
				reducers = i;
				break;
			}
		}

		for (int i=1; i <= reducers; i++) {
			sendCount[i] = (row_size * row_size * row_size * 2)/reducers;
		}

        for (int i=0; i<size; i++)
        {
            if (sendCount[i] != 0)
            {
                printf("Task Reduce assigned to process %d\n", i);
            }
        }
	}

	MPI_Scatter(sendCount, 1, MPI_INT, &blocksize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//printf("Process: %d has block size of %d for reducer\n", rank, blocksize);

	MPI_Barrier(MPI_COMM_WORLD);

	switch(rank)
	{
		case 0:
		{
			int * displacement = (int *) (malloc(sizeof(int) * size));
			displacement[0] = 0;

			for (int i=1; i < size; i++) {
				displacement[i] = displacement[i-1] + sendCount[i-1];
			}

			input_keys_red = (struct key_value *)(malloc(sizeof(struct key_value) * blocksize));
			MPI_Scatterv(final_keys, sendCount, displacement, myStruct, input_keys_red, blocksize, myStruct, 0, MPI_COMM_WORLD);
			break;
		}
		default:
		{
			input_keys_red = (struct key_value *)(malloc(sizeof(struct key_value) * blocksize));
			MPI_Scatterv(NULL, NULL, NULL, myStruct, input_keys_red, blocksize, myStruct, 0, MPI_COMM_WORLD);
			
            if (blocksize != 0)
            {
                printf("Process %d received task Reduce on %s\n", rank, name);
            }
            
            break;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		free(final_keys);
		final_keys = NULL;
	}

	// for (int j=0; j < blocksize; j++) {
	// 	printf("P:%d, i=>%d, k=>%d, m=>%d, j=>%d, v=>%d\n", rank, input_keys_red[j].i, input_keys_red[j].k, input_keys_red[j].matrix, input_keys_red[j].j, input_keys_red[j].value_at_index);
	// }

	// Reduction Code will be written here


	// Mapping Operation

    struct reducer_key_value * Reducer_KV = NULL;
	struct reducer_key_value * reducer_final_keys=NULL;	

	Reducer_KV = reducer(input_keys_red,blocksize,row_size);

	MPI_Barrier(MPI_COMM_WORLD);

	free(input_keys_red);
	input_keys_red = NULL;
	
	
	//gathering of data 
	
	int blr[3] = {1, 1, 1};
	MPI_Datatype types_R[3] = {MPI_INT, MPI_INT,MPI_INT};
	MPI_Aint offsets_R[3] = {
		offsetof(struct reducer_key_value, i),
		offsetof(struct reducer_key_value, k),
		offsetof(struct reducer_key_value, value_at_index),
	};
	MPI_Datatype myStruct_R;
	MPI_Type_create_struct(3, blr, offsets_R,types_R, &myStruct_R);
	MPI_Type_commit(&myStruct_R);
        
	if (rank == 0) {
		reducer_final_keys= (struct reducer_key_value *)(malloc(sizeof(struct reducer_key_value) * (row_size * row_size)));
		
		int * displacement = (int *) (malloc(sizeof(int) * size));

		sendCount[0] = 0;
		displacement[0] = 0;

		for (int i=1; i < size; i++) {
			
			sendCount[i] /= (row_size+row_size);	// as each value has now mapped based on row_size
			displacement[i] = displacement[i-1] + sendCount[i-1];
		}

		MPI_Gatherv(Reducer_KV, blocksize/(row_size+row_size), myStruct_R, reducer_final_keys, sendCount, displacement, myStruct_R, 0, MPI_COMM_WORLD);

		free(Reducer_KV);
		Reducer_KV = NULL;
		MPI_Barrier(MPI_COMM_WORLD);
		
		// for(int j=0;j< row_size * row_size  ;j++) {
		// 	printf("P:%d, i=>%d, k=>%d,v=>%d\n", rank, reducer_final_keys[j].i, reducer_final_keys[j].k, reducer_final_keys[j].value_at_index);
		// }

	}
	else
	{
		MPI_Gatherv(Reducer_KV, blocksize/(row_size+row_size), myStruct_R, NULL, NULL, NULL, myStruct_R, 0, MPI_COMM_WORLD);
		
        if (blocksize != 0)
        {
            printf("Process %d has completed task Reduce\n", rank);
        }

		free(Reducer_KV);
		Reducer_KV = NULL;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// Now the root process has key-value pairs after the reduction process.
	// Saved in 'reducer_final_keys' variable
    

	if (rank == 0)
	{
        printf("Job has been completed!\n");
		// re-arranging of the final keys from the reducer
		struct reducer_key_value * temp_keys_r = (struct reducer_key_value *)(malloc(sizeof(struct reducer_key_value) * row_size * row_size));

		int index = 0;
		for (int i=0; i < (row_size * row_size); i++)
		{
			index = (reducer_final_keys[i].i * row_size) + reducer_final_keys[i].k;

			temp_keys_r[index].i = reducer_final_keys[i].i;
			temp_keys_r[index].k = reducer_final_keys[i].k;
			temp_keys_r[index].value_at_index = reducer_final_keys[i].value_at_index;
		}

		for (int i=0; i < (row_size * row_size); i++)
		{
			reducer_final_keys[i].i = temp_keys_r[i].i;
			reducer_final_keys[i].k = temp_keys_r[i].k;
			reducer_final_keys[i].value_at_index = temp_keys_r[i].value_at_index;
		}

		free(temp_keys_r);
		temp_keys_r = NULL;

		// for(int j=0;j< row_size * row_size  ;j++) {
		// 	printf("P:%d, i=>%d, k=>%d, v=>%d\n", rank, reducer_final_keys[j].i, reducer_final_keys[j].k, reducer_final_keys[j].value_at_index);
		// }

		int ** computed_matrix = (int**)(malloc(sizeof(int*) * row_size));
		int ** actual_matrix = (int**)(malloc(sizeof(int*) * row_size));
		for (int i=0; i < row_size; i++)
		{
			computed_matrix[i] = (int*)(malloc(sizeof(int) * row_size));
			actual_matrix[i] = (int*)(malloc(sizeof(int) * row_size));
		}

		// computed_matrix
		for (int i=0; i<row_size; i++) {
			for (int j=0; j<row_size; j++) {
				computed_matrix[i][j] = reducer_final_keys[(i*row_size)+j].value_at_index;
			}
		}

		free(reducer_final_keys);
		reducer_final_keys = NULL;
		
		// actual_matrix
		for (int i=0; i<row_size; i++) {
			for (int j=0; j<row_size; j++) {
				actual_matrix[i][j] = 0;
				for (int k=0; k<row_size; k++) {
					actual_matrix[i][j] += (matrix[(i*row_size*4)+(k*4)+3] * matrix[(MATRIX_SIZE*4)+(k*row_size*4)+(j*4)+3]);
				}
			}
		}

		free(matrix);
		matrix = NULL;

		// // display
		// for (int i=0; i<row_size; i++) {
		// 	for (int j=0; j<row_size; j++) {
		// 		printf("%d, %d\t\t", computed_matrix[i][j], actual_matrix[i][j]);
		// 	}
		// 	printf("\n");
		// }



		// WRITING RESULTANT MATRIX C INTO A FILE
		fptr = fopen("./MatrixC.txt" ,"w");
		if(!fptr)
		{
			printf("Error Opening file\n");
			exit(0);
		}

		// int row_size = 4;

		for(int i = 0; i < row_size ;i++)
		{
			for(int j = 0; j < row_size; j++)
			{
				fprintf(fptr, "%d", computed_matrix[i][j]);
				if (j != row_size - 1) {
					fprintf(fptr, ",");
				}
			}
			fprintf(fptr, "\n");
		}
		fclose(fptr);


        // comparison
        int equal = 1;
		for (int i=0; i<row_size; i++) {
			for (int j=0; j<row_size; j++) {
			    if (computed_matrix[i][j] != actual_matrix[i][j])
                {
                    equal = 0;
                    break;
                }
			}
            if (equal == 0)
            {
                break;
            }
		}
        
        if (equal == 1)
        {
            printf("Matrix comparison function returned: True\n");
        }
        else
        {
            printf("Matrix comparison function returned: False\n");
        }

		for (int i=0; i < row_size; i++)
		{
			free(computed_matrix[i]);
			free(actual_matrix[i]);
		}
		free(computed_matrix);
		free(actual_matrix);

		free(sendCount);

	}

	MPI_Finalize();

    return EXIT_SUCCESS;
}
