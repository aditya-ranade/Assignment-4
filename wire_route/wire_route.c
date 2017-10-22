#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <assert.h>
#include <stdbool.h>
#include "mpi.h"
#include "wire_route.h"
#include <limits.h>




        // int block_lengths[6] = {1, 1, 1, 1, 1, 1};
        // MPI_Datatype types[6] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
        // MPI_Aint     offsets[6];

        // offsets[0] = offsetof(wire_t, x1);
        // offsets[1] = offsetof(wire_t, x2);
        // offsets[2] = offsetof(wire_t, y1);
        // offsets[3] = offsetof(wire_t, y2);
        // offsets[4] = offsetof(wire_t, bend);
        // offsets[5] = offsetof(wire_t, path);

        // MPI_Type_create_struct(6, block_lengths, offsets, types, &mpi_wire_type);
        // MPI_Type_commit(&mpi_wire_type);

        // int block_lengths_c[4] = {1, 1, costs->rows*costs->cols, costs->rows * costs->cols};
        // MPI_Datatype types_c[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
        // MPI_Aint     offsets_c[4];

        // offsets_c[0] = offsetof(costs_t, array);
        // offsets_c[1] = offsetof(costs_t, transpose);
        // offsets_c[2] = offsetof(costs_t, rows);
        // offsets_c[3] = offsetof(costs_t, cols);

        // MPI_Type_create_struct(4, block_lengths_c, offsets_c, types_c, &mpi_cost_type);
        // MPI_Type_commit(&mpi_cost_type);
    //}
    //MPI_Bcast(wires, num_wires, mpi_wire_type, root, MPI_COMM_WORLD);
    //MPI_Bcast(costs, 1, mpi_cost_type, root, MPI_COMM_WORLD);

typedef struct {
    int x1;
    int x2;
    int y1;
    int y2;
    int bend;
    int path;
} wire_t;

typedef struct {
    int *array;
    int *transpose;
    int rows;
    int cols;
} costs_t;

wire_t *wires;
costs_t *costs;

int num_wires;

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))




static inline void init(int numRows, int numCols, int numWires)
{
    costs = (costs_t *)malloc(sizeof(costs_t));
    wires = (wire_t *)malloc(sizeof(wire_t)*numWires);
    costs->array = (int *)calloc(numRows*numCols, sizeof(int));
    costs->transpose = (int *)calloc(numRows*numCols, sizeof(int));
    costs->rows = numRows;
    costs->cols = numCols;
    num_wires = numWires;
        
}

// Initialize a given wire
static inline void initWire(int wireIndex, int x1, int y1, int x2, int y2)
{
    (&wires[wireIndex])->x1 = x1;
    (&wires[wireIndex])->y1 = y1;
    (&wires[wireIndex])->x2 = x2;
    (&wires[wireIndex])->y2 = y2;
}

// Return number of rows
static inline int getNumRows()
{
    return costs->rows;
}

// Return number of cols
static inline int getNumCols()
{
    return costs->cols;
}

// Return number of wires
static inline int getNumWires()
{
    return num_wires;
}

// Get cost array entry
static inline int getCost(int row, int col)
{
    int dim_x = costs->cols;
    return costs->array[row*dim_x + col];
}

// Get a wire placement. Returns number of points (should be 2-4 for 0-2 bends).
static inline int getWire(int wireIndex, int* x1, int* y1, int* x2, int* y2, int* x3, int* y3, int* x4, int* y4)
{
    int wx1 = wires[wireIndex].x1;
    int wy1 = wires[wireIndex].y1;
    int wx2 = wires[wireIndex].x2;
    int wy2 = wires[wireIndex].y2;
    int path = wires[wireIndex].path;
    int bend = wires[wireIndex].bend;


    if (wx1 == wx2 || wy1 == wy2) {
        *x1 = wx1;
        *y1 = wy1;
        *x2 = wx2;
        *y2 = wy2;
        return 2;
    }

    if (path == 0) {
        if (bend == wx2) {
            *x1 = wx1;
            *y1 = wy1;
            *x2 = wx2;
            *y2 = wy1;
            *x3 = wx2;
            *y3 = wy2;
            return 3;
        }
        else {
            *x1 = wx1;
            *y1 = wy1;
            *x2 = bend;
            *y2 = wy1;
            *x3 = bend;
            *y3 = wy2;
            *x4 = wx2;
            *y4 = wy2;
            return 4;
        }
    }
    else {
        if (bend == wy2) {
            *x1 = wx1;
            *y1 = wy1;
            *x2 = wx1;
            *y2 = wy2;
            *x3 = wx2;
            *y3 = wy2;
            return 3;
        }
        else {
            *x1 = wx1;
            *y1 = wy1;
            *x2 = wx1;
            *y2 = bend;
            *x3 = wx2;
            *y3 = bend;
            *x4 = wx2;
            *y4 = wy2;
            return 4;
        }
    }
}


void increment_cost_x(costs_t *costs, int x1, int x2, int y, int dim_x, int dim_y, int sum) {
    int multiplier = 1;
    if (x1 > x2) multiplier = -1;
    while (x1 != x2) {  
        costs->array[y*dim_x + x1] += sum;
        costs->transpose[x1*dim_y + y] += sum;
        x1 += multiplier;
    }
}


void increment_cost_y(costs_t *costs, int y1, int y2, int x, int dim_x, int dim_y, int sum) {
    int multiplier = 1;
    if (y1 > y2) multiplier = -1;
    while (y1 != y2) {
        costs->array[y1*dim_x + x] += sum;
        costs->transpose[x*dim_y + y1] += sum;
        y1 += multiplier;
    }
}


void traverse_path(int x1, int y1, int x2, int y2, int path, int bend, costs_t *costs, int sum) {
    
    int dim_x = costs->cols;
    int dim_y = costs->rows;
    if (path == 0) {
        increment_cost_x(costs, x1, bend, y1, dim_x, dim_y, sum);
        increment_cost_y(costs, y1, y2, bend, dim_x, dim_y, sum);
        increment_cost_x(costs, bend, x2, y2, dim_x, dim_y, sum);
    }

    else {
        increment_cost_y(costs, y1, bend, x1, dim_x, dim_y, sum);
        increment_cost_x(costs, x1, x2, bend, dim_x, dim_y, sum);
        increment_cost_y(costs, bend, y2, x2, dim_x, dim_y, sum);
    }
    costs->array[y2*dim_x + x2] += sum;
    costs->transpose[x2*dim_y + y2] += sum;
    //fprintf(stdout, "%d \n", costs->array[y2*dim_x + x2])
    
}



void traverse_path1(int x1, int y1, int x2, int y2, int path, int bend, costs_t *costs,
                                    int *min_costs, int* agg_costs, int i) {
    
    int x_multiplier = (x1 > x2) ? -1 : 1;
    int y_multiplier = (y1 > y2) ? -1 : 1;

    int dimx = costs->cols;
    int dimy = costs->rows;
    int *array = costs->array;
    int *transpose = costs->transpose;

    int val1 = 0;
    int val2 = 0;
    int val3 = 0;
    int val4 = 0;
    int aggCost1 = 0;
    int aggCost2 = 0;
    int aggCost3 = 0;
    int aggCost4 = 0;
    int maxCost1 = 0;
    int maxCost2 = 0;
    int maxCost3 = 0;
   
    if (path == 0) {    
        while (x1 != bend) {
            val1 = costs->array[y1*dimx + x1] + 1;
            aggCost1 += val1;
            if (val1 > maxCost1) maxCost1 = val1;       
            x1 += x_multiplier; 
        }
        while (y1 != y2) {
            val2 = costs->transpose[x1*dimy + y1] + 1;
            aggCost2 += val2;
            if (val2 > maxCost2) maxCost2 = val2;
            y1 += y_multiplier;     
        }
        while (x1 != x2) {
            val3 = costs->array[y2*dimx + x1] + 1;
            aggCost3 += val3;
            if (val3 > maxCost3) maxCost3 = val3;
            x1 += x_multiplier;
        }
    }
    else {
        while (y1 != bend) {
            val1 = costs->transpose[x1*dimy + y1] + 1;
            aggCost1 += val1;
            if (val1 > maxCost1) maxCost1 = val1;       
            y1 += y_multiplier; 
        }
        while (x1 != x2) {
            val2 = costs->array[y1*dimx + x1] + 1;
            aggCost2 += val2;
            if (val2 > maxCost2) maxCost2 = val2;
            x1 += x_multiplier;     
        }
        while (y1 != y2) {
            val3 = costs->transpose[x1*dimy + y1] + 1;
            aggCost3 += val3;
            if (val3 > maxCost3) maxCost3 = val3;
            y1 += y_multiplier;
        }
    }
    val4 = array[y2*dimx + x2] + 1;
    aggCost4 = val4;
    int j = MAX(MAX(MAX(aggCost1, aggCost3), aggCost2), aggCost4);
    agg_costs[i] = aggCost1 + aggCost2 + aggCost3 + aggCost4;
}





static inline void set_random_path(wire_t *w) {
    int x_multiplier = (w->x1 > w->x2) ? -1 : 1;
    int y_multiplier = (w->y1 > w->y2) ? -1 : 1;
    //fprintf(stdout, "%d %d ", w->x1, w->x2);
    int num_paths_x = abs(w->x1 - w->x2);
    int num_paths_y = abs(w->y1 - w->y2);

    int number = rand() % (num_paths_x+num_paths_y);
    if (number < num_paths_x) {
        w->path = 0;
        w->bend = w->x1 + x_multiplier*(number + 1);
    }
    else {
        w->path = 1;
        int offset = number - num_paths_x;
        w->bend = w->y1 + y_multiplier*(offset + 1);
    }
}


// Perform computation, including reading/writing output files


// void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations)
// {
    

//     const int root = 0;
//     int tag0 = 0; 
//     int tag1 = 1;
//     MPI_Status status0, status1;
//     int source;
//     int choice;
//     int numPaths;
//     int numPaths_x;
//     int path, bend;
//     int x1, y1, x2, y2;
//     int x_multiplier, y_multiplier;
//     MPI_Datatype mpi_wire_type;
//     MPI_Datatype mpi_cost_type;
//     int mC;

//     readInput(inputFilename);
//     int size = costs->cols * costs->rows;


//     int *min_costs = (int *)malloc(sizeof(int)*costs->cols*costs->rows);
//     int *agg_costs = (int *)malloc(sizeof(int)*costs->cols*costs->rows);


//     int batchSize = num_wires/1000 + 1;
//     for (int b = 0; b < num_wires; b+= batchSize) {
//         int span = (batchSize + nproc - 1)/ nproc;
//         int startIndex = MIN(b + batchSize, b + (procID * span));
//         startIndex = MIN(num_wires, startIndex);
//         int endIndex = MIN(b + batchSize, b + (startIndex + span));
//         endIndex = MIN(num_wires, endIndex);

//         for (int i = startIndex; i < endIndex; i++) {

//             wire_t *w = &wires[i];
//             x1 = w->x1;
//             y1 = w->y1;
//             x2 = w->x2;
//             y2 = w->y2;

//             x_multiplier = (x1 > x2) ? -1 : 1;
//             y_multiplier = (y1 > y2) ? -1 : 1;

//             numPaths = abs(x1 - x2) + abs(y1-y2);
//             numPaths_x = abs(x1 - x2);


//             //traverse_path(w->x1, w->y1, w->x2, w->y2, w->path, w->bend, costs, -1);

//             choice = rand() % 2;

//             if (choice == 0) {
//                 for (int j = 0; j < numPaths; j++) {
//                     if (j < numPaths_x) {
//                         path = 0;
//                         bend = x1 + (x_multiplier * (j + 1));
//                         traverse_path1(x1, y1, x2, y2, path, bend, costs, min_costs, 
//                                                                     agg_costs, j);
//                     }
//                     else {
//                         path = 1;
//                         bend = y1 + (y_multiplier * (j - numPaths_x + 1));
//                         traverse_path1(x1, y1, x2, y2, path, bend, costs, min_costs, 
//                                                                     agg_costs, j);
//                     }
//                 }

//                 int index = 0;
//                 int minCost = INT_MAX;
//                 for (int k = 0; k < numPaths; k++) {
//                     if (min_costs[k] < minCost) {
//                         minCost = min_costs[k];
//                         index = k;
//                     }
//                     if (min_costs[k] == minCost) {
//                         if (agg_costs[k] < agg_costs[index]) index = k;
//                     }
//                 }

                
//                 if (index < numPaths_x) {
//                     w->path = 0;
//                     w->bend = w->x1 + (x_multiplier *(index+ 1));
//                 }
//                 else {
//                     w->path = 1;
//                     w->bend = w->y1 + (y_multiplier * (index - numPaths_x + 1));
//                 }
//             }
//             else {
//                 int number = rand() % (numPaths);
//                 if (number < numPaths_x) {
//                     w->path = 0;
//                     w->bend = w->x1 + x_multiplier*(number + 1);
//                 }
//                 else {
//                     w->path = 1;
//                     int offset = number - numPaths_x;
//                     w->bend = w->y1 + y_multiplier*(offset + 1);
//                 }
//             }

//         }

//         if (procID != root) {
//             MPI_Send(&wires[startIndex], endIndex - startIndex, MPI_INT, root, tag0, MPI_COMM_WORLD);
//         }

//         if (procID == root) {
//             for (int source = 1; source < nproc; source++) {
//                  int startIndex = MIN(b + batchSize, b + (source * span));
//                  int endIndex = MIN(b + batchSize, b + (startIndex + span));
//                  endIndex = MIN(num_wires, endIndex);
//                  startIndex = MIN(num_wires, startIndex);


//                  MPI_Recv(&wires[startIndex], endIndex - startIndex, MPI_INT, source, tag0, MPI_COMM_WORLD, &status0);
//             }
//             int end = b + batchSize;
//             for (int i = b; i < end; i++) {
//                 traverse_path(wires[i].x1, wires[i].y1, wires[i].x2, wires[i].y2, wires[i].path, wires[i].bend, costs, 1);
//             }
//         }
//         MPI_Bcast(costs->array, size, MPI_INT, root, MPI_COMM_WORLD);
//         MPI_Bcast(costs->transpose, size, MPI_INT, root, MPI_COMM_WORLD);
//     }

//     if (procID == root) {
//         int maxCost = 0;
//         for (int i = 0; i < costs->cols*costs->rows; i++) {
//             if (costs->array[i] > maxCost) maxCost = costs->array[i];
//         }
//         fprintf(stdout, "cost = %d \n", maxCost);
//     }

//     free(costs->array);
//     free(costs->transpose);
//     free(costs);
//     free(wires);
//     free(min_costs);
//     free(agg_costs);
// }












void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations)
{
    

    const int root = 0;
    int tag0 = 0; 
    int tag1 = 1;
    MPI_Status status0, status1;
    int source;
    int choice;
    int span, startIndex, endIndex;
    int numPaths;
    int numPaths_x;
    int path, bend;
    int x1, y1, x2, y2;
    int x_multiplier, y_multiplier;
    MPI_Datatype mpi_wire_type;
    MPI_Datatype mpi_cost_type;
    int mC;

    readInput(inputFilename);
    int size = costs->cols * costs->rows;


    int *min_costs = (int *)calloc(costs->cols*costs->rows, sizeof(int));
    int *agg_costs = (int *)calloc(costs->cols*costs->rows, sizeof(int));



    for (int i = 0; i < num_wires; i++) {
        wire_t w = wires[i];
        x1 = w.x1;
        y1 = w.y1;
        x2 = w.x2;
        y2 = w.y2;

        //if (procID != root) fprintf(stdout, "%d %d %d %d \n", x1, y1, x2, y2);

        x_multiplier = (x1 > x2) ? -1 : 1;
        y_multiplier = (y1 > y2) ? -1 : 1;

        numPaths = abs(x1 - x2) + abs(y1-y2);
        numPaths_x = abs(x1 - x2);

        //traverse_path(w->x1, w->y1, w->x2, w->y2, w->path, w->bend, costs, -1);

        if (procID == root) choice = 0;
        MPI_Bcast(&choice, 1, MPI_INT, root, MPI_COMM_WORLD);

        if (choice == 0) {
            span = (numPaths + nproc - 1)/ nproc;
            startIndex = MIN(numPaths, procID * span);
            endIndex = MIN(numPaths, startIndex + span);

            for (int j = startIndex; j < endIndex; j++) {
                if (j < numPaths_x) {
                    path = 0;
                    bend = x1 + (x_multiplier * (j + 1));
                    traverse_path1(x1, y1, x2, y2, path, bend, costs, min_costs, 
                                                                agg_costs, j);
                }
                else {
                    path = 1;
                    bend = y1 + (y_multiplier * (j - numPaths_x + 1));
                    traverse_path1(x1, y1, x2, y2, path, bend, costs, min_costs, 
                                                                agg_costs, j);
                }
            }

            if (procID != root) {
                //fprintf(stdout, "%d %d(%d,%d)\t", numPaths, span, startIndex, endIndex);
                //if (startIndex < endIndex) {
                    MPI_Send(&min_costs[startIndex], endIndex - startIndex, MPI_INT, root, tag0, MPI_COMM_WORLD);
                    MPI_Send(&agg_costs[startIndex], endIndex - startIndex, MPI_INT, root, tag1, MPI_COMM_WORLD);
                //}
            }

            int index = 0;

            if (procID == root) {
                for (source = 1; source < nproc; source++) {
                    int startIndex1 = MIN(numPaths,source * span);
                    int endIndex1 = MIN(numPaths, startIndex1 + span);
                    //if (startIndex < endIndex) {
                    MPI_Recv(&min_costs[startIndex], endIndex1 - startIndex1, MPI_INT, source, tag0, MPI_COMM_WORLD, &status0);
                    MPI_Recv(&agg_costs[startIndex], endIndex1 - startIndex1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status1);
                    //}
                }

                int minCost = INT_MAX;
                for (int k = 0; k < numPaths; k++) {
                    if (min_costs[k] < minCost) {
                        minCost = min_costs[k];
                        index = k;
                    }
                    if (min_costs[k] == minCost) {
                        if (agg_costs[k] < agg_costs[index]) index = k;
                    }
                }
            }


            MPI_Bcast(&index, 1, MPI_INT, root, MPI_COMM_WORLD);

            //fprintf(stdout, "index = %d, procID = %d \n", index, procID);


            
            if (index < numPaths_x) {
                w.path = 0;
                w.bend = w.x1 + (x_multiplier *(index+ 1));
            }
            else {
                w.path = 1;
                w.bend = w.y1 + (y_multiplier * (index - numPaths_x + 1));
            }
        }
        else {
            int number = rand() % numPaths;
            if (number < numPaths_x) {
                w.path = 0;
                w.bend = w.x1 + x_multiplier*(number + 1);
            }
            else {
                w.path = 1;
                int offset = number - numPaths_x;
                w.bend = w.y1 + y_multiplier*(offset + 1);
            }
        }


        traverse_path(x1, y1, x2, y2, w.path, w.bend, costs, 1);   


    }

   // TODO Implement code here
    //TODO Decide which processors should be reading/writing files
    
    if (procID == root) {
       writeCost(inputFilename, nproc);
       writeOutput(inputFilename, nproc);  
    
        mC = 0;
        for (int m; m < size; m++) {
            int currCost = costs->array[m];
            if (currCost > mC) mC = currCost;
        }
        fprintf(stdout, "%d \n", mC);

    }

    free(costs->array);
    free(costs->transpose);
    free(costs);
    free(wires);
    free(min_costs);
    free(agg_costs);
}




// void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations)
// {
    

//     const int root = 0;
//     int tag0 = 0; 
//     int tag1 = 1;
//     MPI_Status status0, status1;
//     int source;
//     int choice;
//     int span, startIndex, endIndex;
//     int numPaths;
//     int numPaths_x;
//     int path, bend;
//     int x1, y1, x2, y2;
//     int x_multiplier, y_multiplier;
//     MPI_Datatype mpi_wire_type;
//     MPI_Datatype mpi_cost_type;
//     int mC;

//     readInput(inputFilename);
//     int size = costs->cols * costs->rows;


//     int *min_costs = (int *)calloc(costs->cols*costs->rows, sizeof(int));
//     int *agg_costs = (int *)calloc(costs->cols*costs->rows, sizeof(int));


//     for (int i = 0; i < num_wires; i++) {
//         wire_t *w = &wires[i];
//         x1 = w->x1;
//         y1 = w->y1;
//         x2 = w->x2;
//         y2 = w->y2;

//         //if (procID != root) fprintf(stdout, "%d %d %d %d \n", x1, y1, x2, y2);

//         x_multiplier = (x1 > x2) ? -1 : 1;
//         y_multiplier = (y1 > y2) ? -1 : 1;

//         numPaths = abs(x1 - x2) + abs(y1-y2);
//         numPaths_x = abs(x1 - x2);

//         //traverse_path(w->x1, w->y1, w->x2, w->y2, w->path, w->bend, costs, -1);

//         choice = 0;

//         if (choice == 0) {
//             span = numPaths/ nproc + 1;
//             startIndex = MIN(procID * span, numPaths);
//             endIndex = MIN(numPaths, startIndex + span);

//             for (int j = startIndex; j < endIndex; j++) {
//                 if (j < numPaths_x) {
//                     path = 0;
//                     bend = x1 + (x_multiplier * (j + 1));
//                     traverse_path1(x1, y1, x2, y2, path, bend, costs, min_costs, 
//                                                                 agg_costs, j);
//                 }
//                 else {
//                     path = 1;
//                     bend = y1 + (y_multiplier * (j - numPaths_x + 1));
//                     traverse_path1(x1, y1, x2, y2, path, bend, costs, min_costs, 
//                                                                 agg_costs, j);
//                 }
//             }

//             if (procID != root) {
//                 fprintf(stdout, "%d %d(%d,%d)\t", numPaths, span, startIndex, endIndex);
//                 if (startIndex < endIndex) {
//                     MPI_Send(&min_costs[startIndex], endIndex - startIndex, MPI_INT, root, tag0, MPI_COMM_WORLD);
//                     MPI_Send(&agg_costs[startIndex], endIndex - startIndex, MPI_INT, root, tag1, MPI_COMM_WORLD);
//                 }
//             }

//             int index = 0;

//             if (procID == root) {
//                 for (source = 1; source < nproc; source++) {
//                     int startIndex1 = MIN(source * span, numPaths);
//                     int endIndex1 = MIN(numPaths, startIndex1 + span);
//                     if (startIndex < endIndex) {
//                     MPI_Recv(&min_costs[startIndex], endIndex1 - startIndex1, MPI_INT, source, tag0, MPI_COMM_WORLD, &status0);
//                     MPI_Recv(&agg_costs[startIndex], endIndex1 - startIndex1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status0);
//                     }
//                 }

             
//                 int minCost = INT_MAX;
//                 for (int k = 0; k < costs->cols*costs->rows; k++) {
//                     if (min_costs[k] < minCost) {
//                         minCost = min_costs[k];
//                         index = k;
//                     }
//                     if (min_costs[k] == minCost) {
//                         if (agg_costs[k] < agg_costs[index]) index = k;
//                     }
//                 }
            
//                 if (index < numPaths_x) {
//                      w->path = 0;
//                      w->bend = w->x1 + (x_multiplier *(index+ 1));
//                 }
//                 else {
//                     w->path = 1;
//                     w->bend = w->y1 + (y_multiplier * (index - numPaths_x + 1));
//                 }
//             }
//         }
//         else {
//             if (procID == root) set_random_path(w);
//         }

//         if (procID == root) traverse_path(x1, y1, x2, y2, w->path, w->bend, costs, 1);   

//         MPI_Bcast(&(w->path), 1, MPI_INT, root, MPI_COMM_WORLD);
//         MPI_Bcast(&(w->bend), 1, MPI_INT, root, MPI_COMM_WORLD);
//         MPI_Bcast(costs->array, size, MPI_INT, root, MPI_COMM_WORLD);
//         MPI_Bcast(costs->transpose, size, MPI_INT, root, MPI_COMM_WORLD);
//         MPI_Barrier(MPI_COMM_WORLD);

//     }

//    // TODO Implement code here
//     //TODO Decide which processors should be reading/writing files
    
//     if (procID == root) {
//        writeCost(inputFilename, nproc);
//        writeOutput(inputFilename, nproc);  
    
//         mC = 0;
//         for (int m; m < size; m++) {
//             int currCost = costs->array[m];
//             if (currCost > mC) mC = currCost;
//         }
//         fprintf(stdout, "%d \n", mC);

//     }

//     free(costs->array);
//     free(costs->transpose);
//     free(costs);
//     free(wires);
//     free(min_costs);
//     free(agg_costs);
// }


// Read input file
void readInput(char* inputFilename)
{
    FILE* fp = fopen(inputFilename, "r");
    int numRows;
    int numCols;
    int numWires;
    int wireIndex;
    int x1;
    int y1;
    int x2;
    int y2;
    if (fp == NULL) {
        fprintf(stderr, "Failed to read input file %s\n", inputFilename);
        exit(-1);
    }
    if (fscanf(fp, "%d %d", &numCols, &numRows) != 2) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if ((numRows <= 0) || (numCols <= 0)) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (fscanf(fp, "%d", &numWires) != 1) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (numWires <= 0) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    init(numRows, numCols, numWires);
    for (wireIndex = 0; wireIndex < numWires; wireIndex++) {
        if (fscanf(fp, "%d %d %d %d", &x1, &y1, &x2, &y2) != 4) {
            fprintf(stderr, "Invalid input file format\n");
            exit(-1);
        }
        initWire(wireIndex, x1, y1, x2, y2);
    }
    fclose(fp);
}

// Write cost array file based on input filename
void writeCost(char* inputFilename, int nproc)
{
    char* dname = strdup(inputFilename);
    char* bname = strdup(inputFilename);
    char* costsFilename = malloc(strlen(inputFilename) + 100);
    FILE* fp;
    int numRows = getNumRows();
    int numCols = getNumCols();
    int row;
    int col;
    assert(costsFilename != NULL);
    sprintf(costsFilename, "%s/costs_%s_%d.txt", dirname(dname), basename(bname), nproc);
    fp = fopen(costsFilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Failed to write costs file %s\n", costsFilename);
        exit(-1);
    }
    fprintf(fp, "%d %d\n", numCols, numRows);
    
    for (row = 0; row < numRows; row++) {
        for (col = 0; col < numCols; col++) {
            fprintf(fp, "%d ", getCost(row, col));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    free(costsFilename);
    free(bname);
    free(dname);
}

// Write wire output file based on input filename
void writeOutput(char* inputFilename, int nproc)
{
    char* dname = strdup(inputFilename);
    char* bname = strdup(inputFilename);
    char* outputFilename = malloc(strlen(inputFilename) + 100);
    FILE* fp;
    int numRows = getNumRows();
    int numCols = getNumCols();
    int numWires = getNumWires();
    int wireIndex;
    int numPoints;
    int x1;
    int y1;
    int x2;
    int y2;
    int x3;
    int y3;
    int x4;
    int y4;
    assert(outputFilename != NULL);
    sprintf(outputFilename, "%s/output_%s_%d.txt", dirname(dname), basename(bname), nproc);
    fp = fopen(outputFilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Failed to write output file %s\n", outputFilename);
        exit(-1);
    }
    fprintf(fp, "%d %d\n", numCols, numRows);
    fprintf(fp, "%d\n", numWires);
    for (wireIndex = 0; wireIndex < numWires; wireIndex++) {
        numPoints = getWire(wireIndex, &x1, &y1, &x2, &y2, &x3, &y3, &x4, &y4);
        switch (numPoints) {
            case 2:
                fprintf(fp, "%d %d %d %d\n", x1, y1, x2, y2);
                break;

            case 3:
                fprintf(fp, "%d %d %d %d %d %d\n", x1, y1, x2, y2, x3, y3);
                break;

            case 4:
                fprintf(fp, "%d %d %d %d %d %d %d %d\n", x1, y1, x2, y2, x3, y3, x4, y4);
                break;

            default:
                assert(0); // invalid number of points
                break;
        }
    }
    fclose(fp);
    free(outputFilename);
    free(bname);
    free(dname);
}
