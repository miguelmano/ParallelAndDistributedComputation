#include "world_gen.h"
#include "life3d-omp.h"

#include "omp.h"
#include <stdio.h>
#include <stdlib.h>

outputs data_output[N_SPECIES];
int CELLS_POWER_TWO = 0;

/// @brief Responsible for reading the data about the current generation, and store in case there is a maximum
/// @param generation Current generation 
/// @param species Array containing the data about the species on this generation
void store_output_data(long long generation, long long species[]) {

    if (generation == 0) {      // First generation, all values are max
        for (int i = 1; i < N_SPECIES + 1; i++) {
            data_output[i-1].generation = generation;
            data_output[i-1].quantity = species[i];
        }
    }
    else {
        for (int i = 1; i < N_SPECIES + 1; i++) {
            if (species[i] > data_output[i-1].quantity) {
                data_output[i-1].generation = generation;
                data_output[i-1].quantity = species[i];
            }
        }
    }
}


/// @brief Simulate current generation and store values
/// @param grid 3D World Cube
/// @param new_grid Grid to be changed during simulation
/// @param n_cells Number of cells per side of the cube 
/// @param generation Current Generation
/// @param max_generations Number of generations to be simulated
void simulation(char*** grid, char*** new_grid, long long n_cells, long long generation, long long max_generations) {   
    long long species[N_SPECIES + 1] = {0};

    if (generation > max_generations) {
        return;
    }

    pair values;
    long long i, j, k;

    #pragma omp parallel for collapse(2) schedule(dynamic) private(j,k,values) shared(grid,new_grid) reduction(+:species[:N_SPECIES+1])
    for (i = 0; i < n_cells; i++) {
        for (j = 0; j < n_cells; j++) {
            for (k = 0; k < n_cells; k++) {
                values = compute_neighbors(grid, n_cells, i, j, k, values);
                species[grid[i][j][k]]++;
                new_grid[i][j][k] = compute_value(values.live_neighbors, grid[i][j][k], values.majority);
            }
        }
    }
    
    store_output_data(generation, species);
    simulation(new_grid, grid, n_cells, generation + 1, max_generations);
}

/// @brief Responsible for receiving the number of live neighbors of a certain cell, and computing the value for its next generation
/// @param live_neighbors number of live neighbors
/// @param value current value of the cell
/// @param majority Majority of Scecies around cell
/// @return value of the cell in the next generation (-1 in case of the cell dying in the next generation)
int compute_value(int live_neighbors, int value, int majority) {

    if (value != 0){
        if (live_neighbors <= 4)
            return 0;
        else if (live_neighbors > 4 && live_neighbors <= 13)
            return value;
        else 
            return 0;
    }
    else if (live_neighbors >= 7 && live_neighbors <= 10){
        return majority;
    }

    return 0;

}


/// @brief Prints the content of the 3D Grid, just for debug purpose
/// @param grid 3D World Cube
/// @param n_cells Number of cells per side of the Cube
void print_grid(char*** grid, long long n_cells) {
    
    for (long long k = 0; k < n_cells; k++) {      
        printf("Layer %lld:\n", k);
        for (long long i = 0; i < n_cells; i++){
            for (long long j = 0; j < n_cells; j++) {
                if (grid[k][i][j] == 0)
                    printf("  ");
                else 
                    printf("%d ", grid[k][i][j]);
            }
            printf("\n");
        }
    }
}

/// @brief Prints output to stdout
void print_result() {

    for (int i = 0; i < N_SPECIES; i++) {
        fprintf(stdout, "%d %lld %lld\n", i+1, data_output[i].quantity, data_output[i].generation);
    }
}

/// @brief calculate the neighbors coordenates for every cell and check the majority of species and returns this informations
/// @param n_cells number of cells per side of the cube
/// @param x number of cells per side of the cube
/// @param y number of cells per side of the cube
/// @param z number of cells per side of the cube
/// @param values pair struct to save the results and output
/// @return pair values struct received as argument and changed 
pair compute_neighbors(char*** grid, long long n_cells, long long x, long long y, long long z, pair values) {

    long long n_x, n_y, n_z, neighbor_value;
    long long frequency[N_SPECIES+1] = {0};
    long long max = 0;
    int majority = 0;
    int i, j, k;

    for (i = -1; i <= 1; i++) {
        for (j = -1; j <= 1; j++) {
            for (k = -1; k <= 1; k++) {

                if (i == 0 && j == 0 && k == 0)
                    continue;

                if (CELLS_POWER_TWO) {
                    n_x = (x + i + n_cells) & (n_cells - 1);
                    n_y = (y + j + n_cells) & (n_cells - 1);
                    n_z = (z + k + n_cells) & (n_cells - 1);
                } else {
                    n_x = x + i;
                    n_y = y + j;
                    n_z = z + k;

                    if (n_x < 0)
                        n_x += n_cells;
                    else if (n_x >= n_cells)
                        n_x -= n_cells;

                    if (n_y < 0)
                        n_y += n_cells;
                    else if (n_y >= n_cells)
                        n_y -= n_cells;

                    if (n_z < 0)
                        n_z += n_cells;
                    else if (n_z >= n_cells)
                        n_z -= n_cells;
                }
                
                neighbor_value = grid[n_x][n_y][n_z];
                frequency[neighbor_value]++;
            }
        }
    }

    values.live_neighbors = N_NEIGHBORS - frequency[0];
    
    if (grid[x][y][z] != 0 || values.live_neighbors > 10 || values.live_neighbors < 7) {
        values.majority = 0;
        return values;
    }

    for (int i = 1; i <= N_SPECIES; i++) {
        if (frequency[i] > max) {
            majority = i;
            max = frequency[i];
        }
    }

    values.majority = majority;

    return values;
}

char ***gen_copy_grid(long long N)
{
    int x,y;
    char*** grid;

    grid = (char ***) malloc(N * sizeof(char **));
    if(grid == NULL) {
        printf("Failed to allocate copy matrix\n");
        exit(1);
    }
    
    for(x = 0; x < N; x++) {
        grid[x] = (char **) malloc(N * sizeof(char *));
        if(grid[x] == NULL) {
            printf("Failed to allocate copy matrix\n");
            exit(1);
        }
        grid[x][0] = (char *) calloc(N * N, sizeof(char));
        if(grid[x][0] == NULL) {
            printf("Failed to allocate copy matrix\n");
            exit(1);
        }
        for (y = 1; y < N; y++)
            grid[x][y] = grid[x][0] + y * N;
    }

    return grid;
}

int is_power_two(long long n_cells) {
    return ((n_cells & (n_cells - 1)) == 0);
}

int main(int argc, char *argv[]) {

    double exec_time;
    char*** grid;
    char*** new_grid;

    if (argc != 5) {
        fprintf(stderr, "Input Missing Format: ./life3d generations n_cells density input_seed\n");
    }

    int generations = atoi(argv[1]);
    long long n_cells = atoll(argv[2]);
    float density = atof(argv[3]);
    int input_seed = atoi(argv[4]);

    CELLS_POWER_TWO = is_power_two(n_cells);

    grid = gen_initial_grid(n_cells, density, input_seed);
    new_grid = gen_copy_grid(n_cells);
    
    exec_time = -omp_get_wtime();

    simulation(grid, new_grid, n_cells, 0, generations);

    exec_time += omp_get_wtime();
    
    fprintf(stderr, "%.1fs\n", exec_time);
    print_result();

    free(grid);
    free(new_grid);
    return 0;
}
