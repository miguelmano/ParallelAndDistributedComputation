#include "world_gen.h"
#include "life3d-mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

outputs data_output[N_SPECIES];
int CELLS_POWER_TWO = 0;

MPI_Comm MPI_COMM_CART;
MPI_Datatype MPI_GHOST_CELL;
MPI_Request send_request_horizontal[2];
MPI_Request send_request_vertical[2];
MPI_Request receive_request_horizontal[2];
MPI_Request receive_request_vertical[2];
MPI_Status status;

int process_up = 0;
int process_down = 0; 
int process_left = 0;
int process_right = 0;

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
void simulation(char*** grid, char*** new_grid, long long n_cells, long long generation, long long max_generations, int x_size, int y_size, int rank) {   
    long long species[N_SPECIES + 1] = {0};
    long long total_species[N_SPECIES + 1] = {0};

    if (generation > max_generations) {
        return;
    }

    pair values;
    long long x, y, z;

    ghost_edge up_ghost_edge;
    ghost_edge down_ghost_edge;
    ghost_edge left_ghost_edge;
    ghost_edge right_ghost_edge;

    /* HORIZONTAL EDGES */

    create_horizontal_ghost_edge(&up_ghost_edge, grid[1], y_size, n_cells);
    create_horizontal_ghost_edge(&down_ghost_edge, grid[(x_size-2)], y_size, n_cells);

    MPI_Isend(up_ghost_edge.cells, up_ghost_edge.alive_cells, MPI_GHOST_CELL, process_up, 0, MPI_COMM_CART, &send_request_horizontal[0]);
    MPI_Isend(down_ghost_edge.cells, down_ghost_edge.alive_cells, MPI_GHOST_CELL, process_down, 1, MPI_COMM_CART, &send_request_horizontal[1]);
    
    ghost_edge received_up_edge;
    ghost_edge received_down_edge;

    MPI_Probe(process_up, 1, MPI_COMM_CART, &status);
    MPI_Get_count(&status, MPI_GHOST_CELL, &received_up_edge.alive_cells);
    MPI_Probe(process_down, 0, MPI_COMM_CART, &status);
    MPI_Get_count(&status, MPI_GHOST_CELL, &received_down_edge.alive_cells);

    received_up_edge.cells = (cell *) calloc(received_up_edge.alive_cells, sizeof(cell));
    received_down_edge.cells = (cell *) calloc(received_down_edge.alive_cells, sizeof(cell));

    MPI_Irecv(received_up_edge.cells, received_up_edge.alive_cells, MPI_GHOST_CELL, process_up, 1, MPI_COMM_CART, &receive_request_horizontal[0]);
    MPI_Irecv(received_down_edge.cells, received_down_edge.alive_cells, MPI_GHOST_CELL, process_down, 0, MPI_COMM_CART, &receive_request_horizontal[1]);
    
    MPI_Waitall(2, receive_request_horizontal, MPI_STATUSES_IGNORE);
  
    add_horizontal_ghost_edge(&received_up_edge, grid[0]);
    add_horizontal_ghost_edge(&received_down_edge, grid[(x_size-1)]);

    /* VERTICAL EDGES */

    create_vertical_ghost_edge(&left_ghost_edge, grid, x_size, 1, n_cells);
    create_vertical_ghost_edge(&right_ghost_edge, grid, x_size, (y_size-2), n_cells);

    MPI_Isend(left_ghost_edge.cells, left_ghost_edge.alive_cells, MPI_GHOST_CELL, process_left, 2, MPI_COMM_CART, &send_request_vertical[0]);
    MPI_Isend(right_ghost_edge.cells, right_ghost_edge.alive_cells, MPI_GHOST_CELL, process_right, 3, MPI_COMM_CART, &send_request_vertical[1]);
    
    ghost_edge received_left_edge;
    ghost_edge received_right_edge;

    MPI_Probe(process_left, 3, MPI_COMM_CART, &status);
    MPI_Get_count(&status, MPI_GHOST_CELL, &received_left_edge.alive_cells);
    MPI_Probe(process_right, 2, MPI_COMM_CART, &status);
    MPI_Get_count(&status, MPI_GHOST_CELL, &received_right_edge.alive_cells);

    received_left_edge.cells = (cell *) calloc(received_left_edge.alive_cells, sizeof(cell));
    received_right_edge.cells = (cell *) calloc(received_right_edge.alive_cells, sizeof(cell));

    MPI_Irecv(received_left_edge.cells, received_left_edge.alive_cells, MPI_GHOST_CELL, process_left, 3, MPI_COMM_CART, &receive_request_vertical[0]);
    MPI_Irecv(received_right_edge.cells, received_right_edge.alive_cells, MPI_GHOST_CELL, process_right, 2, MPI_COMM_CART, &receive_request_vertical[1]);

    MPI_Waitall(2, receive_request_vertical, MPI_STATUSES_IGNORE);

    add_vertical_ghost_edge(&received_left_edge, grid, 0);
    add_vertical_ghost_edge(&received_right_edge, grid, (y_size-1));

    #pragma omp parallel for collapse(2) schedule(dynamic) private(y,z,values) shared(grid,new_grid) reduction(+:species[:N_SPECIES+1])
    for (x = 1; x < (x_size-1); x++) {
        for (y = 1; y < (y_size-1); y++) {
            for (z = 0; z < n_cells; z++) { 
                values = compute_neighbors(grid, n_cells, x, y, z, values);
                species[grid[x][y][z]]++;
                new_grid[x][y][z] = compute_value(values.live_neighbors, grid[x][y][z], values.majority);
            }
        }
    }
    
    // TODO: figure out why these give SegFault
    //free(up_ghost_edge.cells); 
    //free(down_ghost_edge.cells);
    //free(left_ghost_edge.cells);
    //free(right_ghost_edge.cells);
    
    free(received_up_edge.cells);
    free(received_down_edge.cells);
    free(received_left_edge.cells);
    free(received_right_edge.cells);

    clear_horizontal_ghost_edge(grid[0], y_size, n_cells);
    clear_horizontal_ghost_edge(grid[(x_size-1)], y_size, n_cells);
    clear_vertical_ghost_edge(grid, x_size, 0, n_cells);
    clear_vertical_ghost_edge(grid, x_size, (y_size-1), n_cells);

    MPI_Reduce(species, total_species, (N_SPECIES + 1), MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_CART);

    if (rank == 0)
        store_output_data(generation, total_species);
    
    simulation(new_grid, grid, n_cells, generation + 1, max_generations, x_size, y_size, rank);
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
void print_grid(char*** grid, int x_size, int y_size, long long n_cells) {
    
    for (long long k = 0; k < x_size; k++) {      
        printf("Layer %lld:\n", k);
        for (long long i = 0; i < y_size; i++){
            for (long long j = 0; j < n_cells; j++) {
                if (grid[k][i][j] == 0)
                    printf(". ");
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

                n_x = x + i;
                n_y = y + j;

                if (CELLS_POWER_TWO) {
                    n_z = (z + k + n_cells) & (n_cells - 1);
                } else {
                    n_z = z + k;

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

char ***gen_copy_grid(partition* part, long long N) {
    int x,y;
    char*** grid;

    grid = (char ***) malloc(part->x_size * sizeof(char **));
    if(grid == NULL) {
        printf("Failed to allocate copy matrix\n");
        exit(1);
    }
    
    for(x = 0; x < part->x_size; x++) {
        grid[x] = (char **) malloc(part->y_size * sizeof(char *));
        if(grid[x] == NULL) {
            printf("Failed to allocate copy matrix\n");
            exit(1);
        }
        grid[x][0] = (char *) calloc(part->y_size * N, sizeof(char));
        if(grid[x][0] == NULL) {
            printf("Failed to allocate copy matrix\n");
            exit(1);
        }
        for (y = 1; y < part->y_size; y++)
            grid[x][y] = grid[x][0] + y * N;
    }
    return grid;
}

int is_power_two(long long n_cells) {
    return ((n_cells & (n_cells - 1)) == 0);
}

void calculate_partition(partition* part, int n_cells, int* dims, int* rank_coords) {
    part->x_size = n_cells / dims[0];
    part->y_size = n_cells / dims[1];

	int x_remainder = n_cells % dims[0];
	int y_remainder = n_cells % dims[1];
	
    part->first_x = part->x_size * rank_coords[0] + x_remainder;
	part->first_y = part->y_size * rank_coords[1] + y_remainder;
    part->last_x = part->first_x + part->x_size - 1;
    part->last_y = part->first_y + part->y_size - 1;

	if (rank_coords[0] == 0) { // Process x is 0
		part->x_size += x_remainder;
		part->first_x -= x_remainder;
	}
	if (rank_coords[1] == 0) { // Process y is 0
		part->y_size += y_remainder;
		part->first_y -= y_remainder;
	}
	
    // Extra space for ghost edges
    part->x_size += 2;
    part->y_size += 2;
}

void create_horizontal_ghost_edge(ghost_edge* horizontal_edge, char **grid_edge, int y_size, long long n_cells) {
    int y, z;
    long long alive_cells = 0;
    long long ind = 0;

    for (y = 1; y < (y_size - 1); y++) {
        for (z = 0; z < n_cells; z++) {
            if (grid_edge[y][z])
                alive_cells++;
        }
    }

    horizontal_edge->alive_cells = alive_cells;

    if (alive_cells == 0) {
        return;
    }

    horizontal_edge->cells = (cell *) calloc(alive_cells, sizeof(cell));

    for (y = 1; y < (y_size - 1); y++) {
        for (z = 0; z < n_cells; z++) {
            if (grid_edge[y][z]) {
                horizontal_edge->cells[ind].i = y;
                horizontal_edge->cells[ind].z = z;
                horizontal_edge->cells[ind].val = grid_edge[y][z];

                ind++;
            }
        }
    }
}

void create_vertical_ghost_edge(ghost_edge* vertical_edge, char ***grid_edge, int x_size, int y_pos, long long n_cells) {
    int x, z;
    long long alive_cells = 0;
    long long ind = 0;

    for (x = 0; x < x_size; x++) {
        for (z = 0; z < n_cells; z++) {
            if (grid_edge[x][y_pos][z])
                alive_cells++;
        }
    }

    vertical_edge->alive_cells = alive_cells;

    if (alive_cells == 0) {
        return;
    }

    vertical_edge->cells = (cell *) calloc(alive_cells, sizeof(cell));

    for (x = 0; x < x_size; x++) {
        for (z = 0; z < n_cells; z++) {
            if (grid_edge[x][y_pos][z]) {

                vertical_edge->cells[ind].i = x;
                vertical_edge->cells[ind].z = z;
                vertical_edge->cells[ind].val = grid_edge[x][y_pos][z];

                ind++;
            }
        }
    }
}

void add_horizontal_ghost_edge(ghost_edge* horizontal_edge, char **grid) {
    int ind;

    for (ind = 0; ind < horizontal_edge->alive_cells; ind++) {
        grid[horizontal_edge->cells[ind].i][horizontal_edge->cells[ind].z] = horizontal_edge->cells[ind].val;
    }
}

void add_vertical_ghost_edge(ghost_edge* vertical_edge, char ***grid, int y_pos) {
    int ind;

    for (ind = 0; ind < vertical_edge->alive_cells; ind++) {
        grid[vertical_edge->cells[ind].i][y_pos][vertical_edge->cells[ind].z] = vertical_edge->cells[ind].val;
    }
}

void clear_horizontal_ghost_edge(char **grid, int y_size, long long n_cells) {
    int y, z;

    for (y = 0; y < y_size; y++) {
        for (z = 0; z < n_cells; z++) {
            grid[y][z] = 0;
        }
    }
}

void clear_vertical_ghost_edge(char ***grid, int x_size, int y_pos, long long n_cells) {
    int x, z;

    for (x = 1; x < (x_size-1); x++) {
        for (z = 0; z < n_cells; z++) {
            grid[x][y_pos][z] = 0;
        }
    }
}

int main(int argc, char *argv[]) {

    double exec_time;
    char*** grid;
    char*** new_grid;

    int rank, procs;
    int dims[2] = {0, 0};
    int periods[2] = {1, 1};
    int rank_coords[2] = {0, 0}; 
    partition part;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc != 5) {
        fprintf(stderr, "Input Missing Format: ./life3d generations n_cells density input_seed\n");
    }

    int generations = atoi(argv[1]);
    long long n_cells = atoll(argv[2]);
    float density = atof(argv[3]);
    int input_seed = atoi(argv[4]);

    MPI_Type_contiguous(3, MPI_INT, &MPI_GHOST_CELL);
	MPI_Type_commit(&MPI_GHOST_CELL);

    MPI_Dims_create(procs, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &MPI_COMM_CART);
    MPI_Cart_coords(MPI_COMM_CART, rank, 2, rank_coords);
    MPI_Cart_shift(MPI_COMM_CART, 0, 1, &process_up, &process_down);
	MPI_Cart_shift(MPI_COMM_CART, 1, 1, &process_left, &process_right);

    CELLS_POWER_TWO = is_power_two(n_cells);
    
    calculate_partition(&part, n_cells, dims, rank_coords);
    grid = gen_initial_grid(n_cells, density, input_seed, part.first_x, part.first_y, part.last_x, part.last_y, part.x_size, part.y_size);
    new_grid = gen_copy_grid(&part, n_cells);

    exec_time = -omp_get_wtime();
    
    simulation(grid, new_grid, n_cells, 0, generations, part.x_size, part.y_size, rank);

    exec_time += omp_get_wtime();

    if (rank == 0) {
        fprintf(stderr, "%.1fs\n", exec_time);
        print_result();
    }

    free(grid);
    free(new_grid);

    MPI_Finalize();
    return 0;
}
