#ifndef LIFE3D_H
#define LIFE3D_H

#define N_NEIGHBORS 26

typedef struct {
    int live_neighbors;
    int majority;
} pair;

typedef struct {
    long long quantity;
    long long generation;
} outputs;

typedef struct {
    int first_x;
    int first_y;
    int last_x;
    int last_y;
    int x_size;
    int y_size;
} partition;

typedef struct {
    int i;
    int z;
    int val;
} cell;

typedef struct {
    cell *cells;
    int alive_cells;
} ghost_edge;

pair compute_neighbors(char*** grid, long long n_cells, long long x, long long y, long long z, pair values);

void store_output_data(long long generation, long long species[]);

void simulation(char*** grid, char*** new_grid, long long n_cells, long long generation, long long max_generations, int x_size, int y_size, int rank);

int compute_value(int live_neighbors, int value, int majority);

void print_grid(char*** grid, int x_size, int y_size, long long n_cells);
                            
char ***gen_copy_grid(partition* part, long long N);

void print_result();

int is_power_two(long long n_cells);

void calculate_partition(partition* part, int n_cells, int* dims, int* rank_coords);

void create_horizontal_ghost_edge(ghost_edge* horizontal_edge, char **grid_edge, int y_size, long long n_cells);

void create_vertical_ghost_edge(ghost_edge* vertical_edge, char ***grid_edge, int x_size, int y_pos, long long n_cells);

void add_horizontal_ghost_edge(ghost_edge* horizontal_edge, char **grid);

void add_vertical_ghost_edge(ghost_edge* vertical_edge, char ***grid, int y_pos);

void clear_horizontal_ghost_edge(char **grid, int y_size, long long n_cells);

void clear_vertical_ghost_edge(char ***grid, int x_size, int y_pos, long long n_cells);

#endif // LIFE3D_H
