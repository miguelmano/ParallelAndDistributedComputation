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

pair compute_neighbors(char*** grid, long long n_cells, long long x, long long y, long long z, pair values);

void store_output_data(long long generation, long long species[]);

void simulation(char*** grid, char*** new_grid, long long n_cells, long long generation, long long max_generations);

int compute_value(int live_neighbors, int value, int majority);

void print_grid(char*** grid, long long n_cells);

char ***gen_copy_grid(long long N);

void print_result();

int is_power_two(long long n_cells);

#endif // LIFE3D_H
