#ifndef WORLD_GENERATOR_H
#define WORLD_GENERATOR_H

#define N_SPECIES 9

void init_r4uni(int input_seed);

float r4uni();

char ***gen_initial_grid(long long N, float density, int input_seed, int first_x, int first_y, int last_x, int last_y, int x_size, int y_size);

#endif