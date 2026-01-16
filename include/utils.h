#ifndef UTILS_H
#define UTILS_H

#define DS_PI 3.14159265358979323846
#define INDEX(i, j, k) (NY * NZ * (i + 2 - i_start) + NZ * (j) + (k))
#define INDEX_F(i, j, k, p) (NY * NZ * NP * (i + 1 - i_start) + NZ * NP * (j) + NP * (k) + (p))

inline int mod(int x, int n)
{
    if (x < 0)
        return x + n;
    else if (x > n - 1)
        return x - n;
    else
        return x;
}

#define FOR_DOMAIN                        \
    for (int i = i_start; i < i_end; i++) \
        for (int j = 0; j < NY; j++)      \
            for (int k = 0; k < NZ; k++)

#define TIME(name, functions)                                        \
    if ((params->process_rank == 0) && (params->t_log == params->t)) \
    {                                                                \
        printf("%-70s", name);                                       \
    }                                                                \
    MPI_Barrier(MPI_COMM_WORLD);                                     \
    start_substep = MPI_Wtime();                                     \
    functions                                                        \
        MPI_Barrier(MPI_COMM_WORLD);                                 \
    duration_substep = MPI_Wtime() - start_substep;                  \
    if ((params->process_rank == 0) && (params->t_log == params->t)) \
    {                                                                \
        printf("[%7.4fs]\n", duration_substep);                      \
    }

#define TIME_OUTPUT(functions)                                                                \
    if ((params->process_rank == 0) && (params->t_log == params->t))                          \
    {                                                                                         \
        sprintf(output_info, "> \033[0;32mOutput to data_%d.h5...\033[0m", params->n_output); \
        printf("%-81s", output_info);                                                         \
    }                                                                                         \
    MPI_Barrier(MPI_COMM_WORLD);                                                              \
    start_substep = MPI_Wtime();                                                              \
    functions                                                                                 \
        MPI_Barrier(MPI_COMM_WORLD);                                                          \
    duration_substep = MPI_Wtime() - start_substep;                                           \
    if ((params->process_rank == 0) && (params->t_log == params->t))                                                            \
    {                                                                                         \
        printf("[%7.4fs]\n", duration_substep);                                               \
    }

#endif