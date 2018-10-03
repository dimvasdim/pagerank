#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define N 26518  //number of nodes
#define d 0.85
#define e 0.0001 // 10^-4

struct timeval startwtime, endwtime;
double seq_time;

void create_A(double **A, char *filename);

int main(int argc, char **argv)
{
  int i, j;
  double diff = 1;
  int iterations = 0;
  if (argc != 2)
  {
    printf("Need as an argument the txt file of the links between the nodes...\n");
    exit(1);
  }
  //Create A which is the stochastic table for our system.
  double **A = (double **) malloc (N * sizeof(double *));
  if (A == NULL)
  {
    printf("Couldn't allocate memory...\n");
    exit(1);
  }
  for (i=0; i<N; i++)
  {
    A[i] = (double *) malloc (N * sizeof(double));
    if (A[i] == NULL)
    {
      printf("Couldn't allocate memory...\n");
      exit(1);
    }
    //Initialize it as all zero
    for (j=0; j<N; j++)
    {
      A[i][j] = 0;
    }
  }
  create_A(A, argv[1]);

  double *x_prior = (double *) malloc (N * sizeof(double));
  if (x_prior == NULL)
  {
    printf("Couldn't allocate memory...\n");
    exit(1);
  }
  double *x_after = (double *) malloc (N * sizeof(double));
  if (x_after == NULL)
  {
    printf("Couldn't allocate memory...\n");
    exit(1);
  }
  // Initialize x_prior
  for (i=0; i<N; i++)
  {
    x_prior[i] = (double) 1 / N;
  }
  printf("\n");
  double b = (double) (1-d) / N;

  gettimeofday (&startwtime, NULL);
  while (diff > e)
  {
    diff = 0;
    for (i=0; i<N; i++)
    {
      x_after[i] = b;
      for (j=0; j<i; j++)
      {
        x_after[i] -= A[i][j] * x_after[j];
      }
      for (j=(i+1); j<N; j++)
      {
        x_after[i] -= A[i][j] * x_prior[j];
      }
      x_after[i] = (double) x_after[i] / A[i][i];
      diff += pow(x_after[i] - x_prior[i], 2);
      x_prior[i] = x_after[i];
    }
    diff = sqrt(diff);
    iterations++;
  }
  gettimeofday (&endwtime, NULL);
  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
                        + endwtime.tv_sec - startwtime.tv_sec);

  // Show the results
  printf("The result of the PageRank algorithm is:\n");
  for(i=(N-20); i<N; i++)
  {
    printf("%e  ", x_prior[i]);
  }
  double sum = 0;
  for (i=0; i<N; i++)
  {
    sum += x_prior[i];
  }
  printf("\nThe sum of the Pageranks is: %e\n", sum);
  printf("The iterations needed for the PageRank algorithm is: %d\n", iterations);
  printf("Time needed for PageRank algorithm is %f sec\n", seq_time);
  printf("\n \n");

  free(x_prior);
  free(x_after);
  for (i=0; i<N; i++)
  {
    free(A[i]);
  }
  free(A);
  return(0);
}


void create_A(double **A, char *filename)
{
  int i, j;
  FILE *fin;
  fin = fopen(filename, "r");
  int *count_links = (int *) malloc (N * sizeof(int));
  if (count_links == NULL)
  {
    printf("Couldn't allocate memory...\n");
    exit(1);
  }
  for (i=0; i<N; i++)
  {
    count_links[i] = 0;
  }
  while (!feof(fin))
  {
    if (fscanf(fin, "%d %d\n", &i, &j))
    {
      A[j][i] = 1;
      count_links[i]++;
    }
  }
  fclose(fin);
  for (j=0; j<N; j++)
  {
    if (count_links[j] > 1) // Greater than 1 cause if there is only 1 link there is no need to divide
    {
      for (i=0; i<N; i++)
      {
        A[i][j] = A[i][j] / count_links[j];
      }
    }
    else if (count_links[j] == 0)
    {
      for (i=0; i<N; i++)
      {
        A[i][j] = (double) 1 / N;
      }
    }
  }
  free(count_links);
  for (i=0; i<N; i++)
  {
    for (j=0; j<N; j++)
    {
      A[i][j] = -d*A[i][j];
    }
    A[i][i] += 1;
  }
}
