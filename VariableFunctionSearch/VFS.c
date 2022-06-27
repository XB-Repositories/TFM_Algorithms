#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "VFS.h"
#include "list.h"
#include "EvaluateQAP.h"
#include "globals.h"


//Updates the best solution found according to f.
void check_best(strct_optimo* optimo, strct_optimo* actual){
	if(actual->opt_fitness < optimo->opt_fitness){
			optimo->opt_fitness=actual->opt_fitness;
			memcpy(optimo->opt_permu, actual->opt_permu, sizeof(int)*PermuSize);
	}
}

//Tabu exception (better than the best solution found).
int exception(strct_optimo* optimo, long double cost_neigh){
		if(cost_neigh < optimo->opt_fitness){
			return 1;
		}else{
			return 0;
		}
}

//Escapes from a local optimum using the objective functions of the decomposition (f_1,f_2,f_3).
int escape_search(int *evals, int max_evals, int iteration, strct_optimo* actual, strct_optimo* optimo, List* tabu, int tabu_size){
		long double cost_prev[4], cost_neigh[4], cost;
		int i, j, k, best_i, best_j, aux, cont;
		int *prev;

		prev=malloc(PermuSize*sizeof(int));
		memcpy(prev, actual->opt_permu, PermuSize*sizeof(int));

		EvaluateQAP_all(cost_prev, actual->opt_permu);

		cont=0;
		//Exploration of the swap neighborhood.
	  for (i=0; i<PermuSize-1; i++) {
				if(*evals>=max_evals) break;

	      for (j=i+1; j<PermuSize; j++) {
						if(*evals>=max_evals) break;

	          aux=actual->opt_permu[i];
	          actual->opt_permu[i]=actual->opt_permu[j];
	          actual->opt_permu[j]=aux;

						EvaluateQAP_change_all(cost_prev, cost_neigh, prev, actual->opt_permu, i, j);
						*evals+=1;

						//The algorithm moves to the best solution (according to f) that is better than the actual solution in at least one of the objective functions of the ELD.
						//Tabu movements are not considered.
	          if ((cost_neigh[0] < cost_prev[0] || cost_neigh[1] < cost_prev[1] || cost_neigh[2] < cost_prev[2]) && !inlist(*tabu,i,j)) {
								cost = cost_neigh[3];
								if(!cont || cost<actual->opt_fitness){
									best_i=i;
		              best_j=j;
		              actual->opt_fitness=cost;
									cont=1;
								}
	          }

	          aux= actual->opt_permu[j];
	          actual->opt_permu[j]=actual->opt_permu[i];
	          actual->opt_permu[i]=aux;
	      }
	  }
		//The tabu list is updated.
	  if (cont) {
	      aux=actual->opt_permu[best_i];
	      actual->opt_permu[best_i]=actual->opt_permu[best_j];
	      actual->opt_permu[best_j]=aux;
				recalculate_vals(actual->opt_permu,best_i,best_j);
				insertback(tabu,best_i,best_j);
				if(length(*tabu) > tabu_size) popfront(tabu);
	  }
		free(prev);
		return cont;
}

//One step of the main local search (based on f).
int local_search (int *evals, int max_evals, int iteration, strct_optimo* actual, strct_optimo* optimo, List* tabu, int tabu_size){
		long double cost_neigh, cost_prev;
		int i, j, best_i, best_j, aux, cont;
		int * prev;

		prev=malloc(PermuSize*sizeof(int));
		memcpy(prev, actual->opt_permu, PermuSize*sizeof(int));
		cost_prev=actual->opt_fitness;

		cont=0;
		//Exploration of the swap neighborhood.
	  for (i=0; i<PermuSize-1; i++) {
				if(*evals>=max_evals) break;

	      for (j=i+1; j<PermuSize; j++) {
						if(*evals>=max_evals) break;

	          aux=actual->opt_permu[i];
	          actual->opt_permu[i]=actual->opt_permu[j];
	          actual->opt_permu[j]=aux;

	          cost_neigh=EvaluateQAP_change(0, cost_prev, prev, actual->opt_permu, i, j);
						*evals+=1;

						//The algorithm moves to the best solution that is better than the actual solution according to f.
						//Tabu movements are not considered.
	          if (cost_neigh < actual->opt_fitness && (!inlist(*tabu,i,j) || exception(optimo,cost_neigh))) {
	              best_i=i;
	              best_j=j;
	              actual->opt_fitness=cost_neigh;
								cont=1;
	          }

	          aux= actual->opt_permu[j];
	          actual->opt_permu[j]=actual->opt_permu[i];
	          actual->opt_permu[i]=aux;
	      }
	  }
		//The tabu list is updated.
	  if (cont) {
	      aux=actual->opt_permu[best_i];
	      actual->opt_permu[best_i]=actual->opt_permu[best_j];
	      actual->opt_permu[best_j]=aux;
				recalculate_vals(actual->opt_permu,best_i,best_j);
				insertback(tabu,best_i,best_j);
				if(length(*tabu) > tabu_size) popfront(tabu);
				//The best solution found is updated.
				check_best(optimo, actual);
	  }
		free(prev);
		return cont;
}

//Variable Function Search.
void VFS(strct_optimo* optimo, int tabu_size, int max_evals){
		int evals,prev_evals,iteration,cont,cost;
		strct_optimo actual;
		List* tabu;

		tabu = malloc(sizeof(List));
		initlist(tabu);

		actual.opt_permu=malloc(PermuSize*sizeof(int));
		memcpy(actual.opt_permu, optimo->opt_permu, sizeof(int)*PermuSize);
		actual.opt_fitness=optimo->opt_fitness;

		calculate_vals(actual.opt_permu);

		iteration=1;
		evals=1;
		cont=1;

		//Main loop.
		while(evals < max_evals){
				//Main local search based on f.
				while(cont && evals < max_evals){
					prev_evals = evals;
					cont = local_search(&evals, max_evals, iteration, &actual, optimo, tabu, tabu_size);
					if(!cont){
						evals = prev_evals;
					}else{
						iteration++;
					}
				}
				//Escape from the local optimum using f_1, f_2, f_3.
				if(evals < max_evals){
					cont = escape_search(&evals, max_evals, iteration, &actual, optimo, tabu, tabu_size);
					iteration++;
				}
		}
}
