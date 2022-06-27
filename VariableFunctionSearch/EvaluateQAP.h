void initialize_vals();
void calculate_vals(int *genes);
void recalculate_vals(int *genes, int gene_1, int gene_2);
long double EvaluateQAP(int fitnessFunctionId, int * genes);
void EvaluateQAP_all(long double * fitness, int * genes);
long double EvaluateQAP_change(int fitnessFunctionId, long double fitness, int * genes, int * genes_swap, int gene_1, int gene_2);
void EvaluateQAP_change_all(long double * fitness, long double * fitness_new, int * genes, int * genes_swap, int gene_1, int gene_2);
