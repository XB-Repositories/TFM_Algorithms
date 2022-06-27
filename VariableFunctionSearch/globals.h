#ifndef LocalSearch_AllProblems_globals_h
#define LocalSearch_AllProblems_globals_h

extern char * Instance;
extern int PermuSize, Repetition, TabuSize, MaxEvals, Symmetry;
extern int ** flow_matrix;
extern int ** zero_flow_matrix;
extern int ** dist_matrix;
extern int ** zero_dist_matrix;
extern long int m_totalpq;
extern long int ** m_pq;
extern long int **** sum;
extern long int *** sum_exp;
extern int *** its;

#endif
