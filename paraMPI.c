#include <stdio.h>
#include <mpi.h>


double getGauge(double *p, int N0, int N1){
	const int i_gauge = (N0 -1)/2, j_gauge = (N1 -1)/4; /* at (0.5, 0.25) */
	return p[i_gauge + N0 * j_gauge];
}

void setInitialPressure(double *p, int N0, int s0, int N1){
	/* TODO had to make edit here to get j0 variable was: */
	/* const int i0 = (N0-1)/2, i1 = (N1-1)/2; */
	const int i0 = (N0-1)/2, j0 = (N1-1)/2; /* (0.5, 0.5) */
	printf("i0,j0: %d, %d\n", i0, j0);
	int i, j;
	for (j = 0; j < N1; j++) {
		for (i = 0; i < N0; i++){
			if (i == i0 && j == j0) {
				/* TODO why is there s0 instead of N0?? */
				p[i + s0 * j] = 1;  /* p = 1 at the centre of domain */
			} else {
				p[i + s0 * j] = 0; /* p = 0 elsewhere */
			}
		}

	}
}


void copyPressure(double *p_dst, const double *p_src, int N0, int N1){
	int i, j;
	for (j=0; j < N1; j++) {
		for (i=0; i < N0; i++) {
			p_dst[i + N0 * j] = p_src[i + N0 * j];
		}
	}
}



void updatePressure(double *p_new, const double *p, int N0, int N1,
					double a, double b) {
	int i, j;
	for (j = 1; j < N1 - 1; j++) {
		for (i = 1; i < N0 - 1; i++) {
			p_new[i + N0 * j] = 
					a * p[i + N0 * j]	+       	/* centre */
					b * (p[i - 1 + N0 * j] +		/* west */
						 p[i + 1 + N0 * j] + 		/* east */
						 p[i + N0 * (j - 1)] +		/* south */
						 p[i + N0 * (j + 1)]);		/* north */
		}
	}
}

void updateNorthPressure(double *p_new, const double *p, 
						const double *p_north, int N0,
						int N1_loc, double a, double b){
	int i; /* j = N1_loc -1 */
	for (j = 1; j < N0 - 1; j++){
		p_new[i + N0 * (N1_loc -1)] = 
				a * p[i + N0 * (N1_loc -1)]	+       	/* centre */
				b * (p[i - 1 + N0 * (N1_loc -1))] +		/* west */
					p[i + 1 + N0 * (N1_loc -1)] + 		/* east */
					p[i + N0 * (N1_loc - 2)] +		/* south */
					p_north[i]);		/* north */		
	}
}

void updateSouthPressure(double *p_new, const double *p, 
						const double *p_south, int N0,
						int N1_loc, double a, double b){
	int i; /* j = 0 */
	for (j = 1; j < N0 - 1; j++){
		p_new[i] = 
				a * p[i]		+      	/* centre */
				b * (p[i - 1] 	+		/* west */
					p[i + 1] 	+ 		/* east */
					p_south[i] 	+		/* south */
					p[i + N0]);			/* north */		
	}
}


void printGrid(double *p, int N0, int N1, FILE *fp){
	int i, j;
	for (j = 1; j < N1 - 1; j++) {
		for (i = 1; i < N0 - 1; i++) {
			fprintf(fp, "%.2lf ", p[i + N0 *j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}


int main(int argc, char **argv){
	MPI_Init(&argc, &argv);
	const int N0 = 51, N1 = N0; /* grid size, must be odd */
	/* Watch out for Courant condition (x2 resolution, x4 more time steps */
	const double dx = 1./(N0-1), dt = 1e-5, t_end = 0.5,
			a = 1-4 * dt/ (dx * dx), b = dt/(dx * dx);
	/* TODO probably want to malloc this?? */
	double t = 0., p_gauge, p[N0 * N1], p_new[N0 * N1], tmp[N0 * N1]; 
	int n = 0;
	setInitialPressure(p, N0, N0, N1);
	/* TODO Need to add an MPI_requet array thingy */
	int myrank = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status[2];
	MPI_Request requestR[2];
	MPI_Request requestS[2];

	while (t < t_end) {
		/* Sort out sending/reading data then do my own work */
		MPI_Irecv(p_south, N0, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, &requestR[0]);
		MPI_Irecv(p_north, N0, MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, &requestR[1]);
		MPI_Isend(p, N0, MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD, &requestS[0]);
		MPI_Isend(&p[N0 * (N1_loc, - 1)], N0, MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, &requestS[1]);
		
		updatePressure(p_new, p, N0, N1, a, b);
		
		/* Finished own work, need to wait for sides... */
		MPI_Wait(&requestR[0], &status[0]);
		updateSouthPressure(p_new, p, p_south, N0, N1_loc, a, b);

		MPI_Wait(&requestR[1], &status[1]);
		updateNorthPressure(p_new, p, p_north, N0, N1_loc, a, b);

		/* TODO what is this wait doing? waiting for our data to send i guess*/
		MPI_Waitall(2, requestS, status);

		p_gauge = getGauge(p_new, N0, N1);
		n++; t += dt;
		printf("time step: %d, t=%e, p at gauge = %e\n", n, t, p_gauge);

		/* TODO can just swap pointers here instead of copy */
		copyPressure(p, p_new, N0, N1);

	}
	MPI_Finalize();
	return 0;
}
