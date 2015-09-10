#include <stdio.h>


double getGauge(double *p, int N0, int N1){
	const int i_gauge = (N0 -1)/2, j_gauge = (N1 -1)/2; /* at (0.5, 0.25) */
	/* printf("Gauge at: %d, %d\n", i_gauge, j_gauge); */
	return p[i_gauge - 1 + N0 * j_gauge];
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
				p[i + N0 * j] = 1;  /* p = 1 at the centre of domain */
			} else {
				p[i + N0 * j] = 0; /* p = 0 elsewhere */
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
			/*if (p_new[i + N0 * j] > 0.) {
				printf("Setting (%d, %d) to %.2lf",i,j, p_new[i + N0 * j]);
			}*/
		}
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
	const int N0 = 51, N1 = N0; /* grid size, must be odd */
	/* Watch out for Courant condition (x2 resolution, x4 more time steps */
	const double dx = 1./(N0-1), dt = 1e-5, t_end = 0.5,
			a = 1-4 * dt/ (dx * dx), b = dt/(dx * dx);
	/* TODO probably want to malloc this?? */
	double t = 0., p_gauge, p[N0 * N1], p_new[N0 * N1], tmp[N0 * N1]; 
	int n = 0;
	setInitialPressure(p, N0, N0, N1);

	while (t < t_end) {
		updatePressure(p_new, p, N0, N1, a, b);
		p_gauge = getGauge(p_new, N0, N1);
		n++; t += dt;
		if (t < 0.1){
			printf("time step: %d, t=%e, p at gauge = %e\n", n, t, p_gauge);
		}
		/* TODO can just swap pointers here instead of copy */
		copyPressure(p, p_new, N0, N1);

	}

}
