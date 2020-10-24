#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
//#include <boost/math/constants/constants.hpp>

//const double pi = boost::math::constants::pi<double>();

class Frame
{
public:
	Frame() {
		x = y = z = NULL;
	}
	Frame(int n) {
		x = new double[n];
		y = new double[n];
		z = new double[n];
	}
	~Frame() {
		if (!x) delete [] x;
		if (!y) delete [] y;
		if (!z) delete [] z;
	}
	double *x, *y, *z;
};

class Point
{
public:
	double x, y;
};

int main(int argc, char* argv[])
{
	if (argc<3) {
		printf("\n%s input.xyz output.dat [start_frame] [end_frame] [nq] [L]\n\n", argv[0]);
		return -1;
	}

	int i,j,k;
	int n_frames, i_frames, n_atoms, n_lips;
	int nqx, nqy;
	int st_frame, en_frame;
//  double min_wave, lbox;
	double *qxs, *qys, qx, qy, qt, qb, qmod;
	double lx, ly, lz;
	double **F, **N;
	double *xt, *xb, *yt, *yb, *zt, *zb, zm, factor;
	Point et, eb, **etx, **ety, **ebx, **eby, F_one, N_one;
	Frame frame;
	char *inpfile;
	char *outfile;
	
//	min_wave = 0.01; // Minimum wavelength in nm
	nqx = nqy = 100; // Maximum wave numbers
	// Change box length based on the simulation system
	lx = ly = lz = 250.0; // Box lengths
	st_frame = 0; // First frame to process
	en_frame = -1; // Last frame to process. -1 stands for last frame.
		
	inpfile = new char[strlen(argv[1])+1];
	outfile = new char[strlen(argv[2])+1];
	strcpy(inpfile, argv[1]);
	strcpy(outfile, argv[2]);
	if (argc>=4) st_frame = atoi(argv[3]);
	if (argc>=5) en_frame = atoi(argv[4]);
	if (argc>=6) nqx = nqy = atoi(argv[5]);
	if (argc>=7) lx = ly = lz = atof(argv[6]);
	
	// Pregenerate wave numbers
	qxs = new double[nqx+1];
	qys = new double[nqy+1];
	for (i=0;i<=nqx;++i) qxs[i] = 2.0*i*M_PI/lx;
	for (i=0;i<=nqy;++i) qys[i] = 2.0*i*M_PI/ly;
	
	// Creating and initializing arrays
	F = new double*[nqy+1];
	N = new double*[nqy+1];
	for (j=0;j<=nqy;++j) {
		F[j] = new double[nqx+1];
		N[j] = new double[nqx+1];
		for (i=0;i<=nqx;++i) {
			F[j][i] = 0.0;
			N[j][i] = 0.0;
		}
	}
	
	// Create intermediate arrays
	etx = new Point*[nqx];
	ety = new Point*[nqy];
	ebx = new Point*[nqx];
	eby = new Point*[nqy];
	
	// Open input file
	FILE *input = fopen(inpfile, "r");
	if (!input) {
		printf("\nERROR: Input file not found!\n\n");
		return 0;
	}
	
	int nstr, atom_id;
	char line[256], buff[256], *str[5];
	n_frames = 0;
	i_frames = 0;
	n_atoms = 0;
	while ( fgets(line, sizeof(line), input) != NULL ) {
		n_atoms = atoi(line);
		n_lips = n_atoms/2;
		
		// Skip one line
		fgets(line, sizeof(line), input);
		
		if (i_frames>=st_frame && (en_frame==-1 || i_frames<=en_frame)) {
		
			// Allocate memory
			if (frame.x==NULL) {
				frame.x = new double[n_atoms];
				frame.y = new double[n_atoms];
				frame.z = new double[n_atoms];
				
				xt = new double[n_lips];
				xb = new double[n_lips];
				yt = new double[n_lips];
				yb = new double[n_lips];
				zt = new double[n_lips];
				zb = new double[n_lips];
				
				for (i=0;i<=nqx;i++) {
					etx[i] = new Point[n_lips];
					ebx[i] = new Point[n_lips];
				}
				
				for (i=0;i<=nqy;i++) {
					ety[i] = new Point[n_lips];
					eby[i] = new Point[n_lips];
				}
			}
			
			if (i_frames%100==0) printf("Processing frame # %d\n", i_frames);
			
			// Read coordinates
			for (i=0;i<n_atoms;++i) {
				fgets(buff, sizeof(buff), input);
				
				// Split the line string and read the first 4 values
				nstr = 0;
				str[nstr] = strtok(buff," \t\n");
				while ( str[nstr]!=NULL ) {
					nstr++;
					if (nstr>=4) break;
					str[nstr] = strtok(NULL," \t\n");
				}
				if (nstr!=4) {
					printf("ERROR: Error reading XYZ file!\n");
					return 0;
				}
				
				atom_id = atoi(str[0]);
				frame.x[i] = atof(str[1]);
				frame.y[i] = atof(str[2]);
				frame.z[i] = atof(str[3]);
			}
			
			// Fourier transpose
			// Create xt, yt, zt, xb, yb, zb arrays;
			zm = 0.0;
			for (i=0;i<n_lips;i++) {
				xt[i] = frame.x[i];
				xb[i] = frame.x[i+n_lips];
				yt[i] = frame.y[i];
				yb[i] = frame.y[i+n_lips];
				zt[i] = frame.z[i];
				zb[i] = frame.z[i+n_lips];
				zm += zt[i] + zb[i];
			}
			zm /= n_atoms;
			
			for (k=0;k<n_lips;k++) {
				for (i=0;i<=nqx;i++) {
					qx = qxs[i];
					etx[i][k].x = cos(qx*xt[k]);
					etx[i][k].y = sin(qx*xt[k]);
					ebx[i][k].x = cos(qx*xb[k]);
					ebx[i][k].y = sin(qx*xb[k]);
				}
				for (i=0;i<=nqy;i++) {
					qy = qys[i];
					ety[i][k].x = cos(qy*yt[k]);
					ety[i][k].y = sin(qy*yt[k]);
					eby[i][k].x = cos(qy*yb[k]);
					eby[i][k].y = sin(qy*yb[k]);
				}
			}

			// Processing each frame
			factor = 0.25/(n_lips*n_lips);
			double Fx=0.0, Fy=0.0;
			for (i=0;i<=nqx;++i) {
//				qx = qxs[i];
				for (j=0;j<=nqy;++j) {
//					qy = qys[j];
					F_one.x = F_one.y = 0.0;
					N_one.x = N_one.y = 0.0;
					for (k=0;k<n_lips;++k) {
						et.x = etx[i][k].x * ety[j][k].x - etx[i][k].y * ety[j][k].y;
						et.y = -(etx[i][k].x * ety[j][k].y + etx[i][k].y * ety[j][k].x);
						eb.x = ebx[i][k].x * eby[j][k].x - ebx[i][k].y * eby[j][k].y;
						eb.y = -(ebx[i][k].x * eby[j][k].y + ebx[i][k].y * eby[j][k].x);
						
						F_one.x += (zt[k]-zm)*et.x + (zb[k]-zm)*eb.x;
						F_one.y += (zt[k]-zm)*et.y + (zb[k]-zm)*eb.y;

						N_one.x += et.x + eb.x;
						N_one.y += et.y + eb.y;
					}
					F[j][i] += factor*(F_one.x*F_one.x + F_one.y*F_one.y);
					N[j][i] += factor*(N_one.x*N_one.x + N_one.y*N_one.y);
				}
			}
			
			n_frames++;
		} else {
			// Skip n lines
			for (i=0;i<n_atoms;++i) fgets(buff, sizeof(buff), input);
		}
		i_frames++;
	}

	fclose(input);
	
	printf("n_lips=%d\n", n_lips);
	printf("n_frames=%d\n", n_frames);
	
	for (i=0;i<=nqx;++i) {
		for (j=0;j<=nqy;++j) {
			F[j][i] /= n_frames;
			N[j][i] /= n_frames;
		}
	}
	
	FILE *output = fopen(outfile, "w");
	for (i=0;i<=nqx;++i) {
		qx = qxs[i];
		for (j=0;j<=nqy;++j) {
			qy = qys[j];
			qmod = qx*qx + qy*qy;
			
			fprintf(output,"%f %f %f %f %f %f\n", F[j][i], qx, qy, N[j][i], qmod, sqrt(qmod));
		}
	}
	fclose(output);
	
	delete [] qxs;
	delete [] qys;
	for (j=0;j<=nqy;++j) {
		delete [] F[j];
		delete [] N[j];
	}
	for (i=0;i<=nqx;i++) {
		delete [] etx[i];
		delete [] ebx[i];
	}
	for (i=0;i<=nqy;i++) {
		delete [] ety[i];
		delete [] eby[i];
	}
	delete [] F;
	delete [] N;
	delete [] xt;
	delete [] xb;
	delete [] yt;
	delete [] yb;
	delete [] zt;
	delete [] zb;
	delete [] etx;
	delete [] ety;
	delete [] ebx;
	delete [] eby;
	
	return 0;
}
