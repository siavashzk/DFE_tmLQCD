#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

#define TEST_MAIN
#include "su3.h"
#include "global.h"
#include "geometry.h"
#include "hopping.h"

void create_random_input(spinor* s, su3* u);
void create_random_spinor(spinor * s);
void create_random_su3vector(su3_vector *v);
void create_random_su3(su3 *s);
void create_random_complex(float complex *a);
void add_4d_halos_spinor(spinor* with_halos, spinor* orig, int halos, int eo);
void add_4d_halos_gauge(su3* with_halos, su3* orig, int halos, int eo) ;
int AreNotSameComplex(complex float a, float complex b);
int compare_spinor (spinor *a, spinor *b);
void read_spinor(char * filename, spinor *out);
void read_gauge(char * filename, su3 *s, int V);
void reorganize_ueven (su3 *out, su3 *in);
void reorganize_back_ueven (su3 *out, su3 *in);
void devide_gauge_to_oddeven(su3 const * const in, su3 * const even, su3 * const odd);

int main(void) {

	L = 4;
	T = 4;
	VOLUME = L * L * L * T;
	LX = LY = LZ = L;
	//int VOLUMEH1 = (L/2+4) * (L+8) * (L+8) * (T+8);
	//int VOLUMEH2 = (L/2+3) * (L+6) * (L+6) * (T+6);
	//int VOLUMEH3 = (L/2+2) * (L+4) * (L+4) * (T+4);
	int VOLUMEH1 = (L/2+1) * (L+2) * (L+2) * (T+2);

	spinor *in, *out, *out_dfe, *out_expected;
	su3 *uodd, *ueven, *u, *ueven_;

	ka3 = ka2 = ka1 = ka0 = 1;

	in = malloc(VOLUME * sizeof(spinor));
	out = &in[VOLUME / 2];
	out_dfe = malloc(VOLUMEH1 * sizeof(spinor));
	out_expected = malloc(VOLUME/2 * sizeof(spinor));
	ueven = malloc(VOLUME * 4 * sizeof(su3));
	uodd = &ueven[VOLUME / 2 * 4];
	ueven_ = malloc(VOLUME/2 * 4 * sizeof(su3));
	u = malloc (VOLUME / 2 * 8 * sizeof(su3));

	/*for (int i=0 ; i<VOLUME/2 ; i++ ) {
		create_random_spinor(in + i);
	}
	for (int i=0 ; i<VOLUME*2 ; i++ ) {
		create_random_su3(ueven + i);
		create_random_su3(uodd + i);
	}*/
	read_spinor("in_spinor.txt", in);
	read_spinor("out_spinor.txt", out_expected);
	read_gauge("in_gauge0.txt", u, VOLUME/2);

	devide_gauge_to_oddeven(u, ueven_, uodd);
	reorganize_back_ueven(ueven, ueven_);
	//reorganize_ueven(ueven_, ueven);

	spinor *in_halos = malloc(VOLUMEH1 * sizeof(spinor));
	//spinor *out_halos = malloc(VOLUMEH2 * sizeof(spinor));
	su3 *uodd_halos  = malloc(VOLUMEH1 * 4 * sizeof(spinor));
	su3 *ueven_halos_  = malloc(VOLUMEH1 * 4 * sizeof(spinor));
	//su3 *ueven_halos  = malloc(VOLUME * 4 * sizeof(spinor));

	add_4d_halos_spinor(in_halos, in, 1, 0);
	add_4d_halos_gauge(uodd_halos, uodd, 1, 0);
	add_4d_halos_gauge(ueven_halos_, ueven_, 1, 0);
	//add_4d_halos_gauge(ueven_halos, ueven, 3, 1);

	/*int lhx, lhy, lhz, h;
	h = 1;
	lhx = lhy = L+2*h;
	lhz = L/2+h;
	for (int i=0 ; i<2 ; i++ ) {
		spinor *a = in + i;
		printf("%f %f \n", creal(a->s0.c0), cimag(a->s0.c0) );

	}
	for (int i=0 ; i<5 ; i++ ) {
		spinor *a = in_halos + h*lhx*lhy*lhz +
					           h*lhy*lhz +
					           h*lhz + i;
		printf("%f %f \n", creal(a->s0.c0), cimag(a->s0.c0) );

	}*/

	init_geometry();

	//tm_times_Hopping_Matrix(out, in, uodd, ueven, 0.5+I*0, 0);
	//tm_sub_Hopping_Matrix(out, in2, in, uodd, ueven, 0.5+I*0, 1);
	Qtm_pm_psi(out, in, uodd, ueven, 1 + 0*I, 1 + 0*I);

	max_file_t *maxfile = S_LQCD_init();
	max_engine_t *engine = max_load(maxfile, "*");

	max_actions_t* act = max_actions_init(maxfile, "default");

	max_set_param_double(act, "ka", 1);
	max_set_param_double(act, "cfactor", 1);

	max_queue_input(act, "times1kernel_spinor_in", in_halos, VOLUMEH1 * sizeof(spinor));
	max_queue_input(act, "times1kernel_gauge0", uodd_halos, VOLUMEH1 * 4 * sizeof(su3));
	max_queue_input(act, "times1kernel_gauge1", ueven_halos_, VOLUMEH1 * 4 * sizeof(su3));

	max_queue_output(act, "times1kernel_spinor_out", out_dfe,  VOLUME/2 * sizeof(spinor));
	//max_queue_output(act, "sub2kernel_spinor_out", out_dfe,  VOLUME/2 * sizeof(spinor));

	printf("Running on DFE (mode: ComputeWithScalar)...\n");
	max_run(engine, act);
	max_unload(engine);
	printf("Done\n");

	//add_4d_halos_spinor(out_halos, out, 3, 1);
	for (int i=0 ; i<VOLUME/2 ; i++ ) {
		int error = compare_spinor(out+i , out_dfe+i);
		if (error) {
			printf("Wrong %d! %d\n", error,i);
			spinor *a = out + i;
			spinor *b = out_dfe + i;
			printf("%f %f    %f %f ", creal(a->s0.c0), cimag(a->s0.c0),
					                  creal(b->s0.c0), cimag(b->s0.c0));
			return 1;
		}
	}
	printf("Good.\n");

	return 0;
}

void print_cmplx(float * a, int N) {
	for (int i = 0; i < N; i += 2) {
		printf("(%f,%f i)\n", a[i], a[i + 1]);
	}
}

void create_random_input(spinor* s, su3* u) {
	for (int i = 0; i < VOLUME / 2; i++) {
		create_random_spinor(s + i);
	}
	for (int i = 0; i < VOLUME * 4; i++) {
		create_random_su3(u + i);
	}
}

void create_random_spinor(spinor * s) {
	create_random_su3vector(&s->s0);
	create_random_su3vector(&s->s1);
	create_random_su3vector(&s->s2);
	create_random_su3vector(&s->s3);
}

void create_random_su3vector(su3_vector *v) {
	create_random_complex(&v->c0);
	create_random_complex(&v->c1);
	create_random_complex(&v->c2);
}

void create_random_su3(su3 *s) {
	create_random_complex(&s->c00);
	create_random_complex(&s->c01);
	create_random_complex(&s->c02);
	create_random_complex(&s->c10);
	create_random_complex(&s->c11);
	create_random_complex(&s->c12);
	create_random_complex(&s->c20);
	create_random_complex(&s->c21);
	create_random_complex(&s->c22);
}

void create_random_complex(float complex *a) {
	float r[2];
	for (int i = 0; i < 2; i++) {
		r[i] = (float) (rand()) / RAND_MAX * 2;
	}
	*a = r[0] + I * r[1];
}

void add_4d_halos_spinor(spinor* with_halos, spinor* orig, int halos, int eo) {
	int lhx, lhy, lhz;
	lhx = lhy = L+2*halos;
	lhz = L/2+halos;

	eo = eo ^ 1;

	for (int t = -halos ; t < T+halos ; t++ ) {
		for (int x = -halos ; x < L+halos ; x++ ) {
			for (int y = -halos ; y < L+halos ; y++ ) {
				for (int z = -(halos/2) ; z < L/2+(halos+1)/2 ; z++ ) {
					int tt = (t+T)%T;
					int xx = (x+L)%L;
					int yy = (y+L)%L;
					int isOddRow = (tt & 1) ^ (xx & 1) ^ (yy & 1) ^ eo;
					int zz;
					if (halos%2 == 0) {
						zz = (z+L/2)%(L/2);
					} else {
						zz = (z+L*2-isOddRow*halos)%(L/2);
					}

					with_halos[ (t+halos)*lhx*lhy*lhz +
					            (x+halos)*lhy*lhz +
					            (y+halos)*lhz + (z+halos/2)] =
							orig[ tt*L*L*L/2 + xx*L*L/2 + yy*L/2 + zz];
				}
			}
		}
	}
}

void add_4d_halos_gauge(su3* with_halos, su3* orig, int halos, int eo) {
	int lhx, lhy, lhz;
	lhx = lhy = L+2*halos;
	lhz = L/2+halos;

	eo = eo ^ 1;

	for (int t = -halos ; t < T+halos ; t++ ) {
		for (int x = -halos ; x < L+halos ; x++ ) {
			for (int y = -halos ; y < L+halos ; y++ ) {
				for (int z = -(halos/2) ; z < L/2+(halos+1)/2 ; z++ ) {
					int tt = (t+T)%T;
					int xx = (x+L)%L;
					int yy = (y+L)%L;
					int isOddRow = (tt & 1) ^ (xx & 1) ^ (yy & 1) ^ eo;
					int zz;
					if (halos%2 == 0) {
						zz = (z+L/2)%(L/2);
					} else {
						zz = (z+L*2-isOddRow*halos)%(L/2);
					}

					for (int i=0 ; i < 4 ; i++ ){
						with_halos[ (t+halos)*lhx*lhy*lhz*4 +
						            (x+halos)*lhy*lhz*4 +
						            (y+halos)*lhz*4 +
						            (z+halos/2)*4 + i] =
						      orig[ tt*L*L*L/2*4 + xx*L*L/2*4 + yy*L/2*4 + zz*4 + i];
					}
				}
			}
		}
	}
}

int AreNotSameComplex(complex float a, float complex b)
{
    return ( fabs(creal(a) - creal(b)) > 0.001 ) ||
    	   ( cimag(creal(a) - cimag(b)) > 0.001 );
}

int compare_spinor (spinor *a, spinor *b) {
	if (AreNotSameComplex(a->s0.c0 ,b->s0.c0)) return 1;
	if (AreNotSameComplex(a->s0.c1 ,b->s0.c1)) return 2;
	if (AreNotSameComplex(a->s0.c2 ,b->s0.c2)) return 3;
	if (AreNotSameComplex(a->s1.c0 ,b->s1.c0)) return 4;
	if (AreNotSameComplex(a->s1.c1 ,b->s1.c1)) return 5;
	if (AreNotSameComplex(a->s1.c2 ,b->s1.c2)) return 6;
	if (AreNotSameComplex(a->s2.c0 ,b->s2.c0)) return 7;
	if (AreNotSameComplex(a->s2.c1 ,b->s2.c1)) return 8;
	if (AreNotSameComplex(a->s2.c2 ,b->s2.c2)) return 9;
	if (AreNotSameComplex(a->s3.c0 ,b->s3.c0)) return 10;
	if (AreNotSameComplex(a->s3.c1 ,b->s3.c1)) return 11;
	if (AreNotSameComplex(a->s3.c2 ,b->s3.c2)) return 12;
	return 0;
}

void read_spinor(char * filename, spinor *out) {
	FILE * fptr = fopen(filename, "r" );
	char temp[64];


	for (int t=0 ; t<T ; t++ ) {
		for (int z=0 ; z<L ; z++ ) {
			for (int y=0 ; y<L ; y++ ) {
				for (int x=0 ; x<L/2 ; x++ ) {

					int isOddRow = (t & 1) ^ (z & 1) ^ (y & 1);
					int tt = t;               // converting from checkerboarded
					int zz = z/2;             // coordinates of qphix along x-axis
					int yy = y;               // to tmLQCD checkerboarding along
					int xx = (2*x)+isOddRow;  // along y-axis

					spinor *s = &out [ (((((tt*L)+xx)*L)+yy)*L/2)+zz ];

					fgets (temp, 64, fptr);
					float a, b;

					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s0.c0 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s0.c1 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s0.c2 = a + I * b;

					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s1.c0 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s1.c1 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s1.c2 = a + I * b;

					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s2.c0 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s2.c1 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s2.c2 = a + I * b;

					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s3.c0 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s3.c1 = a + I * b;
					fscanf(fptr, "%f\n%f\n", &a, &b);
					s->s3.c2 = a + I * b;

				}
			}
		}
	}
	fclose(fptr);
}

void read_gauge(char * filename, su3 *s, int V) {
	FILE * fptr = fopen(filename, "r" );
	char temp[64];
	for (int i = 0 ; i < V ; i++ ) {
		fgets (temp, 64, fptr);
		float a, b;

		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c00 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c01 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c02 = a + I * b;

		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c10 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c11 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c12 = a + I * b;

		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c20 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c21 = a + I * b;
		fscanf(fptr, "%f\n%f\n", &a, &b);
		s->c22 = a + I * b;

		s++;
	}
	fclose(fptr);
}

void reorganize_ueven (su3 *out, su3 *in) {
	int i = 0;
	for (int t=0 ; t<T ; t++ ) {
		for (int x=0 ; x<L ; x++ ) {
			for (int y=0 ; y<L ; y++ ) {
				for (int z=0 ; z<L/2 ; z++ ) {
					int isOddRow = (t & 1) ^ (x & 1) ^ (y & 1);
					for (int mu=0; mu<4 ; mu++ ) {
						int tt = (mu==0)?t+1:t;
						int xx = (mu==1)?x+1:x;
						int yy = (mu==2)?y+1:y;
						int zz = (mu==3)?( (isOddRow)?(z):z+1 ):z;

						tt = (tt+T) % T;
						xx = (xx+L) % L;
						yy = (yy+L) % L;
						zz = (zz+L) % (L/2);

						out[i] = in [ ((((((tt*L)+xx)*L)+yy)*L/2)+zz)*4+mu ];
						i++;
					}
				}
			}
		}
	}
}

void reorganize_back_ueven (su3 *out, su3 *in) {
	int i = 0;
	for (int t=0 ; t<T ; t++ ) {
		for (int x=0 ; x<L ; x++ ) {
			for (int y=0 ; y<L ; y++ ) {
				for (int z=0 ; z<L/2 ; z++ ) {
					int isOddRow = (t & 1) ^ (x & 1) ^ (y & 1);
					for (int mu=0; mu<4 ; mu++ ) {
						int tt = (mu==0)?t+1:t;
						int xx = (mu==1)?x+1:x;
						int yy = (mu==2)?y+1:y;
						int zz = (mu==3)?( (isOddRow)?(z):z+1 ):z;

						tt = (tt+T) % T;
						xx = (xx+L) % L;
						yy = (yy+L) % L;
						zz = (zz+L) % (L/2);

						out[((((((tt*L)+xx)*L)+yy)*L/2)+zz)*4+mu ] = in [ i ];
						i++;
					}
				}
			}
		}
	}
}

/* Converting from qphix style gauge for the whole lattice, to even/odd separated
 * tmLQCD style gauge fields
 */
void devide_gauge_to_oddeven(su3 const * const in, su3 * const even, su3 * const odd) {
	int i = 0;
	for (int t=0 ; t<T ; t++ ) {
		for (int z=0 ; z<L ; z++ ) {
			for (int y=0 ; y<L ; y++ ) {
				for (int x=0 ; x<L/2 ; x++ ) {
					for (int mu=0; mu<4 ; mu++ ) {
						for (int f=-1; f<=1 ; f+= 2) {
							su3 tmp = in[i];

							int mu_ = (mu+1)%4;

							int isOddRow = (t & 1) ^ (z & 1) ^ (y & 1);
							int t_ = t;               // converting from checkerboarded
							int z_ = z/2;             // coordinates of qphix along x-axis
							int y_ = y;               // to tmLQCD checkerboarding along
							int x_ = (2*x)+isOddRow;  // along y-axis

							if (f == -1) {
								int tt = (mu_==0)?t_-1:t_;
								int xx = (mu_==1)?x_-1:x_;
								int yy = (mu_==2)?y_-1:y_;
								int zz = (mu_==3)?( (isOddRow)?(z_-1):z_ ):z_;
								tt = (tt+T) % T;
								xx = (xx+L) % L;
								yy = (yy+L) % L;
								zz = (zz+L) % (L/2);
								even[ ((((((tt*L)+xx)*L)+yy)*L/2)+zz)*4+mu ] = tmp;
							} else {
								int tt = (mu_==0)?t_+1:t_;
								int xx = (mu_==1)?x_+1:x_;
								int yy = (mu_==2)?y_+1:y_;
								int zz = (mu_==3)?( (isOddRow)?(z_):z_+1 ):z_;
								tt = (tt+T) % T;
								xx = (xx+L) % L;
								yy = (yy+L) % L;
								zz = (zz+L) % (L/2);

								//printf("%d\n", ((((((tt*L)+xx)*L)+yy)*L/2)+zz)*4+mu );
								odd[ ((((((tt*L)+xx)*L)+yy)*L/2)+zz)*4+mu ] = tmp;
							}

							i++;
						}
					}
				}
			}
		}
	}
}

/*int LX, LY, LZ;
LX = LY = L+2;
LZ = L/2+1;
int i = 0;
for (int t = 0 ; t < T ; t++ ) {
		for (int x = 0 ; x < L ; x++ ) {
			for (int y = 0 ; y < L ; y++ ) {
				for (int z = 0 ; z < L/2 ; z++ ) {
					int isOddRow = (t & 1) ^ (x & 1) ^ (y & 1) ^ 1;
					int error = compare_spinor(
							out + i  ,
							out_dfe + (t+1)*LX*LY*LZ +
				            (x+1)*LY*LZ +
				            (y+1)*LZ + z + isOddRow
				            );
					i++;
					if (error) {
						printf("Wrong %d! %d\n", error,i);
						spinor *a = out + i;
						spinor *b = out_dfe + i;
						printf("%f %f    %f %f ", creal(a->s0.c0), cimag(a->s0.c0),
								                  creal(b->s0.c0), cimag(b->s0.c0));
						return 1;
					}
				}
			}
		}
	}*



spinor **** create_4d_spinor_wrapper(spinor *s) {
	spinor**** output;
	output = calloc(T, sizeof(spinor***));
	output[0] = calloc(T * LX, sizeof(spinor**));
	for (int i = 1; i < T; i++)
		output[i] = output[0] + i * LX;

	output[0][0] = calloc(T * LX * LY, sizeof(spinor*));
	for (int i = 1; i < T * LX; i++)
		output[0][i] = output[0][0] + i * LY;

	for (int i = 0; i < T * LX * LY; i++)
		output[0][0][i] = s + i * LX / 2;

	return output;
}*/
