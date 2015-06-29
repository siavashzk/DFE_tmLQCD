
int ix;
spinor * restrict s;
halfspinor * restrict * phi;
spinor rs;
su3_vector psi, chi, psi2, chi2;


s = k;
_prefetch_spinor(s);
_prefetch_su3(U0);

FILE *fptr = fopen ("debug_output.txt","w");

   phi = NBPointer[ieo];
   ix=0;
   for(int i = 0; i < (VOLUME)/2; i++){

	 //printf("cpu %f %f\n", creal(s->s0.c0), cimag(s->s0.c0) );
     _hop_t_p_pre();
     //printf("cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     U0++;
     ix++;
     
     _hop_t_m_pre();
     //printf("cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     ix++;
     
     _hop_x_p_pre();
     //printf("cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     U0++;
     ix++;
     
     _hop_x_m_pre();
     //printf("cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     ix++;
     
     _hop_y_p_pre();
     //printf("cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     U0++;
     ix++;
     
     _hop_y_m_pre();
     //printf("cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     ix++;
     
     _hop_z_p_pre();
     //printf("cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     U0++;
     ix++;
     
     _hop_z_m_pre();
     //printf("cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );

     s++;            
     ix++;
     
   }
   
   s = l;
   _prefetch_su3(U1);
   
   phi = NBPointer[ieo+2];
   
   ix = 0;
   for(int i = 0; i < (VOLUME)/2; i++){

#ifdef _TM_SUB_HOP
     pn=p+i;
#endif

     _hop_t_p_post();
     //fprintf(fptr,"cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     ix++;
     
     _hop_t_m_post();
     //fprintf(fptr,"cpu %f %f\n", creal(psi.c0), cimag(psi.c0) );
     ix++;
     U1++;
     
     _hop_x_p_post();
     //fprintf(fptr,"cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     ix++;
     
     _hop_x_m_post();
     //fprintf(fptr,"cpu %f %f\n", creal(psi.c0), cimag(psi.c0) );
     U1++;
     ix++;
     
     _hop_y_p_post();
     //fprintf(fptr,"cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     ix++;
     
     _hop_y_m_post();
     //fprintf(fptr,"cpu %f %f\n", creal(psi.c0), cimag(psi.c0) );
     U1++;
     ix++;
     
     _hop_z_p_post();
     //fprintf(fptr,"cpu %f %f\n", creal(phi[ix]->s0.c0), cimag(phi[ix]->s0.c0) );
     ix++;
     
     _hop_z_m_post();
     //fprintf(fptr,"cpu %f %f\n", creal(psi.c0), cimag(psi.c0) );
     
#ifdef _MUL_G5_CMPLX
     _hop_mul_g5_cmplx_and_store(s);
#elif defined _TM_SUB_HOP
     //printf("cpu %f %f\n", creal(pn->s0.c0), cimag(pn->s0.c0) );
     _g5_cmplx_sub_hop_and_g5store(s);
#else
     _hop_store_post(s);
#endif

     //printf("cpu %f %f\n", creal(s->s0.c0), cimag(s->s0.c0) );
     U1++;
     ix++;
     s++;
   }

