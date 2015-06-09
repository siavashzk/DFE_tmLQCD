#ifndef DFE_HALFSPINOR_HOPPING_H
#define DFE_HALFSPINOR_HOPPING_H

#define _prefetch_spinor(s)
#define _prefetch_halfspinor(hs)
#define _prefetch_su3(U)

#define _hop_t_p_pre()					\
  _vector_assign(rs.s0, s->s0);				\
  _vector_assign(rs.s1, s->s1);				\
  _vector_assign(rs.s2, s->s2);				\
  _vector_assign(rs.s3, s->s3);				\
  _vector_add(psi, rs.s0, rs.s2);			\
  _vector_add(psi2, rs.s1, rs.s3);		\
  _su3_multiply(chi,(*U0),psi);				\
  _su3_multiply(chi2,(*U0),psi2);			\
  _complex_times_vector(phi[ix]->s0, ka0, chi);		\
  _complex_times_vector(phi[ix]->s1, ka0, chi2);

#define _hop_t_m_pre()				\
  _vector_sub(phi[ix]->s0, rs.s0, rs.s2);	\
  _vector_sub(phi[ix]->s1, rs.s1, rs.s3);

#define _hop_x_p_pre()					\
  _vector_i_add(psi, rs.s0, rs.s3);			\
  _vector_i_add(psi2, rs.s1, rs.s2);			\
  _su3_multiply(chi, (*U0), psi);			\
  _su3_multiply(chi2, (*U0), psi2);			\
  _complex_times_vector(phi[ix]->s0, ka1, chi);		\
  _complex_times_vector(phi[ix]->s1, ka1, chi2);

#define _hop_x_m_pre()					\
  _vector_i_sub(phi[ix]->s0, rs.s0, rs.s3);		\
  _vector_i_sub(phi[ix]->s1, rs.s1, rs.s2);

#define _hop_y_p_pre()					\
  _vector_add(psi, rs.s0, rs.s3);			\
  _vector_sub(psi2, rs.s1, rs.s2);			\
  _su3_multiply(chi,(*U0),psi);				\
  _su3_multiply(chi2,(*U0),psi2);			\
  _complex_times_vector(phi[ix]->s0, ka2, chi);		\
  _complex_times_vector(phi[ix]->s1, ka2, chi2);

#define _hop_y_m_pre()					\
  _vector_sub(phi[ix]->s0, rs.s0, rs.s3);		\
  _vector_add(phi[ix]->s1, rs.s1, rs.s2);

#define _hop_z_p_pre()					\
  _vector_i_add(psi, rs.s0, rs.s2);			\
  _vector_i_sub(psi2, rs.s1, rs.s3);			\
  _su3_multiply(chi, (*U0), psi);			\
  _su3_multiply(chi2,(*U0),psi2);			\
  _complex_times_vector(phi[ix]->s0, ka3, chi);		\
  _complex_times_vector(phi[ix]->s1, ka3, chi2);

#define _hop_z_m_pre()					\
  _vector_i_sub(phi[ix]->s0, rs.s0, rs.s2);		\
  _vector_i_add(phi[ix]->s1, rs.s1, rs.s3);

#define _hop_t_p_post()					\
  _vector_assign(rs.s0, phi[ix]->s0);			\
  _vector_assign(rs.s2, phi[ix]->s0);			\
  _vector_assign(rs.s1, phi[ix]->s1);			\
  _vector_assign(rs.s3, phi[ix]->s1);

#define _hop_t_m_post()					\
  _su3_inverse_multiply(chi,(*U1),phi[ix]->s0);		\
  _su3_inverse_multiply(chi2,(*U1),phi[ix]->s1);		\
  _complexcjg_times_vector(psi,ka0,chi);		\
  _complexcjg_times_vector(psi2,ka0,chi2);		\
  _vector_add_assign(rs.s0, psi);			\
  _vector_sub_assign(rs.s2, psi);			\
  _vector_add_assign(rs.s1, psi2);			\
  _vector_sub_assign(rs.s3, psi2);

#define _hop_x_p_post()					\
  _vector_add_assign(rs.s0, phi[ix]->s0);		\
  _vector_i_sub_assign(rs.s3, phi[ix]->s0);		\
  _vector_add_assign(rs.s1, phi[ix]->s1);		\
  _vector_i_sub_assign(rs.s2, phi[ix]->s1);

#define _hop_x_m_post()					\
  _su3_inverse_multiply(chi,(*U1), phi[ix]->s0);		\
  _su3_inverse_multiply(chi2, (*U1), phi[ix]->s1);	\
  _complexcjg_times_vector(psi,ka1,chi);		\
  _complexcjg_times_vector(psi2,ka1,chi2);		\
  _vector_add_assign(rs.s0, psi);			\
  _vector_i_add_assign(rs.s3, psi);			\
  _vector_add_assign(rs.s1, psi2);			\
  _vector_i_add_assign(rs.s2, psi2);

#define _hop_y_p_post()					\
  _vector_add_assign(rs.s0, phi[ix]->s0);		\
  _vector_add_assign(rs.s3, phi[ix]->s0);		\
  _vector_add_assign(rs.s1, phi[ix]->s1);		\
  _vector_sub_assign(rs.s2, phi[ix]->s1);

#define _hop_y_m_post()					\
  _su3_inverse_multiply(chi,(*U1), phi[ix]->s0);		\
  _su3_inverse_multiply(chi2, (*U1), phi[ix]->s1);	\
  _complexcjg_times_vector(psi,ka2,chi);		\
  _complexcjg_times_vector(psi2,ka2,chi2);		\
  _vector_add_assign(rs.s0, psi);			\
  _vector_sub_assign(rs.s3, psi);			\
  _vector_add_assign(rs.s1, psi2);			\
  _vector_add_assign(rs.s2, psi2);

#define _hop_z_p_post()					\
  _vector_add_assign(rs.s0, phi[ix]->s0);		\
  _vector_i_sub_assign(rs.s2, phi[ix]->s0);		\
  _vector_add_assign(rs.s1, phi[ix]->s1);		\
  _vector_i_add_assign(rs.s3, phi[ix]->s1);

#define _hop_z_m_post()					\
  _su3_inverse_multiply(chi,(*U1), phi[ix]->s0);		\
  _su3_inverse_multiply(chi2, (*U1), phi[ix]->s1);	\
  _complexcjg_times_vector(psi,ka3,chi);		\
  _complexcjg_times_vector(psi2,ka3,chi2);		\
  _vector_add_assign(rs.s0, psi);			\
  _vector_add_assign(rs.s1, psi2);			\
  _vector_i_add_assign(rs.s2, psi);			\
  _vector_i_sub_assign(rs.s3, psi2);

#define _hop_mul_g5_cmplx_and_store(res)			\
  _complex_times_vector((res)->s0, cfactor, rs.s0);		\
  _complex_times_vector((res)->s1, cfactor, rs.s1);		\
  _complexcjg_times_vector((res)->s2, cfactor, rs.s2);	\
  _complexcjg_times_vector((res)->s3, cfactor, rs.s3);

#define _g5_cmplx_sub_hop_and_g5store(res)		\
  _complex_times_vector(psi, cfactor, pn->s0);		\
  _vector_sub((res)->s0, psi, rs.s0);			\
  _complex_times_vector(psi2, cfactor, pn->s1);		\
  _vector_sub((res)->s1, psi2, rs.s1);			\
  _complexcjg_times_vector(psi, cfactor, pn->s2);	\
  _vector_sub((res)->s2, rs.s2, psi);			\
  _complexcjg_times_vector(psi2, cfactor, pn->s3);	\
  _vector_sub((res)->s3, rs.s3, psi2);


#define _hop_store_post(res)		\
  _vector_assign(res->s0, rs.s0);	\
  _vector_assign(res->s1, rs.s1);	\
  _vector_assign(res->s2, rs.s2);	\
  _vector_assign(res->s3, rs.s3);


#define _declare_hregs()			\
  spinor rs;					\
  su3_vector psi, chi, psi2, chi2;

#endif

