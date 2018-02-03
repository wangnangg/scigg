#pragma once
#include <spmatvec.hpp>

namespace scigg
{
// solve general continuous time Markov chain for steady state prob
void general_markov_ss_prob(spmatrix_const_view Q, vector_const_view init_prob,
                            vector_mutable_view ss_prob);
void general_markov_ss_dprob(spmatrix_const_view Q, vector_const_view init_prob,
                             spmatrix_const_view dQ,
                             vector_const_view dinit_prob,
                             vector_const_view ss_prob,
                             vector_mutable_view ss_dprob);
// solve general continuous time Markov chain for steady state cumulative time
void general_markov_ss_cum(spmatrix_const_view Q, vector_const_view init_prob,
                           vector_mutable_view ss_cum);
void general_markov_ss_dcum(spmatrix_const_view Q, vector_const_view init_prob,
                            spmatrix_const_view dQ,
                            vector_const_view dinit_prob,
                            vector_const_view ss_cum,
                            vector_mutable_view ss_dcum);

// transient state...
void general_markov_ts_prob(spmatrix_const_view Q, vector_const_view init_prob,
                            real_t time, vector_mutable_view ss_prob);
void general_markov_ts_cum(spmatrix_const_view Q, vector_const_view init_prob,
                           real_t time, vector_mutable_view ss_cum);
}
