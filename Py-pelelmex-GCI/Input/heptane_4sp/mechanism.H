#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  C
1  H
2  O
3  N
*/

// Species
#define NC7H16_ID 0
#define IC8H18_ID 1
#define O2_ID 2
#define N2_ID 3

#define NUM_ELEMENTS 4
#define NUM_SPECIES 4
#define NUM_IONS 0
#define NUM_REACTIONS 2

#define NUM_FIT 4

//  ALWAYS on CPU stuff -- can have different def depending on if we are CPU or
//  GPU based. Defined in mechanism.cpp
void atomicWeight(amrex::Real* awt);
//  MISC
void CKAWT(amrex::Real* awt);
void CKNCF(int* ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void GET_RMAP(int* _rmap);
void CKINU(const int i, int& nspec, int* ki, int* nu);
void CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real* x,
  amrex::Real* q_f,
  amrex::Real* q_r);
void progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T);
//  SPARSE INFORMATION
void SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS);
void SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP);
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS);
void SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base);

// A few mechanism parameters
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKINDX(int& mm, int& kk, int& ii, int& nfit)
{
  mm = 4;
  kk = 4;
  ii = 2;
  nfit = -1; // Why do you need this anyway ?
}

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[4] = {
  0.0099795419390250, // NC7H16
  0.0087541144337839, // IC8H18
  0.0312519532470779, // O2
  0.0356964374955379, // N2
};
#endif
const amrex::Real h_global_imw[4] = {
  0.0099795419390250, // NC7H16
  0.0087541144337839, // IC8H18
  0.0312519532470779, // O2
  0.0356964374955379, // N2
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[4] = {
  100.205000, // NC7H16
  114.232000, // IC8H18
  31.998000,  // O2
  28.014000,  // N2
};
#endif
const amrex::Real h_global_mw[4] = {
  100.205000, // NC7H16
  114.232000, // IC8H18
  31.998000,  // O2
  28.014000,  // N2
};

//  inverse molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real* imw_new)
{
  imw_new[0] = 0.0099795419390250; // NC7H16
  imw_new[1] = 0.0087541144337839; // IC8H18
  imw_new[2] = 0.0312519532470779; // O2
  imw_new[3] = 0.0356964374955379; // N2
}

//  inverse molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
imw(const int n)
{
#if AMREX_DEVICE_COMPILE
  return global_imw[n];
#else
  return h_global_imw[n];
#endif
}
//  molecular weights
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_mw(amrex::Real* mw_new)
{
  mw_new[0] = 100.205000; // NC7H16
  mw_new[1] = 114.232000; // IC8H18
  mw_new[2] = 31.998000;  // O2
  mw_new[3] = 28.014000;  // N2
}

//  molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
mw(const int n)
{
#if AMREX_DEVICE_COMPILE
  return global_mw[n];
#else
  return h_global_mw[n];
#endif
}

//  Returns R, Rc, Patm
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRP(amrex::Real& ru, amrex::Real& ruc, amrex::Real& pa)
{
  ru = 8.31446261815324e+07;
  ruc = 1.98721558317399615845;
  pa = 1.01325e+06;
}

// compute Cv/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cv_R(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    species[2] = +2.21293640e+00 + 1.12748635e-03 * T - 5.75615047e-07 * T2 +
                 1.31387723e-09 * T3 - 8.76855392e-13 * T4;
    // species 3: N2
    species[3] = +2.29867700e+00 + 1.40824000e-03 * T - 3.96322200e-06 * T2 +
                 5.64151500e-09 * T3 - 2.44485500e-12 * T4;
  } else {
    // species 2: O2
    species[2] = +2.69757819e+00 + 6.13519689e-04 * T - 1.25884199e-07 * T2 +
                 1.77528148e-11 * T3 - 1.13643531e-15 * T4;
    // species 3: N2
    species[3] = +1.92664000e+00 + 1.48797700e-03 * T - 5.68476100e-07 * T2 +
                 1.00970400e-10 * T3 - 6.75335100e-15 * T4;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    species[0] = -2.26836187e+00 + 8.54355820e-02 * T - 5.25346786e-05 * T2 +
                 1.62945721e-08 * T3 - 2.02394925e-12 * T4;
  } else {
    // species 0: NC7H16
    species[0] = +2.12148969e+01 + 3.47675750e-02 * T - 1.18407129e-05 * T2 +
                 1.83298478e-09 * T3 - 1.06130266e-13 * T4;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    species[1] = -5.20868893e+00 + 1.11440581e-01 * T - 7.91346582e-05 * T2 +
                 2.92406242e-08 * T3 - 4.43743191e-12 * T4;
  } else {
    // species 1: IC8H18
    species[1] = +2.61373590e+01 + 3.79004890e-02 * T - 1.29437358e-05 * T2 +
                 2.00760372e-09 * T3 - 1.16400580e-13 * T4;
  }
}

// compute Cp/R at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cp_R(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    species[2] = +3.21293640e+00 + 1.12748635e-03 * T - 5.75615047e-07 * T2 +
                 1.31387723e-09 * T3 - 8.76855392e-13 * T4;
    // species 3: N2
    species[3] = +3.29867700e+00 + 1.40824000e-03 * T - 3.96322200e-06 * T2 +
                 5.64151500e-09 * T3 - 2.44485500e-12 * T4;
  } else {
    // species 2: O2
    species[2] = +3.69757819e+00 + 6.13519689e-04 * T - 1.25884199e-07 * T2 +
                 1.77528148e-11 * T3 - 1.13643531e-15 * T4;
    // species 3: N2
    species[3] = +2.92664000e+00 + 1.48797700e-03 * T - 5.68476100e-07 * T2 +
                 1.00970400e-10 * T3 - 6.75335100e-15 * T4;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    species[0] = -1.26836187e+00 + 8.54355820e-02 * T - 5.25346786e-05 * T2 +
                 1.62945721e-08 * T3 - 2.02394925e-12 * T4;
  } else {
    // species 0: NC7H16
    species[0] = +2.22148969e+01 + 3.47675750e-02 * T - 1.18407129e-05 * T2 +
                 1.83298478e-09 * T3 - 1.06130266e-13 * T4;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    species[1] = -4.20868893e+00 + 1.11440581e-01 * T - 7.91346582e-05 * T2 +
                 2.92406242e-08 * T3 - 4.43743191e-12 * T4;
  } else {
    // species 1: IC8H18
    species[1] = +2.71373590e+01 + 3.79004890e-02 * T - 1.29437358e-05 * T2 +
                 2.00760372e-09 * T3 - 1.16400580e-13 * T4;
  }
}

// compute the g/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
gibbs(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    species[2] = -1.005249020000000e+03 * invT - 2.821801190000000e+00 -
                 3.212936400000000e+00 * logT - 5.637431750000000e-04 * T +
                 9.593584116666666e-08 * T2 - 1.094897691666667e-10 * T3 +
                 4.384276960000000e-14 * T4;
    // species 3: N2
    species[3] = -1.020900000000000e+03 * invT - 6.516950000000001e-01 -
                 3.298677000000000e+00 * logT - 7.041200000000000e-04 * T +
                 6.605369999999999e-07 * T2 - 4.701262500000001e-10 * T3 +
                 1.222427500000000e-13 * T4;
  } else {
    // species 2: O2
    species[2] = -1.233930180000000e+03 * invT + 5.084126000000002e-01 -
                 3.697578190000000e+00 * logT - 3.067598445000000e-04 * T +
                 2.098069983333333e-08 * T2 - 1.479401233333333e-12 * T3 +
                 5.682176550000000e-17 * T4;
    // species 3: N2
    species[3] = -9.227977000000000e+02 * invT - 3.053888000000000e+00 -
                 2.926640000000000e+00 * logT - 7.439885000000000e-04 * T +
                 9.474601666666666e-08 * T2 - 8.414199999999999e-12 * T3 +
                 3.376675500000000e-16 * T4;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    species[0] = -2.565865650000000e+04 * invT - 3.664165307000000e+01 +
                 1.268361870000000e+00 * logT - 4.271779100000000e-02 * T +
                 8.755779766666667e-06 * T2 - 1.357881008333333e-09 * T3 +
                 1.011974625000000e-13 * T4;
  } else {
    // species 0: NC7H16
    species[0] = -3.427600810000000e+04 * invT + 1.145189165000000e+02 -
                 2.221489690000000e+01 * logT - 1.738378750000000e-02 * T +
                 1.973452150000000e-06 * T2 - 1.527487316666667e-10 * T3 +
                 5.306513300000000e-15 * T4;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    species[1] = -2.994468750000000e+04 * invT - 4.916085903000000e+01 +
                 4.208688930000000e+00 * logT - 5.572029050000000e-02 * T +
                 1.318910970000000e-05 * T2 - 2.436718683333333e-09 * T3 +
                 2.218715955000000e-13 * T4;
  } else {
    // species 1: IC8H18
    species[1] = -4.079581770000000e+04 * invT + 1.504148540000000e+02 -
                 2.713735900000000e+01 * logT - 1.895024450000000e-02 * T +
                 2.157289300000000e-06 * T2 - 1.673003100000000e-10 * T3 +
                 5.820029000000000e-15 * T4;
  }
}

// compute the a/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
helmholtz(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    species[2] = -1.00524902e+03 * invT - 3.82180119e+00 -
                 3.21293640e+00 * logT - 5.63743175e-04 * T +
                 9.59358412e-08 * T2 - 1.09489769e-10 * T3 +
                 4.38427696e-14 * T4;
    // species 3: N2
    species[3] = -1.02090000e+03 * invT - 1.65169500e+00 -
                 3.29867700e+00 * logT - 7.04120000e-04 * T +
                 6.60537000e-07 * T2 - 4.70126250e-10 * T3 +
                 1.22242750e-13 * T4;
  } else {
    // species 2: O2
    species[2] = -1.23393018e+03 * invT - 4.91587400e-01 -
                 3.69757819e+00 * logT - 3.06759845e-04 * T +
                 2.09806998e-08 * T2 - 1.47940123e-12 * T3 +
                 5.68217655e-17 * T4;
    // species 3: N2
    species[3] = -9.22797700e+02 * invT - 4.05388800e+00 -
                 2.92664000e+00 * logT - 7.43988500e-04 * T +
                 9.47460167e-08 * T2 - 8.41420000e-12 * T3 +
                 3.37667550e-16 * T4;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    species[0] = -2.56586565e+04 * invT - 3.76416531e+01 +
                 1.26836187e+00 * logT - 4.27177910e-02 * T +
                 8.75577977e-06 * T2 - 1.35788101e-09 * T3 +
                 1.01197462e-13 * T4;
  } else {
    // species 0: NC7H16
    species[0] = -3.42760081e+04 * invT + 1.13518917e+02 -
                 2.22148969e+01 * logT - 1.73837875e-02 * T +
                 1.97345215e-06 * T2 - 1.52748732e-10 * T3 +
                 5.30651330e-15 * T4;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    species[1] = -2.99446875e+04 * invT - 5.01608590e+01 +
                 4.20868893e+00 * logT - 5.57202905e-02 * T +
                 1.31891097e-05 * T2 - 2.43671868e-09 * T3 +
                 2.21871595e-13 * T4;
  } else {
    // species 1: IC8H18
    species[1] = -4.07958177e+04 * invT + 1.49414854e+02 -
                 2.71373590e+01 * logT - 1.89502445e-02 * T +
                 2.15728930e-06 * T2 - 1.67300310e-10 * T3 +
                 5.82002900e-15 * T4;
  }
}

// compute the e/(RT) at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesInternalEnergy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    species[2] = +2.21293640e+00 + 5.63743175e-04 * T - 1.91871682e-07 * T2 +
                 3.28469308e-10 * T3 - 1.75371078e-13 * T4 -
                 1.00524902e+03 * invT;
    // species 3: N2
    species[3] = +2.29867700e+00 + 7.04120000e-04 * T - 1.32107400e-06 * T2 +
                 1.41037875e-09 * T3 - 4.88971000e-13 * T4 -
                 1.02090000e+03 * invT;
  } else {
    // species 2: O2
    species[2] = +2.69757819e+00 + 3.06759845e-04 * T - 4.19613997e-08 * T2 +
                 4.43820370e-12 * T3 - 2.27287062e-16 * T4 -
                 1.23393018e+03 * invT;
    // species 3: N2
    species[3] = +1.92664000e+00 + 7.43988500e-04 * T - 1.89492033e-07 * T2 +
                 2.52426000e-11 * T3 - 1.35067020e-15 * T4 -
                 9.22797700e+02 * invT;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    species[0] = -2.26836187e+00 + 4.27177910e-02 * T - 1.75115595e-05 * T2 +
                 4.07364302e-09 * T3 - 4.04789850e-13 * T4 -
                 2.56586565e+04 * invT;
  } else {
    // species 0: NC7H16
    species[0] = +2.12148969e+01 + 1.73837875e-02 * T - 3.94690430e-06 * T2 +
                 4.58246195e-10 * T3 - 2.12260532e-14 * T4 -
                 3.42760081e+04 * invT;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    species[1] = -5.20868893e+00 + 5.57202905e-02 * T - 2.63782194e-05 * T2 +
                 7.31015605e-09 * T3 - 8.87486382e-13 * T4 -
                 2.99446875e+04 * invT;
  } else {
    // species 1: IC8H18
    species[1] = +2.61373590e+01 + 1.89502445e-02 * T - 4.31457860e-06 * T2 +
                 5.01900930e-10 * T3 - 2.32801160e-14 * T4 -
                 4.07958177e+04 * invT;
  }
}

// compute the h/(RT) at the given temperature (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEnthalpy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    species[2] = +3.21293640e+00 + 5.63743175e-04 * T - 1.91871682e-07 * T2 +
                 3.28469308e-10 * T3 - 1.75371078e-13 * T4 -
                 1.00524902e+03 * invT;
    // species 3: N2
    species[3] = +3.29867700e+00 + 7.04120000e-04 * T - 1.32107400e-06 * T2 +
                 1.41037875e-09 * T3 - 4.88971000e-13 * T4 -
                 1.02090000e+03 * invT;
  } else {
    // species 2: O2
    species[2] = +3.69757819e+00 + 3.06759845e-04 * T - 4.19613997e-08 * T2 +
                 4.43820370e-12 * T3 - 2.27287062e-16 * T4 -
                 1.23393018e+03 * invT;
    // species 3: N2
    species[3] = +2.92664000e+00 + 7.43988500e-04 * T - 1.89492033e-07 * T2 +
                 2.52426000e-11 * T3 - 1.35067020e-15 * T4 -
                 9.22797700e+02 * invT;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    species[0] = -1.26836187e+00 + 4.27177910e-02 * T - 1.75115595e-05 * T2 +
                 4.07364302e-09 * T3 - 4.04789850e-13 * T4 -
                 2.56586565e+04 * invT;
  } else {
    // species 0: NC7H16
    species[0] = +2.22148969e+01 + 1.73837875e-02 * T - 3.94690430e-06 * T2 +
                 4.58246195e-10 * T3 - 2.12260532e-14 * T4 -
                 3.42760081e+04 * invT;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    species[1] = -4.20868893e+00 + 5.57202905e-02 * T - 2.63782194e-05 * T2 +
                 7.31015605e-09 * T3 - 8.87486382e-13 * T4 -
                 2.99446875e+04 * invT;
  } else {
    // species 1: IC8H18
    species[1] = +2.71373590e+01 + 1.89502445e-02 * T - 4.31457860e-06 * T2 +
                 5.01900930e-10 * T3 - 2.32801160e-14 * T4 -
                 4.07958177e+04 * invT;
  }
}

// compute the S/R at the given temperature (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
speciesEntropy(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real logT = log(T);

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    species[2] = +3.21293640e+00 * logT + 1.12748635e-03 * T -
                 2.87807523e-07 * T2 + 4.37959077e-10 * T3 -
                 2.19213848e-13 * T4 + 6.03473759e+00;
    // species 3: N2
    species[3] = +3.29867700e+00 * logT + 1.40824000e-03 * T -
                 1.98161100e-06 * T2 + 1.88050500e-09 * T3 -
                 6.11213750e-13 * T4 + 3.95037200e+00;
  } else {
    // species 2: O2
    species[2] = +3.69757819e+00 * logT + 6.13519689e-04 * T -
                 6.29420995e-08 * T2 + 5.91760493e-12 * T3 -
                 2.84108828e-16 * T4 + 3.18916559e+00;
    // species 3: N2
    species[3] = +2.92664000e+00 * logT + 1.48797700e-03 * T -
                 2.84238050e-07 * T2 + 3.36568000e-11 * T3 -
                 1.68833775e-15 * T4 + 5.98052800e+00;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    species[0] = -1.26836187e+00 * logT + 8.54355820e-02 * T -
                 2.62673393e-05 * T2 + 5.43152403e-09 * T3 -
                 5.05987313e-13 * T4 + 3.53732912e+01;
  } else {
    // species 0: NC7H16
    species[0] = +2.22148969e+01 * logT + 3.47675750e-02 * T -
                 5.92035645e-06 * T2 + 6.10994927e-10 * T3 -
                 2.65325665e-14 * T4 - 9.23040196e+01;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    species[1] = -4.20868893e+00 * logT + 1.11440581e-01 * T -
                 3.95673291e-05 * T2 + 9.74687473e-09 * T3 -
                 1.10935798e-12 * T4 + 4.49521701e+01;
  } else {
    // species 1: IC8H18
    species[1] = +2.71373590e+01 * logT + 3.79004890e-02 * T -
                 6.47186790e-06 * T2 + 6.69201240e-10 * T3 -
                 2.91001450e-14 * T4 - 1.23277495e+02;
  }
}

// compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
dcvpRdT(amrex::Real* species, const amrex::Real T)
{
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    species[2] = +1.12748635e-03 - 1.15123009e-06 * T + 3.94163169e-09 * T2 -
                 3.50742157e-12 * T3;
    // species 3: N2
    species[3] = +1.40824000e-03 - 7.92644400e-06 * T + 1.69245450e-08 * T2 -
                 9.77942000e-12 * T3;
  } else {
    // species 2: O2
    species[2] = +6.13519689e-04 - 2.51768398e-07 * T + 5.32584444e-11 * T2 -
                 4.54574124e-15 * T3;
    // species 3: N2
    species[3] = +1.48797700e-03 - 1.13695220e-06 * T + 3.02911200e-10 * T2 -
                 2.70134040e-14 * T3;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    species[0] = +8.54355820e-02 - 1.05069357e-04 * T + 4.88837163e-08 * T2 -
                 8.09579700e-12 * T3;
  } else {
    // species 0: NC7H16
    species[0] = +3.47675750e-02 - 2.36814258e-05 * T + 5.49895434e-09 * T2 -
                 4.24521064e-13 * T3;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    species[1] = +1.11440581e-01 - 1.58269316e-04 * T + 8.77218726e-08 * T2 -
                 1.77497276e-11 * T3;
  } else {
    // species 1: IC8H18
    species[1] = +3.79004890e-02 - 2.58874716e-05 * T + 6.02281116e-09 * T2 -
                 4.65602320e-13 * T3;
  }
}

// Returns the mean specific heat at CP (Eq. 33)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cpbl)
{
  amrex::Real result = 0;
  amrex::Real cpor[4]; //  temporary storage
  cp_R(cpor, T);

  // perform dot product
  for (int id = 0; id < 4; ++id) {
    result += x[id] * cpor[id];
  }

  cpbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CP (Eq. 34)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cpbs)
{
  amrex::Real result = 0.0;

  // compute Cp/R at the given temperature
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    result += y[2] *
              (+3.21293640e+00 + 1.12748635e-03 * T - 5.75615047e-07 * T2 +
               1.31387723e-09 * T3 - 8.76855392e-13 * T4) *
              0.0312519532470779;
    // species 3: N2
    result += y[3] *
              (+3.29867700e+00 + 1.40824000e-03 * T - 3.96322200e-06 * T2 +
               5.64151500e-09 * T3 - 2.44485500e-12 * T4) *
              0.0356964374955379;
  } else {
    // species 2: O2
    result += y[2] *
              (+3.69757819e+00 + 6.13519689e-04 * T - 1.25884199e-07 * T2 +
               1.77528148e-11 * T3 - 1.13643531e-15 * T4) *
              0.0312519532470779;
    // species 3: N2
    result += y[3] *
              (+2.92664000e+00 + 1.48797700e-03 * T - 5.68476100e-07 * T2 +
               1.00970400e-10 * T3 - 6.75335100e-15 * T4) *
              0.0356964374955379;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    result += y[0] *
              (-1.26836187e+00 + 8.54355820e-02 * T - 5.25346786e-05 * T2 +
               1.62945721e-08 * T3 - 2.02394925e-12 * T4) *
              0.0099795419390250;
  } else {
    // species 0: NC7H16
    result += y[0] *
              (+2.22148969e+01 + 3.47675750e-02 * T - 1.18407129e-05 * T2 +
               1.83298478e-09 * T3 - 1.06130266e-13 * T4) *
              0.0099795419390250;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    result += y[1] *
              (-4.20868893e+00 + 1.11440581e-01 * T - 7.91346582e-05 * T2 +
               2.92406242e-08 * T3 - 4.43743191e-12 * T4) *
              0.0087541144337839;
  } else {
    // species 1: IC8H18
    result += y[1] *
              (+2.71373590e+01 + 3.79004890e-02 * T - 1.29437358e-05 * T2 +
               2.00760372e-09 * T3 - 1.16400580e-13 * T4) *
              0.0087541144337839;
  }

  cpbs = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 35)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBL(const amrex::Real T, const amrex::Real x[], amrex::Real& cvbl)
{
  amrex::Real result = 0;
  amrex::Real cvor[4]; //  temporary storage
  cv_R(cvor, T);

  // perform dot product
  for (int id = 0; id < 4; ++id) {
    result += x[id] * cvor[id];
  }

  cvbl = result * 8.31446261815324e+07;
}

// Returns the mean specific heat at CV (Eq. 36)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBS(const amrex::Real T, const amrex::Real y[], amrex::Real& cvbs)
{
  amrex::Real result = 0.0;
  // compute Cv/R at the given temperature
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    result += y[2] *
              (+2.21293640e+00 + 1.12748635e-03 * T - 5.75615047e-07 * T2 +
               1.31387723e-09 * T3 - 8.76855392e-13 * T4) *
              0.0312519532470779;
    // species 3: N2
    result += y[3] *
              (+2.29867700e+00 + 1.40824000e-03 * T - 3.96322200e-06 * T2 +
               5.64151500e-09 * T3 - 2.44485500e-12 * T4) *
              0.0356964374955379;
  } else {
    // species 2: O2
    result += y[2] *
              (+2.69757819e+00 + 6.13519689e-04 * T - 1.25884199e-07 * T2 +
               1.77528148e-11 * T3 - 1.13643531e-15 * T4) *
              0.0312519532470779;
    // species 3: N2
    result += y[3] *
              (+1.92664000e+00 + 1.48797700e-03 * T - 5.68476100e-07 * T2 +
               1.00970400e-10 * T3 - 6.75335100e-15 * T4) *
              0.0356964374955379;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    result += y[0] *
              (-2.26836187e+00 + 8.54355820e-02 * T - 5.25346786e-05 * T2 +
               1.62945721e-08 * T3 - 2.02394925e-12 * T4) *
              0.0099795419390250;
  } else {
    // species 0: NC7H16
    result += y[0] *
              (+2.12148969e+01 + 3.47675750e-02 * T - 1.18407129e-05 * T2 +
               1.83298478e-09 * T3 - 1.06130266e-13 * T4) *
              0.0099795419390250;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    result += y[1] *
              (-5.20868893e+00 + 1.11440581e-01 * T - 7.91346582e-05 * T2 +
               2.92406242e-08 * T3 - 4.43743191e-12 * T4) *
              0.0087541144337839;
  } else {
    // species 1: IC8H18
    result += y[1] *
              (+2.61373590e+01 + 3.79004890e-02 * T - 1.29437358e-05 * T2 +
               2.00760372e-09 * T3 - 1.16400580e-13 * T4) *
              0.0087541144337839;
  }

  cvbs = result * 8.31446261815324e+07;
}

// Returns the mean enthalpy of the mixture in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBML(const amrex::Real T, const amrex::Real x[], amrex::Real& hbml)
{
  amrex::Real result = 0;
  amrex::Real hml[4];                        //  temporary storage
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesEnthalpy(hml, T);

  // perform dot product
  for (int id = 0; id < 4; ++id) {
    result += x[id] * hml[id];
  }

  hbml = result * RT;
}

// Returns mean enthalpy of mixture in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& hbms)
{
  amrex::Real result = 0.0;
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    result +=
      y[2] *
      (+3.21293640e+00 + 5.63743175e-04 * T - 1.91871682e-07 * T2 +
       3.28469308e-10 * T3 - 1.75371078e-13 * T4 - 1.00524902e+03 * invT) *
      0.0312519532470779;
    // species 3: N2
    result +=
      y[3] *
      (+3.29867700e+00 + 7.04120000e-04 * T - 1.32107400e-06 * T2 +
       1.41037875e-09 * T3 - 4.88971000e-13 * T4 - 1.02090000e+03 * invT) *
      0.0356964374955379;
  } else {
    // species 2: O2
    result +=
      y[2] *
      (+3.69757819e+00 + 3.06759845e-04 * T - 4.19613997e-08 * T2 +
       4.43820370e-12 * T3 - 2.27287062e-16 * T4 - 1.23393018e+03 * invT) *
      0.0312519532470779;
    // species 3: N2
    result +=
      y[3] *
      (+2.92664000e+00 + 7.43988500e-04 * T - 1.89492033e-07 * T2 +
       2.52426000e-11 * T3 - 1.35067020e-15 * T4 - 9.22797700e+02 * invT) *
      0.0356964374955379;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    result +=
      y[0] *
      (-1.26836187e+00 + 4.27177910e-02 * T - 1.75115595e-05 * T2 +
       4.07364302e-09 * T3 - 4.04789850e-13 * T4 - 2.56586565e+04 * invT) *
      0.0099795419390250;
  } else {
    // species 0: NC7H16
    result +=
      y[0] *
      (+2.22148969e+01 + 1.73837875e-02 * T - 3.94690430e-06 * T2 +
       4.58246195e-10 * T3 - 2.12260532e-14 * T4 - 3.42760081e+04 * invT) *
      0.0099795419390250;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    result +=
      y[1] *
      (-4.20868893e+00 + 5.57202905e-02 * T - 2.63782194e-05 * T2 +
       7.31015605e-09 * T3 - 8.87486382e-13 * T4 - 2.99446875e+04 * invT) *
      0.0087541144337839;
  } else {
    // species 1: IC8H18
    result +=
      y[1] *
      (+2.71373590e+01 + 1.89502445e-02 * T - 4.31457860e-06 * T2 +
       5.01900930e-10 * T3 - 2.32801160e-14 * T4 - 4.07958177e+04 * invT) *
      0.0087541144337839;
  }

  const amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  hbms = result * RT;
}

// get mean internal energy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBML(const amrex::Real T, const amrex::Real x[], amrex::Real& ubml)
{
  amrex::Real result = 0;
  amrex::Real uml[4];                        //  temporary energy array
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesInternalEnergy(uml, T);

  // perform dot product
  for (int id = 0; id < 4; ++id) {
    result += x[id] * uml[id];
  }

  ubml = result * RT;
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBMS(const amrex::Real T, const amrex::Real y[], amrex::Real& ubms)
{
  amrex::Real result = 0.0;
  const amrex::Real T2 = T * T;
  const amrex::Real T3 = T * T * T;
  const amrex::Real T4 = T * T * T * T;
  const amrex::Real invT = 1.0 / T;

  // species with midpoint at T=1000 kelvin
  if (T < 1000) {
    // species 2: O2
    result +=
      y[2] *
      (+2.21293640e+00 + 5.63743175e-04 * T - 1.91871682e-07 * T2 +
       3.28469308e-10 * T3 - 1.75371078e-13 * T4 - 1.00524902e+03 * invT) *
      0.0312519532470779;
    // species 3: N2
    result +=
      y[3] *
      (+2.29867700e+00 + 7.04120000e-04 * T - 1.32107400e-06 * T2 +
       1.41037875e-09 * T3 - 4.88971000e-13 * T4 - 1.02090000e+03 * invT) *
      0.0356964374955379;
  } else {
    // species 2: O2
    result +=
      y[2] *
      (+2.69757819e+00 + 3.06759845e-04 * T - 4.19613997e-08 * T2 +
       4.43820370e-12 * T3 - 2.27287062e-16 * T4 - 1.23393018e+03 * invT) *
      0.0312519532470779;
    // species 3: N2
    result +=
      y[3] *
      (+1.92664000e+00 + 7.43988500e-04 * T - 1.89492033e-07 * T2 +
       2.52426000e-11 * T3 - 1.35067020e-15 * T4 - 9.22797700e+02 * invT) *
      0.0356964374955379;
  }

  // species with midpoint at T=1391 kelvin
  if (T < 1391) {
    // species 0: NC7H16
    result +=
      y[0] *
      (-2.26836187e+00 + 4.27177910e-02 * T - 1.75115595e-05 * T2 +
       4.07364302e-09 * T3 - 4.04789850e-13 * T4 - 2.56586565e+04 * invT) *
      0.0099795419390250;
  } else {
    // species 0: NC7H16
    result +=
      y[0] *
      (+2.12148969e+01 + 1.73837875e-02 * T - 3.94690430e-06 * T2 +
       4.58246195e-10 * T3 - 2.12260532e-14 * T4 - 3.42760081e+04 * invT) *
      0.0099795419390250;
  }

  // species with midpoint at T=1396 kelvin
  if (T < 1396) {
    // species 1: IC8H18
    result +=
      y[1] *
      (-5.20868893e+00 + 5.57202905e-02 * T - 2.63782194e-05 * T2 +
       7.31015605e-09 * T3 - 8.87486382e-13 * T4 - 2.99446875e+04 * invT) *
      0.0087541144337839;
  } else {
    // species 1: IC8H18
    result +=
      y[1] *
      (+2.61373590e+01 + 1.89502445e-02 * T - 4.31457860e-06 * T2 +
       5.01900930e-10 * T3 - 2.32801160e-14 * T4 - 4.07958177e+04 * invT) *
      0.0087541144337839;
  }

  const amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  ubms = result * RT;
}

// get mixture entropy in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBML(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real& sbml)
{
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(P / 1013250.0);
  amrex::Real sor[4]; //  temporary storage
  speciesEntropy(sor, T);

  // Compute Eq 42
  for (int id = 0; id < 4; ++id) {
    result += x[id] * (sor[id] - log((x[id] + 1e-100)) - logPratio);
  }

  sbml = result * 8.31446261815324e+07;
}

// get mixture entropy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSBMS(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real& sbms)
{
  amrex::Real result = 0;
  // Log of normalized pressure in cgs units dynes/cm^2 by Patm
  amrex::Real logPratio = log(P / 1013250.0);
  amrex::Real sor[4];  //  temporary storage
  amrex::Real x[4];    //  need a ytx conversion
  amrex::Real YOW = 0; // See Eq 4, 6 in CK Manual

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 4; i++) {
    YOW += y[i] * imw(i);
  }
  // Now compute y to x conversion
  x[0] = y[0] / (100.205000 * YOW);
  x[1] = y[1] / (114.232000 * YOW);
  x[2] = y[2] / (31.998000 * YOW);
  x[3] = y[3] / (28.014000 * YOW);
  speciesEntropy(sor, T);
  // Perform computation in Eq 42 and 43
  for (int i = 0; i < 4; i++) {
    result += x[i] * (sor[i] - log((x[i] + 1e-100)) - logPratio);
  }
  // Scale by R/W
  sbms = result * 8.31446261815324e+07 * YOW;
}

//  get temperature given internal energy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_EY(
  const amrex::Real e, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real e1, emin, emax, cv, t1, dt;
  CKUBMS(tmin, y, emin);
  CKUBMS(tmax, y, emax);
  if (e < emin) {
    // Linear Extrapolation below tmin
    CKCVBS(tmin, y, cv);
    t = tmin - (emin - e) / cv;
    ierr = 1;
    return;
  }
  if (e > emax) {
    // Linear Extrapolation above tmax
    CKCVBS(tmax, y, cv);
    t = tmax - (emax - e) / cv;
    ierr = 1;
    return;
  }
  t1 = t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (emax - emin) * (e - emin);
  }
  for (int i = 0; i < maxiter; ++i) {
    CKUBMS(t1, y, e1);
    CKCVBS(t1, y, cv);
    dt = (e - e1) / cv;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  t = t1;
  ierr = 0;
}

//  get temperature given enthalpy in mass units and mass fracs
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_HY(
  const amrex::Real h, const amrex::Real y[], amrex::Real& t, int& ierr)
{
#ifdef CONVERGENCE
  const int maxiter = 5000;
  const amrex::Real tol = 1.e-12;
#else
  const int maxiter = 200;
  const amrex::Real tol = 1.e-6;
#endif
  amrex::Real tmin = 90;   // max lower bound for thermo def
  amrex::Real tmax = 4000; // min upper bound for thermo def
  amrex::Real h1, hmin, hmax, cp, t1, dt;
  CKHBMS(tmin, y, hmin);
  CKHBMS(tmax, y, hmax);
  if (h < hmin) {
    // Linear Extrapolation below tmin
    CKCPBS(tmin, y, cp);
    t = tmin - (hmin - h) / cp;
    ierr = 1;
    return;
  }
  if (h > hmax) {
    // Linear Extrapolation above tmax
    CKCPBS(tmax, y, cp);
    t = tmax - (hmax - h) / cp;
    ierr = 1;
    return;
  }
  t1 = t;
  if (t1 < tmin || t1 > tmax) {
    t1 = tmin + (tmax - tmin) / (hmax - hmin) * (h - hmin);
  }
  for (int i = 0; i < maxiter; ++i) {
    CKHBMS(t1, y, h1);
    CKCPBS(t1, y, cp);
    dt = (h - h1) / cp;
    if (dt > 100.) {
      dt = 100.;
    } else if (dt < -100.) {
      dt = -100.;
    } else if (fabs(dt) < tol) {
      break;
    }
    t1 += dt;
  }
  t = t1;
  ierr = 0;
}

// Compute P = rhoRT/W(x)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPX(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real& P)
{
  amrex::Real XW = 0;                      //  To hold mean molecular wt
  XW += x[0] * 100.205000;                 // NC7H16
  XW += x[1] * 114.232000;                 // IC8H18
  XW += x[2] * 31.998000;                  // O2
  XW += x[3] * 28.014000;                  // N2
  P = rho * 8.31446261815324e+07 * T / XW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(y)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPY(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real& P)
{
  amrex::Real YOW = 0; //  for computing mean MW

  for (int i = 0; i < 4; i++) {
    YOW += y[i] * imw(i);
  }
  P = rho * 8.31446261815324e+07 * T * YOW; // P = rho*R*T/W
}

// Compute P = rhoRT/W(c)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPC(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real c[],
  amrex::Real& P)
{
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 100.205000; // NC7H16
  W += c[1] * 114.232000; // IC8H18
  W += c[2] * 31.998000;  // O2
  W += c[3] * 28.014000;  // N2

  for (int id = 0; id < 4; ++id) {
    sumC += c[id];
  }
  P = rho * 8.31446261815324e+07 * T * sumC / W; // P = rho*R*T/W
}

// Compute rho = PW(x)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOX(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real& rho)
{
  amrex::Real XW = 0;                        //  To hold mean molecular wt
  XW += x[0] * 100.205000;                   // NC7H16
  XW += x[1] * 114.232000;                   // IC8H18
  XW += x[2] * 31.998000;                    // O2
  XW += x[3] * 28.014000;                    // N2
  rho = P * XW / (8.31446261815324e+07 * T); // rho = P*W/(R*T)
}

// Compute rho = P*W(y)/RT
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOY(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real& rho)
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 4; i++) {
    YOW += y[i] * imw(i);
  }

  rho = P / (8.31446261815324e+07 * T * YOW); // rho = P*W/(R*T)
}

// Compute rho = P*W(c)/(R*T)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOC(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real c[],
  amrex::Real& rho)
{
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 100.205000; // NC7H16
  W += c[1] * 114.232000; // IC8H18
  W += c[2] * 31.998000;  // O2
  W += c[3] * 28.014000;  // N2

  for (int id = 0; id < 4; ++id) {
    sumC += c[id];
  }
  rho = P * W / (sumC * T * 8.31446261815324e+07); // rho = PW/(R*T)
}

// get molecular weight for all species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWT(amrex::Real wt[])
{
  get_mw(wt);
}

// given y[species]: mass fractions
// s mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWY(const amrex::Real y[], amrex::Real& wtm)
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 4; i++) {
    YOW += y[i] * imw(i);
  }

  wtm = 1.0 / YOW;
}

// given x[species]: mole fractions
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWX(const amrex::Real x[], amrex::Real& wtm)
{
  amrex::Real XW = 0;      //  see Eq 4 in CK Manual
  XW += x[0] * 100.205000; // NC7H16
  XW += x[1] * 114.232000; // IC8H18
  XW += x[2] * 31.998000;  // O2
  XW += x[3] * 28.014000;  // N2
  wtm = XW;
}

// given c[species]: molar concentration
// returns mean molecular weight (gm/mole)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWC(const amrex::Real c[], amrex::Real& wtm)
{
  // See Eq 5 in CK Manual
  amrex::Real W = 0;
  amrex::Real sumC = 0;
  W += c[0] * 100.205000; // NC7H16
  W += c[1] * 114.232000; // IC8H18
  W += c[2] * 31.998000;  // O2
  W += c[3] * 28.014000;  // N2

  for (int id = 0; id < 4; ++id) {
    sumC += c[id];
  }
  //  CK provides no guard against division by zero
  wtm = W / sumC;
}

// get Cp/R as a function of T
// for all species (Eq 19)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPOR(const amrex::Real T, amrex::Real cpor[])
{
  cp_R(cpor, T);
}

// get H/RT as a function of T
// for all species (Eq 20)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHORT(const amrex::Real T, amrex::Real hort[])
{
  speciesEnthalpy(hort, T);
}

// get S/R as a function of T
// for all species (Eq 21)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSOR(const amrex::Real T, amrex::Real sor[])
{
  speciesEntropy(sor, T);
}

// convert y[species] (mass fracs) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTX(const amrex::Real y[], amrex::Real x[])
{
  amrex::Real YOW = 0;

  for (int i = 0; i < 4; i++) {
    YOW += y[i] * imw(i);
  }

  amrex::Real YOWINV = 1.0 / YOW;

  for (int i = 0; i < 4; i++) {
    x[i] = y[i] * imw(i) * YOWINV;
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real c[])
{
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 4; i++) {
    c[i] = y[i] * imw(i);
  }
  for (int i = 0; i < 4; i++) {
    YOW += c[i];
  }

  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // Now compute conversion

  for (int i = 0; i < 4; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }
}

// convert y[species] (mass fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCR(
  const amrex::Real rho,
  amrex::Real /*T*/,
  const amrex::Real y[],
  amrex::Real c[])
{

  for (int i = 0; i < 4; i++) {
    c[i] = rho * y[i] * imw(i);
  }
}

// convert x[species] (mole fracs) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(const amrex::Real x[], amrex::Real y[])
{
  amrex::Real XW = 0; // See Eq 4, 9 in CK Manual
  // Compute mean molecular wt first
  XW += x[0] * 100.205000; // NC7H16
  XW += x[1] * 114.232000; // IC8H18
  XW += x[2] * 31.998000;  // O2
  XW += x[3] * 28.014000;  // N2
  // Now compute conversion
  amrex::Real XWinv = 1.0 / XW;
  y[0] = x[0] * 100.205000 * XWinv;
  y[1] = x[1] * 114.232000 * XWinv;
  y[2] = x[2] * 31.998000 * XWinv;
  y[3] = x[3] * 28.014000 * XWinv;
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real c[])
{
  amrex::Real PORT = P / (8.31446261815324e+07 * T); // P/RT

  // Compute conversion, see Eq 10
  for (int id = 0; id < 4; ++id) {
    c[id] = x[id] * PORT;
  }
}

// convert x[species] (mole fracs) to c[species] (molar conc)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTCR(
  const amrex::Real rho,
  const amrex::Real /*T*/,
  const amrex::Real x[],
  amrex::Real c[])
{
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 100.205000; // NC7H16
  XW += x[1] * 114.232000; // IC8H18
  XW += x[2] * 31.998000;  // O2
  XW += x[3] * 28.014000;  // N2
  ROW = rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 4; ++id) {
    c[id] = x[id] * ROW;
  }
}

// convert c[species] (molar conc) to x[species] (mole fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTX(const amrex::Real c[], amrex::Real x[])
{
  amrex::Real sumC = 0;

  // compute sum of c
  for (int id = 0; id < 4; ++id) {
    sumC += c[id];
  }

  //  See Eq 13
  amrex::Real sumCinv = 1.0 / sumC;
  for (int id = 0; id < 4; ++id) {
    x[id] = c[id] * sumCinv;
  }
}

// convert c[species] (molar conc) to y[species] (mass fracs)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCTY(const amrex::Real c[], amrex::Real y[])
{
  amrex::Real CW = 0; // See Eq 12 in CK Manual
  // compute denominator in eq 12 first
  CW += c[0] * 100.205000; // NC7H16
  CW += c[1] * 114.232000; // IC8H18
  CW += c[2] * 31.998000;  // O2
  CW += c[3] * 28.014000;  // N2
  // Now compute conversion
  amrex::Real CWinv = 1.0 / CW;
  y[0] = c[0] * 100.205000 * CWinv;
  y[1] = c[1] * 114.232000 * CWinv;
  y[2] = c[2] * 31.998000 * CWinv;
  y[3] = c[3] * 28.014000 * CWinv;
}

// get specific heat at constant volume as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVML(const amrex::Real T, amrex::Real cvml[])
{
  cv_R(cvml, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    cvml[id] *= 8.31446261815324e+07;
  }
}

// get specific heat at constant pressure as a
// function of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPML(const amrex::Real T, amrex::Real cpml[])
{
  cp_R(cpml, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    cpml[id] *= 8.31446261815324e+07;
  }
}

// get internal energy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUML(const amrex::Real T, amrex::Real uml[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesInternalEnergy(uml, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    uml[id] *= RT;
  }
}

// get enthalpy as a function
// of T for all species (molar units)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHML(const amrex::Real T, amrex::Real hml[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T
  speciesEnthalpy(hml, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    hml[id] *= RT;
  }
}

// Returns the standard-state entropies in molar units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSML(const amrex::Real T, amrex::Real sml[])
{
  speciesEntropy(sml, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    sml[id] *= 8.31446261815324e+07;
  }
}

// Returns the specific heats at constant volume
// in mass units (Eq. 29)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVMS(const amrex::Real T, amrex::Real cvms[])
{
  cv_R(cvms, T);
  // multiply by R/molecularweight
  cvms[0] *= 8.297452839831585e+05; // NC7H16
  cvms[1] *= 7.278575721473177e+05; // IC8H18
  cvms[2] *= 2.598431970171023e+06; // O2
  cvms[3] *= 2.967966951578939e+06; // N2
}

// Returns the specific heats at constant pressure
// in mass units (Eq. 26)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPMS(const amrex::Real T, amrex::Real cpms[])
{
  cp_R(cpms, T);
  // multiply by R/molecularweight
  cpms[0] *= 8.297452839831585e+05; // NC7H16
  cpms[1] *= 7.278575721473177e+05; // IC8H18
  cpms[2] *= 2.598431970171023e+06; // O2
  cpms[3] *= 2.967966951578939e+06; // N2
}

// Returns internal energy in mass units (Eq 30.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(const amrex::Real T, amrex::Real ums[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesInternalEnergy(ums, T);

  for (int i = 0; i < 4; i++) {
    ums[i] *= RT * imw(i);
  }
}

// Returns enthalpy in mass units (Eq 27.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHMS(const amrex::Real T, amrex::Real hms[])
{
  amrex::Real RT = 8.31446261815324e+07 * T; // R*T

  speciesEnthalpy(hms, T);

  for (int i = 0; i < 4; i++) {
    hms[i] *= RT * imw(i);
  }
}

// Returns the entropies in mass units (Eq 28.)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKSMS(const amrex::Real T, amrex::Real sms[])
{
  speciesEntropy(sms, T);
  // multiply by R/molecularweight
  sms[0] *= 8.297452839831585e+05; // NC7H16
  sms[1] *= 7.278575721473177e+05; // IC8H18
  sms[2] *= 2.598431970171023e+06; // O2
  sms[3] *= 2.967966951578939e+06; // N2
}

// GPU version of productionRate: no more use of thermo namespace vectors
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
comp_qfqr(
  amrex::Real* qf,
  amrex::Real* qr,
  const amrex::Real* sc,
  const amrex::Real* /*sc_qss*/,
  const amrex::Real T,
  const amrex::Real invT,
  const amrex::Real logT)
{

  // reaction 0: N2 + O2 <=> N2 + O2
  qf[0] = sc[2] * sc[3];
  qr[0] = sc[2] * sc[3];

  // reaction 1: IC8H18 + NC7H16 <=> IC8H18 + NC7H16
  qf[1] = sc[0] * sc[1];
  qr[1] = sc[0] * sc[1];

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 4; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[4];
  gibbs(g_RT, T);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 * invT;
  amrex::Real refCinv = 1 / refC;

  // Evaluate the kfs
  amrex::Real k_f, Corr;

  // reaction 0:  N2 + O2 <=> N2 + O2
  k_f = 6000000000 * exp(-(23147.9666270171) * invT);
  qf[0] *= k_f;
  qr[0] *= k_f * exp(-(g_RT[2] - g_RT[2] + g_RT[3] - g_RT[3]));
  // reaction 1:  IC8H18 + NC7H16 <=> IC8H18 + NC7H16
  k_f = 6000000000 * exp(-(23147.9666270171) * invT);
  qf[1] *= k_f;
  qr[1] *= k_f * exp(-(g_RT[0] - g_RT[0] + g_RT[1] - g_RT[1]));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
productionRate(amrex::Real* wdot, const amrex::Real* sc, const amrex::Real T)
{
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  const amrex::Real refC = 101325 / 8.31446 * invT;
  const amrex::Real refCinv = 1 / refC;

  for (int i = 0; i < 4; ++i) {
    wdot[i] = 0.0;
  }

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int i = 0; i < 4; ++i) {
    mixture += sc[i];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[4];
  gibbs(g_RT, T);

  {
    // reaction 0:  N2 + O2 <=> N2 + O2
    const amrex::Real k_f = 6000000000 * exp(-(23147.9666270171) * invT);
    const amrex::Real qf = k_f * (sc[2] * sc[3]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[2] - g_RT[2] + g_RT[3] - g_RT[3])) * (sc[2] * sc[3]);
    const amrex::Real qdot = qf - qr;
    wdot[2] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[3] += qdot;
  }

  {
    // reaction 1:  IC8H18 + NC7H16 <=> IC8H18 + NC7H16
    const amrex::Real k_f = 6000000000 * exp(-(23147.9666270171) * invT);
    const amrex::Real qf = k_f * (sc[0] * sc[1]);
    const amrex::Real qr =
      k_f * exp(-(g_RT[0] - g_RT[0] + g_RT[1] - g_RT[1])) * (sc[0] * sc[1]);
    const amrex::Real qdot = qf - qr;
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[1] += qdot;
  }
}

// compute the production rate for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real T, amrex::Real C[], amrex::Real wdot[])
{

  // convert to SI
  for (int id = 0; id < 4; ++id) {
    C[id] *= 1.0e6;
  }

  // convert to chemkin units
  productionRate(wdot, C, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    C[id] *= 1.0e-6;
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real wdot[])
{
  amrex::Real c[4]; // temporary storage
  amrex::Real YOW = 0;
  amrex::Real PWORT;

  // Compute inverse of mean molecular wt first
  for (int i = 0; i < 4; i++) {
    YOW += y[i] * imw(i);
  }
  // PW/RT (see Eq. 7)
  PWORT = P / (YOW * 8.31446261815324e+07 * T);
  // multiply by 1e6 so c goes to SI
  PWORT *= 1e6;
  // Now compute conversion (and go to SI)
  for (int i = 0; i < 4; i++) {
    c[i] = PWORT * y[i] * imw(i);
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given P, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXP(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real wdot[])
{
  amrex::Real c[4]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 4; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mass fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWYR(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real y[],
  amrex::Real wdot[])
{
  amrex::Real c[4]; // temporary storage

  // See Eq 8 with an extra 1e6 so c goes to SI
  for (int i = 0; i < 4; i++) {
    c[i] = 1e6 * rho * y[i] * imw(i);
  }

  // call productionRate
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

// Returns the molar production rate of species
// Given rho, T, and mole fractions
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWXR(
  const amrex::Real rho,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real wdot[])
{
  amrex::Real c[4];   // temporary storage
  amrex::Real XW = 0; // See Eq 4, 11 in CK Manual
  amrex::Real ROW;
  // Compute mean molecular wt first
  XW += x[0] * 100.205000; // NC7H16
  XW += x[1] * 114.232000; // IC8H18
  XW += x[2] * 31.998000;  // O2
  XW += x[3] * 28.014000;  // N2
  // Extra 1e6 factor to take c to SI
  ROW = 1e6 * rho / XW;

  // Compute conversion, see Eq 11
  for (int id = 0; id < 4; ++id) {
    c[id] = x[id] * ROW;
  }

  // convert to chemkin units
  productionRate(wdot, c, T);

  // convert to chemkin units
  for (int id = 0; id < 4; ++id) {
    wdot[id] *= 1.0e-6;
  }
}

//  species unit charge number
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRG(int kcharge[])
{
  kcharge[0] = 0; // NC7H16
  kcharge[1] = 0; // IC8H18
  kcharge[2] = 0; // O2
  kcharge[3] = 0; // N2
}

//  species charge per unit mass
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCHRGMASS(amrex::Real zk[])
{

  int kchrg[4];
  CKCHRG(kchrg);

  for (int id = 0; id < 4; ++id) {
    zk[id] = 6.02214076e+23 * 1.60217663e-19 * kchrg[id] * imw(id);
  }
}

// compute an approx to the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian_precond(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real T, const int HP)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 25; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[4];
  for (auto& val : wdot) {
    val = 0.0;
  }

  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real logT = log(T);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 / T;
  amrex::Real refCinv = 1.0 / refC;

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int k = 0; k < 4; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[4];
  gibbs(g_RT, T);

  // compute the species enthalpy
  amrex::Real h_RT[4];
  speciesEnthalpy(h_RT, T);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[4];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 0: N2 + O2 <=> N2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[3];
  k_f = 6000000000 * exp(-(23147.9666270171) * invT);
  dlnkfdT = (23147.9666270171) * invT2;
  // reverse
  phi_r = sc[2] * sc[3];
  Kc = exp(g_RT[2] - g_RT[2] + g_RT[3] - g_RT[3]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[2] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  // d()/d[O2]
  dqdci = +k_f * sc[3] - k_r * sc[3];
  // d()/d[N2]
  dqdci = +k_f * sc[2] - k_r * sc[2];
  // d()/dT

  // reaction 1: IC8H18 + NC7H16 <=> IC8H18 + NC7H16
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[1];
  k_f = 6000000000 * exp(-(23147.9666270171) * invT);
  dlnkfdT = (23147.9666270171) * invT2;
  // reverse
  phi_r = sc[0] * sc[1];
  Kc = exp(g_RT[0] - g_RT[0] + g_RT[1] - g_RT[1]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[1]) + (h_RT[0] + h_RT[1]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  // d()/d[NC7H16]
  dqdci = +k_f * sc[1] - k_r * sc[1];
  // d()/d[IC8H18]
  dqdci = +k_f * sc[0] - k_r * sc[0];
  // d()/dT

  amrex::Real c_R[4], dcRdT[4], e_RT[4];
  amrex::Real* eh_RT;
  if (HP == 1) {
    cp_R(c_R, T);
    dcvpRdT(dcRdT, T);
    eh_RT = &h_RT[0];
  } else {
    cv_R(c_R, T);
    dcvpRdT(dcRdT, T);
    speciesInternalEnergy(e_RT, T);
    eh_RT = &e_RT[0];
  }

  amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT = 0.0, dehmixdT = 0.0;
  for (int k = 0; k < 4; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[20 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 4; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 4; ++m) {
      dehmixdc += eh_RT[m] * J[k * 5 + m];
    }
    J[k * 5 + 4] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[24] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
#else
  amrex::Abort();
#endif
}

// compute an approx to the reaction Jacobian (for preconditioning)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT_SIMPLIFIED(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real* Tp, const int* HP)
{
  amrex::Real c[4];

  for (int k = 0; k < 4; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian_precond(J, c, *Tp, *HP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 4; k++) {
    J[20 + k] *= 1.e-6;
    J[k * 5 + 4] *= 1.e6;
  }
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
aJacobian(
  amrex::Real* J, const amrex::Real* sc, const amrex::Real T, const int consP)
{

#if defined(PELE_COMPILE_AJACOBIAN) || !defined(AMREX_USE_HIP)
  for (int i = 0; i < 25; i++) {
    J[i] = 0.0;
  }

  amrex::Real wdot[4];
  for (auto& val : wdot) {
    val = 0.0;
  }

  const amrex::Real invT = 1.0 / T;
  const amrex::Real invT2 = invT * invT;
  const amrex::Real logT = log(T);

  // reference concentration: P_atm / (RT) in inverse mol/m^3
  amrex::Real refC = 101325 / 8.31446 / T;
  amrex::Real refCinv = 1.0 / refC;

  // compute the mixture concentration
  amrex::Real mixture = 0.0;
  for (int k = 0; k < 4; ++k) {
    mixture += sc[k];
  }

  // compute the Gibbs free energy
  amrex::Real g_RT[4];
  gibbs(g_RT, T);

  // compute the species enthalpy
  amrex::Real h_RT[4];
  speciesEnthalpy(h_RT, T);

  amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
  amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
  amrex::Real dqdci, dcdc_fac, dqdc[4];
  amrex::Real Pr, fPr, F, k_0, logPr;
  amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
  amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
  amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
  amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
  const amrex::Real ln10 = log(10.0);
  const amrex::Real log10e = 1.0 / log(10.0);
  // reaction 0: N2 + O2 <=> N2 + O2
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[2] * sc[3];
  k_f = 6000000000 * exp(-(23147.9666270171) * invT);
  dlnkfdT = (23147.9666270171) * invT2;
  // reverse
  phi_r = sc[2] * sc[3];
  Kc = exp(g_RT[2] - g_RT[2] + g_RT[3] - g_RT[3]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[2] + h_RT[3]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  // d()/d[O2]
  dqdci = +k_f * sc[3] - k_r * sc[3];
  // d()/d[N2]
  dqdci = +k_f * sc[2] - k_r * sc[2];
  // d()/dT

  // reaction 1: IC8H18 + NC7H16 <=> IC8H18 + NC7H16
  // a non-third-body and non-pressure-fall-off reaction
  // forward
  phi_f = sc[0] * sc[1];
  k_f = 6000000000 * exp(-(23147.9666270171) * invT);
  dlnkfdT = (23147.9666270171) * invT2;
  // reverse
  phi_r = sc[0] * sc[1];
  Kc = exp(g_RT[0] - g_RT[0] + g_RT[1] - g_RT[1]);
  k_r = k_f / Kc;
  dlnKcdT = invT * (-(h_RT[0] + h_RT[1]) + (h_RT[0] + h_RT[1]));
  dkrdT = (dlnkfdT - dlnKcdT) * k_r;
  // rate of progress
  q = k_f * phi_f - k_r * phi_r;
  dqdT = (dlnkfdT * k_f * phi_f - dkrdT * phi_r);
  // update wdot
  // d()/d[NC7H16]
  dqdci = +k_f * sc[1] - k_r * sc[1];
  // d()/d[IC8H18]
  dqdci = +k_f * sc[0] - k_r * sc[0];
  // d()/dT

  amrex::Real c_R[4], dcRdT[4], e_RT[4];
  amrex::Real* eh_RT;
  if (consP == 1) {
    cp_R(c_R, T);
    dcvpRdT(dcRdT, T);
    eh_RT = &h_RT[0];
  } else {
    cv_R(c_R, T);
    dcvpRdT(dcRdT, T);
    speciesInternalEnergy(e_RT, T);
    eh_RT = &e_RT[0];
  }

  amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT = 0.0, dehmixdT = 0.0;
  for (int k = 0; k < 4; ++k) {
    cmix += c_R[k] * sc[k];
    dcmixdT += dcRdT[k] * sc[k];
    ehmix += eh_RT[k] * wdot[k];
    dehmixdT += invT * (c_R[k] - eh_RT[k]) * wdot[k] + eh_RT[k] * J[20 + k];
  }

  amrex::Real cmixinv = 1.0 / cmix;
  amrex::Real tmp1 = ehmix * cmixinv;
  amrex::Real tmp3 = cmixinv * T;
  amrex::Real tmp2 = tmp1 * tmp3;
  amrex::Real dehmixdc;
  // dTdot/d[X]
  for (int k = 0; k < 4; ++k) {
    dehmixdc = 0.0;
    for (int m = 0; m < 4; ++m) {
      dehmixdc += eh_RT[m] * J[k * 5 + m];
    }
    J[k * 5 + 4] = tmp2 * c_R[k] - tmp3 * dehmixdc;
  }
  // dTdot/dT
  J[24] = -tmp1 + tmp2 * dcmixdT - tmp3 * dehmixdT;
#else
  amrex::Abort();
#endif
}

// compute the reaction Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT(
  amrex::Real* J,
  const amrex::Real* sc,
  const amrex::Real* Tp,
  const int* consP)
{
  amrex::Real c[4];

  for (int k = 0; k < 4; k++) {
    c[k] = 1.e6 * sc[k];
  }

  aJacobian(J, c, *Tp, *consP);

  // dwdot[k]/dT
  // dTdot/d[X]
  for (int k = 0; k < 4; k++) {
    J[20 + k] *= 1.e-6;
    J[k * 5 + 4] *= 1.e6;
  }
}

// Transport function declarations

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENIMC(int* LENIMC)
{
  *LENIMC = 16;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* LENRMC)
{
  *LENRMC = 412;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNO(int* NO)
{
  *NO = 4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKK(int* KK)
{
  *KK = 4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLITE(int* NLITE)
{
  *NLITE = 0;
}

// Patm in ergs/cm3
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPATM(amrex::Real* PATM)
{
  *PATM = 0.1013250000000000E+07;
}

// the molecular weights in g/mol
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetWT(amrex::Real* WT)
{
  WT[0] = 1.00205000E+02;
  WT[1] = 1.14232000E+02;
  WT[2] = 3.19980000E+01;
  WT[3] = 2.80140000E+01;
}

// the lennard-jones potential well depth eps/kb in K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* EPS)
{
  EPS[0] = 4.59600000E+02;
  EPS[1] = 4.58500000E+02;
  EPS[2] = 1.07400000E+02;
  EPS[3] = 9.75300000E+01;
}

// the lennard-jones collision diameter in Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* SIG)
{
  SIG[0] = 6.25300000E+00;
  SIG[1] = 6.41400000E+00;
  SIG[2] = 3.45800000E+00;
  SIG[3] = 3.62100000E+00;
}

// the dipole moment in Debye
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetDIP(amrex::Real* DIP)
{
  DIP[0] = 0.00000000E+00;
  DIP[1] = 0.00000000E+00;
  DIP[2] = 0.00000000E+00;
  DIP[3] = 0.00000000E+00;
}

// the polarizability in cubic Angstroms
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* POL)
{
  POL[0] = 0.00000000E+00;
  POL[1] = 0.00000000E+00;
  POL[2] = 1.60000000E+00;
  POL[3] = 1.76000000E+00;
}

// the rotational relaxation collision number at 298 K
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* ZROT)
{
  ZROT[0] = 1.00000000E+00;
  ZROT[1] = 1.00000000E+00;
  ZROT[2] = 3.80000000E+00;
  ZROT[3] = 4.00000000E+00;
}

// 0: monoatomic, 1: linear, 2: nonlinear
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* NLIN)
{
  NLIN[0] = 2;
  NLIN[1] = 2;
  NLIN[2] = 1;
  NLIN[3] = 1;
}

// Poly fits for the viscosities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real* COFETA)
{
  COFETA[0] = -2.43263310E+01;
  COFETA[1] = 4.51130748E+00;
  COFETA[2] = -4.35645301E-01;
  COFETA[3] = 1.63053081E-02;
  COFETA[4] = -2.43281428E+01;
  COFETA[5] = 4.52000399E+00;
  COFETA[6] = -4.37006790E-01;
  COFETA[7] = 1.63725669E-02;
  COFETA[8] = -1.60066453E+01;
  COFETA[9] = 2.16753735E+00;
  COFETA[10] = -1.97226850E-01;
  COFETA[11] = 8.50065468E-03;
  COFETA[12] = -1.55270224E+01;
  COFETA[13] = 1.92766908E+00;
  COFETA[14] = -1.66518287E-01;
  COFETA[15] = 7.19100649E-03;
}

// Poly fits for the conductivities, dim NO*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* COFLAM)
{
  COFLAM[0] = -2.37495362E+01;
  COFLAM[1] = 9.84951804E+00;
  COFLAM[2] = -9.67095199E-01;
  COFLAM[3] = 3.34019439E-02;
  COFLAM[4] = -2.64785178E+01;
  COFLAM[5] = 1.10165416E+01;
  COFLAM[6] = -1.12955906E+00;
  COFLAM[7] = 4.08248190E-02;
  COFLAM[8] = -2.11868696E+00;
  COFLAM[9] = 2.98568648E+00;
  COFLAM[10] = -2.86879126E-01;
  COFLAM[11] = 1.23850874E-02;
  COFLAM[12] = 7.60995874E+00;
  COFLAM[13] = -1.18418483E+00;
  COFLAM[14] = 3.03558391E-01;
  COFLAM[15] = -1.54159450E-02;
}

// Poly fits for the diffusion coefficients, dim NO*KK*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* COFD)
{
  COFD[0] = -2.32987245E+01;
  COFD[1] = 5.16842641E+00;
  COFD[2] = -3.96171273E-01;
  COFD[3] = 1.48789792E-02;
  COFD[4] = -2.33591632E+01;
  COFD[5] = 5.17059469E+00;
  COFD[6] = -3.96543827E-01;
  COFD[7] = 1.48982973E-02;
  COFD[8] = -2.03379496E+01;
  COFD[9] = 4.78258755E+00;
  COFD[10] = -3.88556217E-01;
  COFD[11] = 1.61092477E-02;
  COFD[12] = -1.99122802E+01;
  COFD[13] = 4.63932958E+00;
  COFD[14] = -3.71390161E-01;
  COFD[15] = 1.54205723E-02;
  COFD[16] = -2.33591632E+01;
  COFD[17] = 5.17059469E+00;
  COFD[18] = -3.96543827E-01;
  COFD[19] = 1.48982973E-02;
  COFD[20] = -2.34214262E+01;
  COFD[21] = 5.17276715E+00;
  COFD[22] = -3.96917346E-01;
  COFD[23] = 1.49176769E-02;
  COFD[24] = -2.03757148E+01;
  COFD[25] = 4.77900446E+00;
  COFD[26] = -3.88127029E-01;
  COFD[27] = 1.60920270E-02;
  COFD[28] = -1.99500592E+01;
  COFD[29] = 4.63656074E+00;
  COFD[30] = -3.71069499E-01;
  COFD[31] = 1.54081752E-02;
  COFD[32] = -2.03379496E+01;
  COFD[33] = 4.78258755E+00;
  COFD[34] = -3.88556217E-01;
  COFD[35] = 1.61092477E-02;
  COFD[36] = -2.03757148E+01;
  COFD[37] = 4.77900446E+00;
  COFD[38] = -3.88127029E-01;
  COFD[39] = 1.60920270E-02;
  COFD[40] = -1.47079535E+01;
  COFD[41] = 3.10657376E+00;
  COFD[42] = -1.85922460E-01;
  COFD[43] = 7.92680827E-03;
  COFD[44] = -1.44285962E+01;
  COFD[45] = 2.99858376E+00;
  COFD[46] = -1.72232643E-01;
  COFD[47] = 7.34804765E-03;
  COFD[48] = -1.99122802E+01;
  COFD[49] = 4.63932958E+00;
  COFD[50] = -3.71390161E-01;
  COFD[51] = 1.54205723E-02;
  COFD[52] = -1.99500592E+01;
  COFD[53] = 4.63656074E+00;
  COFD[54] = -3.71069499E-01;
  COFD[55] = 1.54081752E-02;
  COFD[56] = -1.44285962E+01;
  COFD[57] = 2.99858376E+00;
  COFD[58] = -1.72232643E-01;
  COFD[59] = 7.34804765E-03;
  COFD[60] = -1.42056777E+01;
  COFD[61] = 2.91297621E+00;
  COFD[62] = -1.61544771E-01;
  COFD[63] = 6.90271324E-03;
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* /*KTDIF*/)
{
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* /*COFTD*/)
{
}

// compute the critical parameters for each species
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS(
  amrex::Real* Tci, amrex::Real* ai, amrex::Real* bi, amrex::Real* acentric_i)
{

  amrex::Real EPS[4];
  amrex::Real SIG[4];
  amrex::Real wt[4];
  amrex::Real Rcst = 83.144598; // in bar [CGS] !
  amrex::Real avogadro = 6.02214199e23;
  amrex::Real boltzmann = 1.3806503e-16; // we work in CGS

  egtransetEPS(EPS);
  egtransetSIG(SIG);
  get_mw(wt);

  // species 0: NC7H16
  Tci[0] = 1.316 * EPS[0];
  ai[0] = (5.55 * avogadro * avogadro * EPS[0] * boltzmann * 1e-24 * SIG[0] *
           SIG[0] * SIG[0]) /
          (wt[0] * wt[0]);
  bi[0] = 0.855 * avogadro * 1e-24 * SIG[0] * SIG[0] * SIG[0] / (wt[0]);
  acentric_i[0] = 0.0;

  // species 1: IC8H18
  Tci[1] = 1.316 * EPS[1];
  ai[1] = (5.55 * avogadro * avogadro * EPS[1] * boltzmann * 1e-24 * SIG[1] *
           SIG[1] * SIG[1]) /
          (wt[1] * wt[1]);
  bi[1] = 0.855 * avogadro * 1e-24 * SIG[1] * SIG[1] * SIG[1] / (wt[1]);
  acentric_i[1] = 0.0;

  // species 2: O2
  // Imported from NIST
  Tci[2] = 154.581000;
  ai[2] = 1e6 * 0.42748 * Rcst * Rcst * Tci[2] * Tci[2] /
          (31.998800 * 31.998800 * 50.430466);
  bi[2] = 0.08664 * Rcst * Tci[2] / (31.998800 * 50.430466);
  acentric_i[2] = 0.022200;

  // species 3: N2
  // Imported from NIST
  Tci[3] = 126.192000;
  ai[3] = 1e6 * 0.42748 * Rcst * Rcst * Tci[3] * Tci[3] /
          (28.013400 * 28.013400 * 33.958000);
  bi[3] = 0.08664 * Rcst * Tci[3] / (28.013400 * 33.958000);
  acentric_i[3] = 0.037200;
}

// compute the critical parameter quantities for each species for SRK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS_SRK(
  amrex::Real* sqrtOneOverTc,
  amrex::Real* sqrtAsti,
  amrex::Real* Bi,
  amrex::Real* Fomega)
{

  // species 0: NC7H16
  // Computed from Lennard-Jones
  sqrtOneOverTc[0] = 4.0661373572610e-02;
  sqrtAsti[0] = 5.5766397379575e+04;
  Bi[0] = 1.2562989577857e+00;
  Fomega[0] = 4.8508000000000e-01;

  // species 1: IC8H18
  // Computed from Lennard-Jones
  sqrtOneOverTc[1] = 4.0710120264334e-02;
  sqrtAsti[1] = 5.0759186273762e+04;
  Bi[1] = 1.1893678422487e+00;
  Fomega[1] = 4.8508000000000e-01;

  // species 2: O2
  // Imported from NIST
  sqrtOneOverTc[2] = 8.0430717653241e-02;
  sqrtAsti[2] = 3.6980081924740e+04;
  Bi[2] = 6.9005216444999e-01;
  Fomega[2] = 5.1945301904908e-01;

  // species 3: N2
  // Imported from NIST
  sqrtOneOverTc[3] = 8.9019282240563e-02;
  sqrtAsti[3] = 4.2022983929562e+04;
  Bi[3] = 9.5560052221563e-01;
  Fomega[3] = 5.4259343186608e-01;
}

// gauss-jordan solver external routine
// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/)
{
  amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

// Replace this routine with the one generated by the Gauss Jordan solver of DW
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve_simplified(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/)
{
  amrex::Abort(
    "sgjsolve_simplified not implemented, choose a different solver ");
}

#endif
