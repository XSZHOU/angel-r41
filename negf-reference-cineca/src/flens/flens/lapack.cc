/*
 *   Copyright (c) 2007, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cassert>
#include <complex>
#include <flens/lapack.h>

#ifdef VECLIB
#    include <Accelerate/Accelerate.h>
#elif defined MKL
#    ifdef MAC
#        include <Intel_MKL/mkl_lapack.h>
#    else
#        include <mkl_lapack.h>
#    endif
#else
#    include <gsl_cblas.h>
//#    include <cblas.h>
#endif

#ifdef NOUNDERSCORE
	#define F77_Name(name) name
#else 	
	#define F77_Name(name)  name##_
#endif


namespace flens {

extern "C" {

    void
    F77_Name(sgetrf)(int *m, int *n, float *a, int *lda, int *ipiv, int *info);

    void
    F77_Name(dgetrf)(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

    void
    F77_Name(zgetrf)(int *m, int *n, std::complex<double> *a, int *lda,
            int *ipiv, int *info);

    void
    F77_Name(sgbtrf)(int *m, int *n, int *kl, int *ku,
            float *ab, int *ldab, int *ipiv, int *info);

    void
    F77_Name(dgbtrf)(int *m, int *n, int *kl, int *ku,
            double *ab, int *ldab, int *ipiv, int *info);

    void
    F77_Name(sgetri)(int *n, float *a, int *lda, const int *ipiv, float *work,
            int *lwork, int *info);

    void
    F77_Name(dgetri)(int *n, double *a, int *lda, const int *ipiv, double *work,
            int *lwork, int *info);

    void
    F77_Name(zgetri)(int *n, std::complex<double> *a, int *lda, const int *ipiv,
            std::complex<double> *work, int *lwork, int *info);

    void
    F77_Name(sgetrs)(char *trans, int *n, int *nrhs, const float *a, int *lda,
            const int *ipiv, float *b, int *ldb, int *info);

    void
    F77_Name(dgetrs)(char *trans, int *n, int *nrhs, const double *a, int *lda,
            const int *ipiv, double *b, int *ldb, int *info);

    void
    F77_Name(sgbtrs)(char *trans, int *n, int *kl, int *ku, int *nrhs, const float *ab,
            int *ldab, const int *ipiv, float *b, int *ldb, int *info);

    void
    F77_Name(dgbtrs)(char *trans, int *n, int *kl, int *ku, int *nrhs, const double *ab,
            int *ldab, const int *ipiv, double *b, int *ldb, int *info);

    void
    F77_Name(sgesv)(int *n, int *nrhs, float *a, int *lda,
           int *ipiv, float *b, int *ldb, int *info);

    void
    F77_Name(dgesv)(int *n, int *nrhs, double *a, int *lda,
           int *ipiv, double *b, int *ldb, int *info);

    void
    F77_Name(zgesv)(int *n, int *nrhs, std::complex<double> *a, int *lda,
           int *ipiv, std::complex<double> *b, int *ldb, int *info);

    void
    F77_Name(sgbsv)(int *n, int *kl, int *ku, int *nrhs, float *ab, int *ldab, int *ipiv,
           float *b, int *ldb, int *info);

    void
    F77_Name(dgbsv)(int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab,
           int *ipiv, double *b, int *ldb, int *info);

    void
    F77_Name(zgbsv)(int *n, int *kl, int *ku, int *nrhs, std::complex<double> *ab,
           int *ldab, int *ipiv, std::complex<double> *b, int *ldb, int *info);

    void
    F77_Name(strtrs)(char *uplo, char *trans, char *diag, int *n, int *nrhs,
            const float *a, int *lda, float *b, int *ldb, int *info);

    void
    F77_Name(dtrtrs)(char *uplo, char *trans, char *diag, int *n, int *nrhs,
            const double *a, int *lda, double *b, int *ldb, int *info);

    void
    F77_Name(sgeqrf)(int *m, int *n, float *a, int *lda, float *tau,
            float *work, int *lwork, int *info);

    void
    F77_Name(dgeqrf)(int *m, int *n, double *a, int *lda, double *tau,
            double *work, int *lwork, int *info);

    void
    F77_Name(sorgqr)(int *m, int *n, int *k, float *a, int *lda, const float *tau,
            float *work, int *lwork, int *info);

    void
    F77_Name(dorgqr)(int *m, int *n, int *k, double *a, int *lda, const double *tau,
            double *work, int *lwork, int *info);

    void
    F77_Name(sormqr)(char *side, char *trans, int *m, int *n, int *k,
            const float *a, int *lda, const float *tau, float *c, int *ldc,
            float *work, int *lwork, int *info);

    void
    F77_Name(dormqr)(char *side, char *trans, int *m, int *n, int *k,
            const double *a, int *lda, const double *tau, double *c, int *ldc,
            double *work, int *lwork, int *info);

    void
    F77_Name(sgels)(char *trans, int *m, int *n, int *nrhs, float *a, int *lda,
           float *b, int *ldb, float *work, int *lwork, int *info);

    void
    F77_Name(dgels)(char *trans, int *m, int *n, int *nrhs, double *a, int *lda,
           double *b, int *ldb, double *work, int *lwork, int *info);

    void
    F77_Name(zgels)(char *trans, int *m, int *n, int *nrhs, std::complex<double> *a,
           int *lda, std::complex<double> *b, int *ldb,
           std::complex<double> *work, int *lwork, int *info);

    void
    F77_Name(sgelss)(int *m, int *n, int *nrhs,
            float *a, int *lda, float *b, int *ldb,
            float *s, float *rcond, int *rank,
            float *work, int *lwork, int *info);

    void
    F77_Name(dgelss)(int *m, int *n, int *nrhs,
            double *a, int *lda, double *b, int *ldb,
            double *s, double *rcond, int *rank,
            double *work, int *lwork, int *info);

    void
    F77_Name(zgelss)(int *m, int *n, int *nrhs,
            std::complex<double> *a, int *lda,
            std::complex<double> *b, int *ldb,
            std::complex<double> *s, std::complex<double> *rcond, int *rank,
            std::complex<double> *work, int *lwork, int *info);

    void
    F77_Name(sgeev)(char *jobvl, char *jobvr, int *n, float *a, int *lda,
           float *wr, float *wi,
           float *vl, int *ldvl,
           float *vr, int *ldvr,
           float *work, int *lwork, int *info);

    void
    F77_Name(dgeev)(char *jobvl, char *jobvr, int *n, double *a, int *lda,
           double *wr, double *wi,
           double *vl, int *ldvl,
           double *vr, int *ldvr,
           double *work, int *lwork, int *info);

    void
    F77_Name(cgeev)(char *jobvl, char *jobvr, int *n, std::complex<float> *a, int *lda,
           std::complex<float> *w,
           std::complex<float> *vl, int *ldvl,
           std::complex<float> *vr, int *ldvr,
           std::complex<float> *work,int *lwork,float *rwork,int *info);

    void
    F77_Name(zgeev)(char *jobvl, char *jobvr, int *n, std::complex<double> *a, int *lda,
           std::complex<double> *w,
           std::complex<double> *vl, int *ldvl,
           std::complex<double> *vr, int *ldvr,
           std::complex<double> *work,int *lwork,double *rwork,int *info);

    void
    F77_Name(sgesvd)(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda,
            float *s, float *u, int *ldu, float *vt, int *ldvt,
            float *work, int *lwork, int *info);

    void
    F77_Name(dgesvd)(char *jobu, char *jobvt, int *m, int *n, double *a, int *lda,
            double *s, double *u, int *ldu, double *vt, int *ldvt,
            double *work, int *lwork, int *info);
    void
    F77_Name(dgecon)(char *norm, int *n, double *a, int *lda, double *anorm, 
            double *rcond, double *work, int *iwork, int *info);
}

int
getrf(int m, int n, float *a, int lda, int *ipiv)
{
    int info;
    F77_Name(sgetrf)(&m, &n, a, &lda, ipiv, &info);
    return info;
}

int
getrf(int m, int n, double *a, int lda, int *ipiv)
{
    int info;
    F77_Name(dgetrf)(&m, &n, a, &lda, ipiv, &info);
    return info;
}

int
getrf(int m, int n, std::complex<double> *a, int lda, int *ipiv)
{
    int info;
    F77_Name(zgetrf)(&m, &n, a, &lda, ipiv, &info);
    return info;
}

int
gbtrf(int m, int n, int kl, int ku, float *ab, int ldab, int *ipiv)
{
    int info;
    F77_Name(sgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

int
gbtrf(int m, int n, int kl, int ku, double *ab, int ldab, int *ipiv)
{
    int info;
    F77_Name(dgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}

int
getri(int n, float *a, int lda, const int *ipiv,
      float *work, int lwork)
{
    int info;
    F77_Name(sgetri)(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

int
getri(int n, double *a, int lda, const int *ipiv,
      double *work, int lwork)
{
    int info;
    F77_Name(dgetri)(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

int
getri(int n, std::complex<double> *a, int lda, const int *ipiv,
      std::complex<double> *work, int lwork)
{
    int info;
    F77_Name(zgetri)(&n, a, &lda, ipiv, work, &lwork, &info);
    return info;
}

int
getrs(Transpose trans, int n, int nrhs, const float *a, int lda,
      const int *ipiv, float *b, int ldb)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    F77_Name(sgetrs)(&_trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
getrs(Transpose trans, int n, int nrhs, const double *a, int lda,
      const int *ipiv, double *b, int ldb)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    F77_Name(dgetrs)(&_trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
gbtrs(Transpose trans, int n, int kl, int ku, int nrhs,
      const float *ab, int ldab, const int *ipiv, float *b, int ldb)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';
    F77_Name(sgbtrs)(&_trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
gbtrs(Transpose trans, int n, int kl, int ku, int nrhs,
      const double *ab, int ldab, const int *ipiv, double *b, int ldb)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';
    F77_Name(dgbtrs)(&_trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
gesv(int n, int nrhs, float *a, int lda, int *ipiv, float *b, int ldb)
{
    int info;
    F77_Name(sgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
gesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
    int info;
    F77_Name(dgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
gesv(int n, int nrhs, std::complex<double> *a, int lda, int *ipiv,
     std::complex<double> *b, int ldb)
{
    int info;
    F77_Name(zgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
    return info;
}

int
gbsv(int n, int kl, int ku, int nrhs, float *ab, int ldab,
     int *ipiv, float *b, int ldb)
{
    int info;
    F77_Name(sgbsv)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
gbsv(int n, int kl, int ku, int nrhs, double *ab, int ldab,
     int *ipiv, double *b, int ldb)
{
    int info;
    F77_Name(dgbsv)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
gbsv(int n, int kl, int ku, int nrhs, std::complex<double> *ab, int ldab,
     int *ipiv, std::complex<double> *b, int ldb)
{
    int info;
    F77_Name(zgbsv)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
    return info;
}

int
trtrs(StorageUpLo upLo, Transpose trans, UnitDiag diag, int n, int nrhs,
      const float *a, int lda, float *b, int ldb)
{
    int info;
    char _upLo = (upLo==Upper) ? 'U' : 'L';
    char _trans = (trans==NoTrans) ? 'N' : 'T';
    char _diag = (diag==Unit) ? 'U' : 'N';

    F77_Name(strtrs)(&_upLo, &_trans, &_diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

int
trtrs(StorageUpLo upLo, Transpose trans, UnitDiag diag, int n, int nrhs,
      const double *a, int lda, double *b, int ldb)
{
    int info;
    char _upLo = (upLo==Upper) ? 'U' : 'L';
    char _trans = (trans==NoTrans) ? 'N' : 'T';
    char _diag = (diag==Unit) ? 'U' : 'N';

    F77_Name(dtrtrs)(&_upLo, &_trans, &_diag, &n, &nrhs, a, &lda, b, &ldb, &info);
    return info;
}

int
geqrf(int m, int n, float *a, int lda, float *tau, float *work, int lwork)
{
    int info;
    F77_Name(sgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

int
geqrf(int m, int n, double *a, int lda, double *tau, double *work, int lwork)
{
    int info;
    F77_Name(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
    return info;
}

int
orgqr(int m, int n, int k, float *a, int lda, const float *tau,
      float *work, int lwork)
{
    int info;
    F77_Name(sorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

int
orgqr(int m, int n, int k, double *a, int lda, const double *tau,
      double *work, int lwork)
{
    int info;
    F77_Name(dorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
    return info;
}

int
ormqr(BlasSide side, Transpose trans, int m, int n, int k,
      const float *a, int lda, const float *tau, float *c, int ldc,
      float *work, int lwork)
{
    int info;
    char _side = (side==Left) ? 'L' : 'R';
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    F77_Name(sormqr)(&_side, &_trans, &m, &n, &k, a, &lda, tau,
            c, &ldc, work, &lwork, &info);
    return info;
}

int
ormqr(BlasSide side, Transpose trans, int m, int n, int k,
      const double *a, int lda, const double *tau, double *c, int ldc,
      double *work, int lwork)
{
    int info;
    char _side = (side==Left) ? 'L' : 'R';
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    F77_Name(dormqr)(&_side, &_trans, &m, &n, &k, a, &lda, tau,
            c, &ldc, work, &lwork, &info);
    return info;
}

int
gels(Transpose trans, int m, int n, int nrhs, float *a, int lda,
     float *b, int ldb, float *work, int lwork)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    F77_Name(sgels)(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

int
gels(Transpose trans, int m, int n, int nrhs, double *a, int lda,
     double *b, int ldb, double *work, int lwork)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    F77_Name(dgels)(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

int
gels(Transpose trans, int m, int n, int nrhs, std::complex<double> *a, int lda,
     std::complex<double> *b, int ldb, std::complex<double> *work, int lwork)
{
    int info;
    char _trans = (trans==NoTrans) ? 'N' : 'T';

    F77_Name(zgels)(&_trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
    return info;
}

int
gelss(int m, int n, int nrhs, float *a, int lda, float *b, int ldb,
     float *s, float rcond, int rank, float *work, int lwork)
{
    int info;
    F77_Name(sgelss)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

int
gelss(int m, int n, int nrhs, double *a, int lda, double *b, int ldb,
     double *s, double rcond, int rank, double *work, int lwork)
{
    int info;
    F77_Name(dgelss)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

int
gelss(int m, int n, int nrhs, std::complex<double> *a, int lda,
      std::complex<double> *b, int ldb, std::complex<double> *s,
      std::complex<double> rcond, int rank, std::complex<double> *work,
      int lwork)
{
    int info;
    F77_Name(zgelss)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank,
            work, &lwork, &info);
    return info;
}

int
geev(bool jobvl, bool jobvr, int n, float *a, int lda,
     float *wr, float *wi,
     float *vl, int ldvl,
     float *vr, int ldvr,
     float *work, int lwork)
{
    int info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    F77_Name(sgeev)(&_jobvl, &_jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
           work, &lwork, &info);
    return info;
}

int
geev(bool jobvl, bool jobvr, int n, double *a, int lda,
     double *wr, double *wi,
     double *vl, int ldvl,
     double *vr, int ldvr,
     double *work, int lwork)
{
    int info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    F77_Name(dgeev)(&_jobvl, &_jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    return info;
}

int
geev(bool jobvl, bool jobvr, int n, std::complex<float> *a, int lda,
     std::complex<float> *w,
     std::complex<float> *vl, int ldvl,
     std::complex<float> *vr, int ldvr,
     std::complex<float> *work, int lwork, float *rwork)
{
    int info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    F77_Name(cgeev)(&_jobvl, &_jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
           work, &lwork, rwork, &info);
    return info;
}

int
geev(bool jobvl, bool jobvr, int n, std::complex<double> *a, int lda,
     std::complex<double> *w,
     std::complex<double> *vl, int ldvl,
     std::complex<double> *vr, int ldvr,
     std::complex<double> *work, int lwork, double *rwork)
{
    int info;
    char _jobvl = (jobvl==true) ? 'V' : 'N';
    char _jobvr = (jobvr==true) ? 'V' : 'N';

    F77_Name(zgeev)(&_jobvl, &_jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
           work, &lwork, rwork, &info);
    return info;
}

int
gesvd(char jobu, char jobvt,
      int m, int n, float *a, int lda,           // A
      float *s,                                  // singular values
      float *u, int ldu,                         // left singular vectors
      float *vt, int ldvt,                       // right singular vectors
      float *work, int lwork)
{
    int info;
    F77_Name(sgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, &info);
    return info;
}

int
gesvd(char jobu, char jobvt,
      int m, int n, double *a, int lda,           // A
      double *s,                                  // singular values
      double *u, int ldu,                         // left singular vectors
      double *vt, int ldvt,                       // right singular vectors
      double *work, int lwork)
{
    int info;
    F77_Name(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
            work, &lwork, &info);
    return info;
}

int
gecon(char norm,
      int n, double *a, int lda,                  // A
      double *anorm,                              // the norm of A
      double *rcond,                              // reciprocal condition number
      double *work, int *iwork)                  
{
    int info;
    F77_Name(dgecon)(&norm, &n, a, &lda, anorm, rcond, work, iwork, &info);
    return info;
}


} // namespace flens
