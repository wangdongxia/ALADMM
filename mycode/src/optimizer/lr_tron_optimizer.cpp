#include <cmath>
#include <cstring>
#include <algorithm>

#include "logging/simple_logging.h"
#include "math/simple_algebra.h"
#include "optimizer/lr_tron_optimizer.h"

int daxpy_(int *n, double *sa, double *sx, int *incx, double *sy, int *incy) {
    long int i, m, ix, iy, nn, iincx, iincy;
    register double ssa;

    /* constant times a vector plus a vector.
       uses unrolled loop for increments equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*) */

    /* Dereference inputs */
    nn = *n;
    ssa = *sa;
    iincx = *incx;
    iincy = *incy;

    if (nn > 0 && ssa != 0.0) {
        if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
        {
            m = nn - 3;
            for (i = 0; i < m; i += 4) {
                sy[i] += ssa * sx[i];
                sy[i + 1] += ssa * sx[i + 1];
                sy[i + 2] += ssa * sx[i + 2];
                sy[i + 3] += ssa * sx[i + 3];
            }
            for (; i < nn; ++i) /* clean-up loop */
                sy[i] += ssa * sx[i];
        } else /* code for unequal increments or equal increments not equal to 1 */
        {
            ix = iincx >= 0 ? 0 : (1 - nn) * iincx;
            iy = iincy >= 0 ? 0 : (1 - nn) * iincy;
            for (i = 0; i < nn; i++) {
                sy[iy] += ssa * sx[ix];
                ix += iincx;
                iy += iincy;
            }
        }
    }

    return 0;
} /* daxpy_ */

double ddot_(int *n, double *sx, int *incx, double *sy, int *incy) {
    long int i, m, nn, iincx, iincy;
    double stemp;
    long int ix, iy;

    /* forms the dot product of two vectors.
       uses unrolled loops for increments equal to one.
       jack dongarra, linpack, 3/11/78.
       modified 12/3/93, array(1) declarations changed to array(*) */

    /* Dereference inputs */
    nn = *n;
    iincx = *incx;
    iincy = *incy;

    stemp = 0.0;
    if (nn > 0) {
        if (iincx == 1 && iincy == 1) /* code for both increments equal to 1 */
        {
            m = nn - 4;
            for (i = 0; i < m; i += 5)
                stemp += sx[i] * sy[i] + sx[i + 1] * sy[i + 1] + sx[i + 2] * sy[i + 2] +
                         sx[i + 3] * sy[i + 3] + sx[i + 4] * sy[i + 4];

            for (; i < nn; i++)        /* clean-up loop */
                stemp += sx[i] * sy[i];
        } else /* code for unequal increments or equal increments not equal to 1 */
        {
            ix = 0;
            iy = 0;
            if (iincx < 0)
                ix = (1 - nn) * iincx;
            if (iincy < 0)
                iy = (1 - nn) * iincy;
            for (i = 0; i < nn; i++) {
                stemp += sx[ix] * sy[iy];
                ix += iincx;
                iy += iincy;
            }
        }
    }

    return stemp;
} /* ddot_ */

double dnrm2_(int *n, double *x, int *incx) {
    long int ix, nn, iincx;
    double norm, scale, absxi, ssq, temp;

/*  DNRM2 returns the euclidean norm of a vector via the function
    name, so that

       DNRM2 := sqrt( x'*x )

    -- This version written on 25-October-1982.
       Modified on 14-October-1993 to inline the call to SLASSQ.
       Sven Hammarling, Nag Ltd.   */

    /* Dereference inputs */
    nn = *n;
    iincx = *incx;

    if (nn > 0 && iincx > 0) {
        if (nn == 1) {
            norm = fabs(x[0]);
        } else {
            scale = 0.0;
            ssq = 1.0;

            /* The following loop is equivalent to this call to the LAPACK
               auxiliary routine:   CALL SLASSQ( N, X, INCX, SCALE, SSQ ) */

            for (ix = (nn - 1) * iincx; ix >= 0; ix -= iincx) {
                if (x[ix] != 0.0) {
                    absxi = fabs(x[ix]);
                    if (scale < absxi) {
                        temp = scale / absxi;
                        ssq = ssq * (temp * temp) + 1.0;
                        scale = absxi;
                    } else {
                        temp = absxi / scale;
                        ssq += temp * temp;
                    }
                }
            }
            norm = scale * sqrt(ssq);
        }
    } else
        norm = 0.0;

    return norm;

} /* dnrm2_ */

int dscal_(int *n, double *sa, double *sx, int *incx) {
    long int i, m, nincx, nn, iincx;
    double ssa;

    /* scales a vector by a constant.
       uses unrolled loops for increment equal to 1.
       jack dongarra, linpack, 3/11/78.
       modified 3/93 to return if incx .le. 0.
       modified 12/3/93, array(1) declarations changed to array(*) */

    /* Dereference inputs */
    nn = *n;
    iincx = *incx;
    ssa = *sa;

    if (nn > 0 && iincx > 0) {
        if (iincx == 1) /* code for increment equal to 1 */
        {
            m = nn - 4;
            for (i = 0; i < m; i += 5) {
                sx[i] = ssa * sx[i];
                sx[i + 1] = ssa * sx[i + 1];
                sx[i + 2] = ssa * sx[i + 2];
                sx[i + 3] = ssa * sx[i + 3];
                sx[i + 4] = ssa * sx[i + 4];
            }
            for (; i < nn; ++i) /* clean-up loop */
                sx[i] = ssa * sx[i];
        } else /* code for increment not equal to 1 */
        {
            nincx = nn * iincx;
            for (i = 0; i < nincx; i += iincx)
                sx[i] = ssa * sx[i];
        }
    }

    return 0;
} /* dscal_ */

double uTMv(int n, double *u, double *M, double *v) {
    const int m = n - 4;
    double res = 0;
    int i;
    for (i = 0; i < m; i += 5)
        res += u[i] * M[i] * v[i] + u[i + 1] * M[i + 1] * v[i + 1] + u[i + 2] * M[i + 2] * v[i + 2] +
               u[i + 3] * M[i + 3] * v[i + 3] + u[i + 4] * M[i + 4] * v[i + 4];
    for (; i < n; i++)
        res += u[i] * M[i] * v[i];
    return res;
}

LRTronOptimizer::LRTronOptimizer(const double *y, const double *z, int dimension, double rho, int max_iterations,
                                 double epsilon, double cg_epsilon, SparseDataset *dataset) :
        y_(y), z_(z), dimension_(dimension), rho_(rho),
        max_iterations_(max_iterations), epsilon_(epsilon),
        cg_epsilon_(cg_epsilon), dataset_(dataset) {
    int sample_num = dataset_->GetSampleNumber();
    D_ = new double[sample_num];
    int pos = 0;
    int neg = 0;
    for (int i = 0; i < sample_num; i++)
        if (dataset_->GetLabel(i) > 0)
            ++pos;
    neg = sample_num - pos;
    tron_epsilon_ = epsilon_ * std::max(std::min(pos, neg), 1) / sample_num;
}

LRTronOptimizer::~LRTronOptimizer() {
    delete[] D_;
}

void LRTronOptimizer::Optimize(double *x) {
    // Parameters for updating the iterates.
    double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;
    // Parameters for updating the trust region size delta.
    double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;

    int i;
    double delta = 0, sMnorm, one = 1.0;
    double alpha, f, fnew, prered, actred, gs;
    int search = 1, iter = 1, inc = 1;
    double *s = new double[dimension_];
    double *r = new double[dimension_];
    double *g = new double[dimension_];

    const double alpha_pcg = 0.01;
    double *M = new double[dimension_];

    // calculate gradient norm at x=0 for stopping condition.
    double *x0 = new double[dimension_];
    for (i = 0; i < dimension_; i++) {
        x0[i] = 0;
    }
    function_value(x0);
    gradient(x0, g);
    double gnorm0 = dnrm2_(&dimension_, g, &inc);
    delete[] x0;

    f = function_value(x);
    gradient(x, g);
    double gnorm = dnrm2_(&dimension_, g, &inc);

    if (gnorm <= tron_epsilon_ * gnorm0)
        search = 0;

    get_diag_preconditioner(M);
    for (i = 0; i < dimension_; i++)
        M[i] = (1 - alpha_pcg) + alpha_pcg * M[i];
    delta = sqrt(uTMv(dimension_, g, M, g));

    double *x_new = new double[dimension_];
    bool reach_boundary;
    bool delta_adjusted = false;
    while (iter <= max_iterations_ && search) {
        trpcg(delta, g, M, s, r, &reach_boundary);

        memcpy(x_new, x, sizeof(double) * dimension_);
        daxpy_(&dimension_, &one, s, &inc, x_new, &inc);

        gs = ddot_(&dimension_, g, &inc, s, &inc);
        prered = -0.5 * (gs - ddot_(&dimension_, s, &inc, r, &inc));
        fnew = function_value(x_new);

        // Compute the actual reduction.
        actred = f - fnew;

        // On the first iteration, adjust the initial step bound.
        sMnorm = sqrt(uTMv(dimension_, s, M, s));
        if (iter == 1 && !delta_adjusted) {
            delta = std::min(delta, sMnorm);
            delta_adjusted = true;
        }

        // Compute prediction alpha*sMnorm of the step.
        if (fnew - f - gs <= 0)
            alpha = sigma3;
        else
            alpha = std::max(sigma1, -0.5 * (gs / (fnew - f - gs)));

        // Update the trust region bound according to the ratio of actual to predicted reduction.
        if (actred < eta0 * prered)
            delta = std::min(alpha * sMnorm, sigma2 * delta);
        else if (actred < eta1 * prered)
            delta = std::max(sigma1 * delta, std::min(alpha * sMnorm, sigma2 * delta));
        else if (actred < eta2 * prered)
            delta = std::max(sigma1 * delta, std::min(alpha * sMnorm, sigma3 * delta));
        else {
            if (reach_boundary)
                delta = sigma3 * delta;
            else
                delta = std::max(delta, std::min(alpha * sMnorm, sigma3 * delta));
        }

        if (actred > eta0 * prered) {
            ++iter;
            memcpy(x, x_new, sizeof(double) * dimension_);
            f = fnew;
            gradient(x, g);
            get_diag_preconditioner(M);
            for (i = 0; i < dimension_; i++)
                M[i] = (1 - alpha_pcg) + alpha_pcg * M[i];

            gnorm = dnrm2_(&dimension_, g, &inc);
            if (gnorm <= tron_epsilon_ * gnorm0)
                break;
        }
        if (f < -1.0e+32) {
            LOG(WARNING) << "WARNING: f < -1.0e+32";
            break;
        }
        if (prered <= 0) {
            LOG(WARNING) << "WARNING: prered <= 0";
            break;
        }
        if (fabs(actred) <= 1.0e-12 * fabs(f) && fabs(prered) <= 1.0e-12 * fabs(f)) {
            LOG(WARNING) << "WARNING: actred and prered too small";
            break;
        }
    }

    delete[] g;
    delete[] r;
    delete[] x_new;
    delete[] s;
    delete[] M;
}

int LRTronOptimizer::trpcg(double delta, double *g, double *M, double *s, double *r, bool *reach_boundary) {
    int i, inc = 1;
    double one = 1;
    double *d = new double[dimension_];
    double *Hd = new double[dimension_];
    double zTr, znewTrnew, alpha, beta, cgtol;
    double *z = new double[dimension_];

    *reach_boundary = false;
    for (i = 0; i < dimension_; i++) {
        s[i] = 0;
        r[i] = -g[i];
        z[i] = r[i] / M[i];
        d[i] = z[i];
    }

    zTr = ddot_(&dimension_, z, &inc, r, &inc);
    cgtol = cg_epsilon_ * sqrt(zTr);
    int cg_iter = 0;
    int max_cg_iter = std::max(dimension_, 5);

    while (cg_iter < max_cg_iter) {
        if (sqrt(zTr) <= cgtol)
            break;
        cg_iter++;
        Hv(d, Hd);

        alpha = zTr / ddot_(&dimension_, d, &inc, Hd, &inc);
        daxpy_(&dimension_, &alpha, d, &inc, s, &inc);

        double sMnorm = sqrt(uTMv(dimension_, s, M, s));
        if (sMnorm > delta) {
            //std::cerr << "cg reaches trust region boundary" << std::endl;
            *reach_boundary = true;
            alpha = -alpha;
            daxpy_(&dimension_, &alpha, d, &inc, s, &inc);

            double sTMd = uTMv(dimension_, s, M, d);
            double sTMs = uTMv(dimension_, s, M, s);
            double dTMd = uTMv(dimension_, d, M, d);
            double dsq = delta * delta;
            double rad = sqrt(sTMd * sTMd + dTMd * (dsq - sTMs));
            if (sTMd >= 0)
                alpha = (dsq - sTMs) / (sTMd + rad);
            else
                alpha = (rad - sTMd) / dTMd;
            daxpy_(&dimension_, &alpha, d, &inc, s, &inc);
            alpha = -alpha;
            daxpy_(&dimension_, &alpha, Hd, &inc, r, &inc);
            break;
        }
        alpha = -alpha;
        daxpy_(&dimension_, &alpha, Hd, &inc, r, &inc);

        for (i = 0; i < dimension_; i++)
            z[i] = r[i] / M[i];
        znewTrnew = ddot_(&dimension_, z, &inc, r, &inc);
        beta = znewTrnew / zTr;
        dscal_(&dimension_, &beta, d, &inc);
        daxpy_(&dimension_, &one, z, &inc, d, &inc);
        zTr = znewTrnew;
    }

    if (cg_iter == max_cg_iter)
        LOG(WARNING) << "WARNING: reaching maximal number of CG steps";

    delete[] d;
    delete[] Hd;
    delete[] z;

    return cg_iter;
}

double LRTronOptimizer::function_value(const double *x) {
    double sum = 0.0;
    for (int i = 0; i < dimension_; ++i) {
        double temp = x[i] - z_[i] + y_[i] / rho_;
        sum += temp * temp;
    }
    sum *= (rho_ / 2.0);
    int sample_num = dataset_->GetSampleNumber();
    for (int i = 0; i < sample_num; ++i) {
        sum += std::log(1 + exp(-dataset_->GetLabel(i) * Dot(x, dataset_->GetSample(i))));
    }
    return sum;

}

void LRTronOptimizer::gradient(const double *x, double *g) {
    for (int i = 0; i < dimension_; ++i) {
        g[i] = y_[i] + rho_ * (x[i] - z_[i]);
    }
    int sample_num = dataset_->GetSampleNumber();
    for (int i = 0; i < sample_num; ++i) {
        double temp = Sigmoid(dataset_->GetLabel(i) * Dot(x, dataset_->GetSample(i)));
        D_[i] = temp * (1 - temp);
        temp = (temp - 1) * dataset_->GetLabel(i);
        const Feature *sample = dataset_->GetSample(i);
        while (sample->index != -1) {
            g[sample->index] += (sample->value * temp);
            ++sample;
        }
    }
}

void LRTronOptimizer::get_diag_preconditioner(double *M) {
    int sample_num = dataset_->GetSampleNumber();
    for (int i = 0; i < dimension_; i++) {
        M[i] = 1;
    }
    for (int i = 0; i < sample_num; i++) {
        const Feature *sample = dataset_->GetSample(i);
        while (sample->index != -1) {
            M[sample->index] += sample->value * sample->value * D_[i];
            ++sample;
        }
    }
}

void LRTronOptimizer::Hv(double *s, double *Hs) {
    int sample_num = dataset_->GetSampleNumber();
    for (int i = 0; i < dimension_; i++) {
        Hs[i] = 0;
    }
    for (int i = 0; i < sample_num; i++) {
        const Feature *sample = dataset_->GetSample(i);
        double xTs = Dot(s, dataset_->GetSample(i));
        xTs = D_[i] * xTs;
        while (sample->index != -1) {
            Hs[sample->index] += sample->value * xTs;
            ++sample;
        }
    }
    for (int i = 0; i < dimension_; i++)
        Hs[i] += rho_ * s[i];
}

void LRTronOptimizer::SetDataset(SparseDataset *dataset) {
    dataset_ = dataset;
    int sample_num = dataset_->GetSampleNumber();
    delete[] D_;
    D_ = new double[sample_num];
    int pos = 0;
    int neg = 0;
    for (int i = 0; i < sample_num; i++)
        if (dataset_->GetLabel(i) > 0)
            ++pos;
    neg = sample_num - pos;
    tron_epsilon_ = epsilon_ * std::max(std::min(pos, neg), 1) / sample_num;
}
