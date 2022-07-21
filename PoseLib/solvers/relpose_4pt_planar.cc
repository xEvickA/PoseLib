#include "relpose_4pt_planar.h"
#include <Eigen/Dense>
#include <math.h>
#include <stdio.h>
#include <PoseLib/misc/essential.h>

#define RELERROR      1.0e-12   /* smallest relative error we want */
//#define MAXPOW        0        /* max power of 10 we wish to search to */
#define MAXIT         800       /* max number of iterations */
#define SMALL_ENOUGH  1.0e-12   /* a coefficient smaller than SMALL_ENOUGH 
* is considered to be zero (0.0). */
#ifndef MAX_DEG
#define MAX_DEG      64
#endif

/* structure type for representing a polynomial */
typedef struct p {
	int ord;
	double coef[MAX_DEG + 1];
} poly;

/*---------------------------------------------------------------------------
* evalpoly
*
* evaluate polynomial defined in coef returning its value.
*--------------------------------------------------------------------------*/

double evalpoly(int ord, double *coef, double x)
{
	double *fp = &coef[ord];
	double f = *fp;

	for (fp--; fp >= coef; fp--)
		f = x * f + *fp;

	return(f);
}

int modrf_pos(int ord, double *coef, double a, double b,
	double *val, int invert)
{
	int  its;
	double fx, lfx;
	double *fp;
	double *scoef = coef;
	double *ecoef = &coef[ord];
	double fa, fb;

	// Invert the interval if required
	if (invert)
	{
		double temp = a;
		a = 1.0 / b;
		b = 1.0 / temp;
	}

	// Evaluate the polynomial at the end points
	if (invert)
	{
		fb = fa = *scoef;
		for (fp = scoef + 1; fp <= ecoef; fp++)
		{
			fa = a * fa + *fp;
			fb = b * fb + *fp;
		}
	}
	else
	{
		fb = fa = *ecoef;
		for (fp = ecoef - 1; fp >= scoef; fp--)
		{
			fa = a * fa + *fp;
			fb = b * fb + *fp;
		}
	}

	// if there is no sign difference the method won't work
	if (fa * fb > 0.0)
		return(0);

	// Return if the values are close to zero already
	if (fabs(fa) < RELERROR)
	{
		*val = invert ? 1.0 / a : a;
		return(1);
	}

	if (fabs(fb) < RELERROR)
	{
		*val = invert ? 1.0 / b : b;
		return(1);
	}

	lfx = fa;

	for (its = 0; its < MAXIT; its++)
	{
		// Assuming straight line from a to b, find zero
		double x = (fb * a - fa * b) / (fb - fa);

		// Evaluate the polynomial at x
		if (invert)
		{
			fx = *scoef;
			for (fp = scoef + 1; fp <= ecoef; fp++)
				fx = x * fx + *fp;
		}
		else
		{
			fx = *ecoef;
			for (fp = ecoef - 1; fp >= scoef; fp--)
				fx = x * fx + *fp;
		}

		// Evaluate two stopping conditions
		if (fabs(x) > RELERROR && fabs(fx / x) < RELERROR)
		{
			*val = invert ? 1.0 / x : x;
			return(1);
		}
		else if (fabs(fx) < RELERROR)
		{
			*val = invert ? 1.0 / x : x;
			return(1);
		}

		// Subdivide region, depending on whether fx has same sign as fa or fb
		if ((fa * fx) < 0)
		{
			b = x;
			fb = fx;
			if ((lfx * fx) > 0)
				fa /= 2;
		}
		else
		{
			a = x;
			fa = fx;
			if ((lfx * fx) > 0)
				fb /= 2;
		}


		// Return if the difference between a and b is very small
		if (fabs(b - a) < fabs(RELERROR * a))
		{
			*val = invert ? 1.0 / a : a;
			return(1);
		}

		lfx = fx;
	}

	//==================================================================
	// This is debugging in case something goes wrong.
	// If we reach here, we have not converged -- give some diagnostics
	//==================================================================

	fprintf(stderr, "modrf overflow on interval %f %f\n", a, b);
	fprintf(stderr, "\t b-a = %12.5e\n", b - a);
	fprintf(stderr, "\t fa  = %12.5e\n", fa);
	fprintf(stderr, "\t fb  = %12.5e\n", fb);
	fprintf(stderr, "\t fx  = %12.5e\n", fx);

	// Evaluate the true values at a and b
	if (invert)
	{
		fb = fa = *scoef;
		for (fp = scoef + 1; fp <= ecoef; fp++)
		{
			fa = a * fa + *fp;
			fb = b * fb + *fp;
		}
	}
	else
	{
		fb = fa = *ecoef;
		for (fp = ecoef - 1; fp >= scoef; fp--)
		{
			fa = a * fa + *fp;
			fb = b * fb + *fp;
		}
	}

	fprintf(stderr, "\t true fa = %12.5e\n", fa);
	fprintf(stderr, "\t true fb = %12.5e\n", fb);
	fprintf(stderr, "\t gradient= %12.5e\n", (fb - fa) / (b - a));

	// Print out the polynomial
	fprintf(stderr, "Polynomial coefficients\n");
	for (fp = ecoef; fp >= scoef; fp--)
		fprintf(stderr, "\t%12.5e\n", *fp);

	return(0);
}

/*---------------------------------------------------------------------------
* modrf
*
* uses the modified regula-falsi method to evaluate the root
* in interval [a,b] of the polynomial described in coef. The
* root is returned is returned in *val. The routine returns zero
* if it can't converge.
*--------------------------------------------------------------------------*/

int modrf(int ord, double *coef, double a, double b, double *val)
{
	// This is an interfact to modrf that takes account of different cases
	// The idea is that the basic routine works badly for polynomials on
	// intervals that extend well beyond [-1, 1], because numbers get too large

	double *fp;
	double *scoef = coef;
	double *ecoef = &coef[ord];
	const int invert = 1;

	double fp1 = 0.0, fm1 = 0.0; // Values of function at 1 and -1
	double fa = 0.0, fb = 0.0; // Values at end points

	// We assume that a < b
	if (a > b)
	{
		double temp = a;
		a = b;
		b = temp;
	}

	// The normal case, interval is inside [-1, 1]
	if (b <= 1.0 && a >= -1.0) return modrf_pos(ord, coef, a, b, val, !invert);

	// The case where the interval is outside [-1, 1]
	if (a >= 1.0 || b <= -1.0)
		return modrf_pos(ord, coef, a, b, val, invert);

	// If we have got here, then the interval includes the points 1 or -1.
	// In this case, we need to evaluate at these points

	// Evaluate the polynomial at the end points
	for (fp = ecoef - 1; fp >= scoef; fp--)
	{
		fp1 = *fp + fp1;
		fm1 = *fp - fm1;
		fa = a * fa + *fp;
		fb = b * fb + *fp;
	}

	// Then there is the case where the interval contains -1 or 1
	if (a < -1.0 && b > 1.0)
	{
		// Interval crosses over 1.0, so cut
		if (fa * fm1 < 0.0)      // The solution is between a and -1
			return modrf_pos(ord, coef, a, -1.0, val, invert);
		else if (fb * fp1 < 0.0) // The solution is between 1 and b
			return modrf_pos(ord, coef, 1.0, b, val, invert);
		else                     // The solution is between -1 and 1
			return modrf_pos(ord, coef, -1.0, 1.0, val, !invert);
	}
	else if (a < -1.0)
	{
		// Interval crosses over 1.0, so cut
		if (fa * fm1 < 0.0)      // The solution is between a and -1
			return modrf_pos(ord, coef, a, -1.0, val, invert);
		else                     // The solution is between -1 and b
			return modrf_pos(ord, coef, -1.0, b, val, !invert);
	}
	else  // b > 1.0
	{
		if (fb * fp1 < 0.0) // The solution is between 1 and b
			return modrf_pos(ord, coef, 1.0, b, val, invert);
		else                     // The solution is between a and 1
			return modrf_pos(ord, coef, a, 1.0, val, !invert);
	}
}

/*---------------------------------------------------------------------------
* modp
*
*  calculates the modulus of u(x) / v(x) leaving it in r, it
*  returns 0 if r(x) is a constant.
*  note: this function assumes the leading coefficient of v is 1 or -1
*--------------------------------------------------------------------------*/

static int modp(poly *u, poly *v, poly *r)
{
	int j, k;  /* Loop indices */

	double *nr = r->coef;
	double *end = &u->coef[u->ord];

	double *uc = u->coef;
	while (uc <= end)
		*nr++ = *uc++;

	if (v->coef[v->ord] < 0.0) {

		for (k = u->ord - v->ord - 1; k >= 0; k -= 2)
			r->coef[k] = -r->coef[k];

		for (k = u->ord - v->ord; k >= 0; k--)
			for (j = v->ord + k - 1; j >= k; j--)
				r->coef[j] = -r->coef[j] - r->coef[v->ord + k]
				* v->coef[j - k];
	}
	else {
		for (k = u->ord - v->ord; k >= 0; k--)
			for (j = v->ord + k - 1; j >= k; j--)
				r->coef[j] -= r->coef[v->ord + k] * v->coef[j - k];
	}

	k = v->ord - 1;
	while (k >= 0 && fabs(r->coef[k]) < SMALL_ENOUGH) {
		r->coef[k] = 0.0;
		k--;
	}

	r->ord = (k < 0) ? 0 : k;

	return(r->ord);
}

/*---------------------------------------------------------------------------
* buildsturm
*
* build up a sturm sequence for a polynomial in smat, returning
* the number of polynomials in the sequence
*--------------------------------------------------------------------------*/

int buildsturm(int ord, poly *sseq)
{
	sseq[0].ord = ord;
	sseq[1].ord = ord - 1;

	/* calculate the derivative and normalise the leading coefficient */
	{
		int i;    // Loop index
		poly *sp;
		double f = fabs(sseq[0].coef[ord] * ord);
		double *fp = sseq[1].coef;
		double *fc = sseq[0].coef + 1;

		for (i = 1; i <= ord; i++)
			*fp++ = *fc++ * i / f;

		/* construct the rest of the Sturm sequence */
		for (sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++) {

			/* reverse the sign and normalise */
			f = -fabs(sp->coef[sp->ord]);
			for (fp = &sp->coef[sp->ord]; fp >= sp->coef; fp--)
				*fp /= f;
		}

		sp->coef[0] = -sp->coef[0]; /* reverse the sign */

		return(sp - sseq);
	}
}

/*---------------------------------------------------------------------------
* numchanges
*
* return the number of sign changes in the Sturm sequence in
* sseq at the value a.
*--------------------------------------------------------------------------*/

int numchanges(int np, poly *sseq, double a)
{
	int changes = 0;

	double lf = evalpoly(sseq[0].ord, sseq[0].coef, a);

	poly *s;
	for (s = sseq + 1; s <= sseq + np; s++) {
		double f = evalpoly(s->ord, s->coef, a);
		if (lf == 0.0 || lf * f < 0)
			changes++;
		lf = f;
	}

	return(changes);
}

/*---------------------------------------------------------------------------
* numroots
*
* return the number of distinct real roots of the polynomial described in sseq.
*--------------------------------------------------------------------------*/

int numroots(int np, poly *sseq, int *atneg, int *atpos, bool non_neg)
{
	int atposinf = 0;
	int atneginf = 0;

	/* changes at positive infinity */
	double f;
	double lf = sseq[0].coef[sseq[0].ord];

	poly *s;
	for (s = sseq + 1; s <= sseq + np; s++) {
		f = s->coef[s->ord];
		if (lf == 0.0 || lf * f < 0)
			atposinf++;
		lf = f;
	}

	// changes at negative infinity or zero
	if (non_neg)
		atneginf = numchanges(np, sseq, 0.0);

	else
	{
		if (sseq[0].ord & 1)
			lf = -sseq[0].coef[sseq[0].ord];
		else
			lf = sseq[0].coef[sseq[0].ord];

		for (s = sseq + 1; s <= sseq + np; s++) {
			if (s->ord & 1)
				f = -s->coef[s->ord];
			else
				f = s->coef[s->ord];
			if (lf == 0.0 || lf * f < 0)
				atneginf++;
			lf = f;
		}
	}

	*atneg = atneginf;
	*atpos = atposinf;

	return(atneginf - atposinf);
}


/*---------------------------------------------------------------------------
* sbisect
*
* uses a bisection based on the sturm sequence for the polynomial
* described in sseq to isolate intervals in which roots occur,
* the roots are returned in the roots array in order of magnitude.
*--------------------------------------------------------------------------*/

int sbisect(int np, poly *sseq,
	double min, double max,
	int atmin, int atmax,
	double *roots)
{
	double mid;
	int atmid;
	int its;
	int  n1 = 0, n2 = 0;
	int nroot = atmin - atmax;

	if (nroot == 1) {

		/* first try a less expensive technique.  */
		if (modrf(sseq->ord, sseq->coef, min, max, &roots[0]))
			return 1;

		/*
		* if we get here we have to evaluate the root the hard
		* way by using the Sturm sequence.
		*/
		for (its = 0; its < MAXIT; its++) {
			mid = (double)((min + max) / 2);
			atmid = numchanges(np, sseq, mid);

			if (fabs(mid) > RELERROR) {
				if (fabs((max - min) / mid) < RELERROR) {
					roots[0] = mid;
					return 1;
				}
			}
			else if (fabs(max - min) < RELERROR) {
				roots[0] = mid;
				return 1;
			}

			if ((atmin - atmid) == 0)
				min = mid;
			else
				max = mid;
		}

		if (its == MAXIT) {
			fprintf(stderr, "sbisect: overflow min %f max %f\
							                         diff %e nroot %d n1 %d n2 %d\n",
													 min, max, max - min, nroot, n1, n2);
			roots[0] = mid;
		}

		return 1;
	}

	/* more than one root in the interval, we have to bisect */
	for (its = 0; its < MAXIT; its++) {

		mid = (double)((min + max) / 2);
		atmid = numchanges(np, sseq, mid);

		n1 = atmin - atmid;
		n2 = atmid - atmax;

		if (n1 != 0 && n2 != 0) {
			sbisect(np, sseq, min, mid, atmin, atmid, roots);
			sbisect(np, sseq, mid, max, atmid, atmax, &roots[n1]);
			break;
		}

		if (n1 == 0)
			min = mid;
		else
			max = mid;
	}

	if (its == MAXIT) {
		fprintf(stderr, "sbisect: roots too close together\n");
		fprintf(stderr, "sbisect: overflow min %f max %f diff %e\
						                      nroot %d n1 %d n2 %d\n",
											  min, max, max - min, nroot, n1, n2);
		for (n1 = atmax; n1 < atmin; n1++)
			roots[n1 - atmax] = mid;
	}

	return 1;
}

int find_real_roots_sturm(
	double *p, int order, double *roots, int *nroots, int maxpow, bool non_neg)
{
	/*
	* finds the roots of the input polynomial.  They are returned in roots.
	* It is assumed that roots is already allocated with space for the roots.
	*/

	poly sseq[MAX_DEG + 1];
	double  min, max;
	int  i, nchanges, np, atmin, atmax;

	// Copy the coefficients from the input p.  Normalize as we go
	double norm = 1.0 / p[order];
	for (i = 0; i <= order; i++)
		sseq[0].coef[i] = p[i] * norm;

	// Now, also normalize the other terms
	double val0 = fabs(sseq[0].coef[0]);
	double fac = 1.0; // This will be a factor for the roots
	if (val0 > 10.0)  // Do this in case there are zero roots
	{
		fac = pow(val0, -1.0 / order);
		double mult = fac;
		for (int i = order - 1; i >= 0; i--)
		{
			sseq[0].coef[i] *= mult;
			mult = mult * fac;
		}
	}

	/* build the Sturm sequence */
	np = buildsturm(order, sseq);

#ifdef RH_DEBUG
	{
		int i, j;

		printf("Sturm sequence for:\n");
		for (i = order; i >= 0; i--)
			printf("%lf ", sseq[0].coef[i]);
		printf("\n\n");

		for (i = 0; i <= np; i++) {
			for (j = sseq[i].ord; j >= 0; j--)
				printf("%10f ", sseq[i].coef[j]);
			printf("\n");
		}

		printf("\n");
	}
#endif

	// get the number of real roots
	*nroots = numroots(np, sseq, &atmin, &atmax, non_neg);

	if (*nroots == 0) {
		// fprintf(stderr, "solve: no real roots\n");
		return 0;
	}

	/* calculate the bracket that the roots live in */
	if (non_neg) min = 0.0;
	else
	{
		min = -1.0;
		nchanges = numchanges(np, sseq, min);
		for (i = 0; nchanges != atmin && i != maxpow; i++) {
			min *= 10.0;
			nchanges = numchanges(np, sseq, min);
		}

		if (nchanges != atmin) {
			//printf("solve: unable to bracket all negative roots\n");
			atmin = nchanges;
		}
	}

	max = 1.0;
	nchanges = numchanges(np, sseq, max);
	for (i = 0; nchanges != atmax && i != maxpow; i++) {
		max *= 10.0;
		nchanges = numchanges(np, sseq, max);
	}

	if (nchanges != atmax) {
		//printf("solve: unable to bracket all positive roots\n");
		atmax = nchanges;
	}

	*nroots = atmin - atmax;

	/* perform the bisection */
	sbisect(np, sseq, min, max, atmin, atmax, roots);

	/* Finally, reorder the roots */
	for (i = 0; i<*nroots; i++)
		roots[i] /= fac;

#ifdef RH_DEBUG

	/* write out the roots */
	printf("Number of roots = %d\n", *nroots);
	for (i = 0; i<*nroots; i++)
		printf("%12.5e\n", roots[i]);

#endif

	return 1;
}


template <typename Derived> void charpoly_danilevsky_piv(Eigen::MatrixBase<Derived> &A, double *p) {
    int n = A.rows();

    for (int i = n - 1; i > 0; i--) {

        int piv_ind = i - 1;
        double piv = std::abs(A(i, i - 1));

        // Find largest pivot
        for (int j = 0; j < i - 1; j++) {
            if (std::abs(A(i, j)) > piv) {
                piv = std::abs(A(i, j));
                piv_ind = j;
            }
        }
        if (piv_ind != i - 1) {
            // Perform permutation
            A.row(i - 1).swap(A.row(piv_ind));
            A.col(i - 1).swap(A.col(piv_ind));
        }
        piv = A(i, i - 1);

        Eigen::VectorXd v = A.row(i);
        A.row(i - 1) = v.transpose() * A;

        Eigen::VectorXd vinv = (-1.0) * v;
        vinv(i - 1) = 1;
        vinv /= piv;
        vinv(i - 1) -= 1;
        Eigen::VectorXd Acol = A.col(i - 1);
        for (int j = 0; j <= i; j++)
            A.row(j) = A.row(j) + Acol(j) * vinv.transpose();

        A.row(i) = Eigen::VectorXd::Zero(n);
        A(i, i - 1) = 1;
    }
    p[n] = 1;
    for (int i = 0; i < n; i++)
        p[i] = -A(0, n - i - 1);
}

void fast_eigenvector_solver(double *eigv, int neig, Eigen::Matrix<double, 30, 30> &AM,
                             Eigen::Matrix<std::complex<double>, 4, 30> &sols) {
    static const int ind[] = {2, 4, 6, 9, 12, 15, 18, 22, 25, 29};
    // Truncated action matrix containing non-trivial rows
    Eigen::Matrix<double, 10, 30> AMs;
    double zi[5];

    for (int i = 0; i < 10; i++) {
        AMs.row(i) = AM.row(ind[i]);
    }
    for (int i = 0; i < neig; i++) {
        zi[0] = eigv[i];
        for (int j = 1; j < 5; j++) {
            zi[j] = zi[j - 1] * eigv[i];
        }
        Eigen::Matrix<double, 10, 10> AA;
        AA.col(0) = AMs.col(2);
        AA.col(1) = AMs.col(3) + zi[0] * AMs.col(4);
        AA.col(2) = AMs.col(5) + zi[0] * AMs.col(6);
        AA.col(3) = AMs.col(1) + zi[0] * AMs.col(7) + zi[1] * AMs.col(8) + zi[2] * AMs.col(9);
        AA.col(4) = AMs.col(11) + zi[0] * AMs.col(12);
        AA.col(5) = AMs.col(13) + zi[0] * AMs.col(14) + zi[1] * AMs.col(15);
        AA.col(6) = AMs.col(10) + zi[0] * AMs.col(16) + zi[1] * AMs.col(17) + zi[2] * AMs.col(18);
        AA.col(7) = AMs.col(20) + zi[0] * AMs.col(21) + zi[1] * AMs.col(22);
        AA.col(8) = AMs.col(19) + zi[0] * AMs.col(23) + zi[1] * AMs.col(24) + zi[2] * AMs.col(25);
        AA.col(9) = AMs.col(0) + zi[0] * AMs.col(26) + zi[1] * AMs.col(27) + zi[2] * AMs.col(28) + zi[3] * AMs.col(29);
        AA(0, 0) = AA(0, 0) - zi[0];
        AA(1, 1) = AA(1, 1) - zi[1];
        AA(2, 2) = AA(2, 2) - zi[1];
        AA(3, 3) = AA(3, 3) - zi[3];
        AA(4, 4) = AA(4, 4) - zi[1];
        AA(5, 5) = AA(5, 5) - zi[2];
        AA(6, 6) = AA(6, 6) - zi[3];
        AA(7, 7) = AA(7, 7) - zi[2];
        AA(8, 8) = AA(8, 8) - zi[3];
        AA(9, 9) = AA(9, 9) - zi[4];

        Eigen::Matrix<double, 9, 1> s = AA.leftCols(9).colPivHouseholderQr().solve(-AA.col(9));
        sols(0, i) = s(3);
        sols(1, i) = s(6);
        sols(2, i) = s(8);
        sols(3, i) = zi[0];
    }
}

int get_sols_relpose_4pt_planar(const Eigen::VectorXd &data, Eigen::Matrix<std::complex<double>, 4, 30> &sols) {
    // Compute coefficients
    const double *d = data.data();
    Eigen::VectorXd coeffs(385);
    coeffs[0] = -d[2] * d[4] * d[6] + d[1] * d[5] * d[6] + d[2] * d[3] * d[7] - d[0] * d[5] * d[7] -
                d[1] * d[3] * d[8] + d[0] * d[4] * d[8];
    coeffs[1] = -d[5] * d[7] * d[9] + d[4] * d[8] * d[9] + d[5] * d[6] * d[10] - d[3] * d[8] * d[10] -
                d[4] * d[6] * d[11] + d[3] * d[7] * d[11] + d[2] * d[7] * d[12] - d[1] * d[8] * d[12] -
                d[2] * d[6] * d[13] + d[0] * d[8] * d[13] + d[1] * d[6] * d[14] - d[0] * d[7] * d[14] -
                d[2] * d[4] * d[15] + d[1] * d[5] * d[15] + d[2] * d[3] * d[16] - d[0] * d[5] * d[16] -
                d[1] * d[3] * d[17] + d[0] * d[4] * d[17];
    coeffs[2] = -d[8] * d[10] * d[12] + d[7] * d[11] * d[12] + d[8] * d[9] * d[13] - d[6] * d[11] * d[13] -
                d[7] * d[9] * d[14] + d[6] * d[10] * d[14] + d[5] * d[10] * d[15] - d[4] * d[11] * d[15] -
                d[2] * d[13] * d[15] + d[1] * d[14] * d[15] - d[5] * d[9] * d[16] + d[3] * d[11] * d[16] +
                d[2] * d[12] * d[16] - d[0] * d[14] * d[16] + d[4] * d[9] * d[17] - d[3] * d[10] * d[17] -
                d[1] * d[12] * d[17] + d[0] * d[13] * d[17];
    coeffs[3] = -d[11] * d[13] * d[15] + d[10] * d[14] * d[15] + d[11] * d[12] * d[16] - d[9] * d[14] * d[16] -
                d[10] * d[12] * d[17] + d[9] * d[13] * d[17];
    coeffs[4] = -d[5] * d[7] * d[18] + d[4] * d[8] * d[18] + d[5] * d[6] * d[19] - d[3] * d[8] * d[19] -
                d[4] * d[6] * d[20] + d[3] * d[7] * d[20] + d[2] * d[7] * d[21] - d[1] * d[8] * d[21] -
                d[2] * d[6] * d[22] + d[0] * d[8] * d[22] + d[1] * d[6] * d[23] - d[0] * d[7] * d[23] -
                d[2] * d[4] * d[24] + d[1] * d[5] * d[24] + d[2] * d[3] * d[25] - d[0] * d[5] * d[25] -
                d[1] * d[3] * d[26] + d[0] * d[4] * d[26];
    coeffs[5] = d[8] * d[13] * d[18] - d[7] * d[14] * d[18] - d[5] * d[16] * d[18] + d[4] * d[17] * d[18] -
                d[8] * d[12] * d[19] + d[6] * d[14] * d[19] + d[5] * d[15] * d[19] - d[3] * d[17] * d[19] +
                d[7] * d[12] * d[20] - d[6] * d[13] * d[20] - d[4] * d[15] * d[20] + d[3] * d[16] * d[20] -
                d[8] * d[10] * d[21] + d[7] * d[11] * d[21] + d[2] * d[16] * d[21] - d[1] * d[17] * d[21] +
                d[8] * d[9] * d[22] - d[6] * d[11] * d[22] - d[2] * d[15] * d[22] + d[0] * d[17] * d[22] -
                d[7] * d[9] * d[23] + d[6] * d[10] * d[23] + d[1] * d[15] * d[23] - d[0] * d[16] * d[23] +
                d[5] * d[10] * d[24] - d[4] * d[11] * d[24] - d[2] * d[13] * d[24] + d[1] * d[14] * d[24] -
                d[5] * d[9] * d[25] + d[3] * d[11] * d[25] + d[2] * d[12] * d[25] - d[0] * d[14] * d[25] +
                d[4] * d[9] * d[26] - d[3] * d[10] * d[26] - d[1] * d[12] * d[26] + d[0] * d[13] * d[26];
    coeffs[6] = -d[14] * d[16] * d[18] + d[13] * d[17] * d[18] + d[14] * d[15] * d[19] - d[12] * d[17] * d[19] -
                d[13] * d[15] * d[20] + d[12] * d[16] * d[20] + d[11] * d[16] * d[21] - d[10] * d[17] * d[21] -
                d[11] * d[15] * d[22] + d[9] * d[17] * d[22] + d[10] * d[15] * d[23] - d[9] * d[16] * d[23] -
                d[11] * d[13] * d[24] + d[10] * d[14] * d[24] + d[11] * d[12] * d[25] - d[9] * d[14] * d[25] -
                d[10] * d[12] * d[26] + d[9] * d[13] * d[26];
    coeffs[7] = -d[8] * d[19] * d[21] + d[7] * d[20] * d[21] + d[8] * d[18] * d[22] - d[6] * d[20] * d[22] -
                d[7] * d[18] * d[23] + d[6] * d[19] * d[23] + d[5] * d[19] * d[24] - d[4] * d[20] * d[24] -
                d[2] * d[22] * d[24] + d[1] * d[23] * d[24] - d[5] * d[18] * d[25] + d[3] * d[20] * d[25] +
                d[2] * d[21] * d[25] - d[0] * d[23] * d[25] + d[4] * d[18] * d[26] - d[3] * d[19] * d[26] -
                d[1] * d[21] * d[26] + d[0] * d[22] * d[26];
    coeffs[8] = -d[17] * d[19] * d[21] + d[16] * d[20] * d[21] + d[17] * d[18] * d[22] - d[15] * d[20] * d[22] -
                d[16] * d[18] * d[23] + d[15] * d[19] * d[23] + d[14] * d[19] * d[24] - d[13] * d[20] * d[24] -
                d[11] * d[22] * d[24] + d[10] * d[23] * d[24] - d[14] * d[18] * d[25] + d[12] * d[20] * d[25] +
                d[11] * d[21] * d[25] - d[9] * d[23] * d[25] + d[13] * d[18] * d[26] - d[12] * d[19] * d[26] -
                d[10] * d[21] * d[26] + d[9] * d[22] * d[26];
    coeffs[9] = -d[20] * d[22] * d[24] + d[19] * d[23] * d[24] + d[20] * d[21] * d[25] - d[18] * d[23] * d[25] -
                d[19] * d[21] * d[26] + d[18] * d[22] * d[26];
    coeffs[10] = -d[5] * d[7] * d[27] + d[4] * d[8] * d[27] + d[5] * d[6] * d[28] - d[3] * d[8] * d[28] -
                 d[4] * d[6] * d[29] + d[3] * d[7] * d[29] + d[2] * d[7] * d[30] - d[1] * d[8] * d[30] -
                 d[2] * d[6] * d[31] + d[0] * d[8] * d[31] + d[1] * d[6] * d[32] - d[0] * d[7] * d[32] -
                 d[2] * d[4] * d[33] + d[1] * d[5] * d[33] + d[2] * d[3] * d[34] - d[0] * d[5] * d[34] -
                 d[1] * d[3] * d[35] + d[0] * d[4] * d[35];
    coeffs[11] = d[8] * d[13] * d[27] - d[7] * d[14] * d[27] - d[5] * d[16] * d[27] + d[4] * d[17] * d[27] -
                 d[8] * d[12] * d[28] + d[6] * d[14] * d[28] + d[5] * d[15] * d[28] - d[3] * d[17] * d[28] +
                 d[7] * d[12] * d[29] - d[6] * d[13] * d[29] - d[4] * d[15] * d[29] + d[3] * d[16] * d[29] -
                 d[8] * d[10] * d[30] + d[7] * d[11] * d[30] + d[2] * d[16] * d[30] - d[1] * d[17] * d[30] +
                 d[8] * d[9] * d[31] - d[6] * d[11] * d[31] - d[2] * d[15] * d[31] + d[0] * d[17] * d[31] -
                 d[7] * d[9] * d[32] + d[6] * d[10] * d[32] + d[1] * d[15] * d[32] - d[0] * d[16] * d[32] +
                 d[5] * d[10] * d[33] - d[4] * d[11] * d[33] - d[2] * d[13] * d[33] + d[1] * d[14] * d[33] -
                 d[5] * d[9] * d[34] + d[3] * d[11] * d[34] + d[2] * d[12] * d[34] - d[0] * d[14] * d[34] +
                 d[4] * d[9] * d[35] - d[3] * d[10] * d[35] - d[1] * d[12] * d[35] + d[0] * d[13] * d[35];
    coeffs[12] = -d[14] * d[16] * d[27] + d[13] * d[17] * d[27] + d[14] * d[15] * d[28] - d[12] * d[17] * d[28] -
                 d[13] * d[15] * d[29] + d[12] * d[16] * d[29] + d[11] * d[16] * d[30] - d[10] * d[17] * d[30] -
                 d[11] * d[15] * d[31] + d[9] * d[17] * d[31] + d[10] * d[15] * d[32] - d[9] * d[16] * d[32] -
                 d[11] * d[13] * d[33] + d[10] * d[14] * d[33] + d[11] * d[12] * d[34] - d[9] * d[14] * d[34] -
                 d[10] * d[12] * d[35] + d[9] * d[13] * d[35];
    coeffs[13] = d[8] * d[22] * d[27] - d[7] * d[23] * d[27] - d[5] * d[25] * d[27] + d[4] * d[26] * d[27] -
                 d[8] * d[21] * d[28] + d[6] * d[23] * d[28] + d[5] * d[24] * d[28] - d[3] * d[26] * d[28] +
                 d[7] * d[21] * d[29] - d[6] * d[22] * d[29] - d[4] * d[24] * d[29] + d[3] * d[25] * d[29] -
                 d[8] * d[19] * d[30] + d[7] * d[20] * d[30] + d[2] * d[25] * d[30] - d[1] * d[26] * d[30] +
                 d[8] * d[18] * d[31] - d[6] * d[20] * d[31] - d[2] * d[24] * d[31] + d[0] * d[26] * d[31] -
                 d[7] * d[18] * d[32] + d[6] * d[19] * d[32] + d[1] * d[24] * d[32] - d[0] * d[25] * d[32] +
                 d[5] * d[19] * d[33] - d[4] * d[20] * d[33] - d[2] * d[22] * d[33] + d[1] * d[23] * d[33] -
                 d[5] * d[18] * d[34] + d[3] * d[20] * d[34] + d[2] * d[21] * d[34] - d[0] * d[23] * d[34] +
                 d[4] * d[18] * d[35] - d[3] * d[19] * d[35] - d[1] * d[21] * d[35] + d[0] * d[22] * d[35];
    coeffs[14] = d[17] * d[22] * d[27] - d[16] * d[23] * d[27] - d[14] * d[25] * d[27] + d[13] * d[26] * d[27] -
                 d[17] * d[21] * d[28] + d[15] * d[23] * d[28] + d[14] * d[24] * d[28] - d[12] * d[26] * d[28] +
                 d[16] * d[21] * d[29] - d[15] * d[22] * d[29] - d[13] * d[24] * d[29] + d[12] * d[25] * d[29] -
                 d[17] * d[19] * d[30] + d[16] * d[20] * d[30] + d[11] * d[25] * d[30] - d[10] * d[26] * d[30] +
                 d[17] * d[18] * d[31] - d[15] * d[20] * d[31] - d[11] * d[24] * d[31] + d[9] * d[26] * d[31] -
                 d[16] * d[18] * d[32] + d[15] * d[19] * d[32] + d[10] * d[24] * d[32] - d[9] * d[25] * d[32] +
                 d[14] * d[19] * d[33] - d[13] * d[20] * d[33] - d[11] * d[22] * d[33] + d[10] * d[23] * d[33] -
                 d[14] * d[18] * d[34] + d[12] * d[20] * d[34] + d[11] * d[21] * d[34] - d[9] * d[23] * d[34] +
                 d[13] * d[18] * d[35] - d[12] * d[19] * d[35] - d[10] * d[21] * d[35] + d[9] * d[22] * d[35];
    coeffs[15] = -d[23] * d[25] * d[27] + d[22] * d[26] * d[27] + d[23] * d[24] * d[28] - d[21] * d[26] * d[28] -
                 d[22] * d[24] * d[29] + d[21] * d[25] * d[29] + d[20] * d[25] * d[30] - d[19] * d[26] * d[30] -
                 d[20] * d[24] * d[31] + d[18] * d[26] * d[31] + d[19] * d[24] * d[32] - d[18] * d[25] * d[32] -
                 d[20] * d[22] * d[33] + d[19] * d[23] * d[33] + d[20] * d[21] * d[34] - d[18] * d[23] * d[34] -
                 d[19] * d[21] * d[35] + d[18] * d[22] * d[35];
    coeffs[16] = -d[8] * d[28] * d[30] + d[7] * d[29] * d[30] + d[8] * d[27] * d[31] - d[6] * d[29] * d[31] -
                 d[7] * d[27] * d[32] + d[6] * d[28] * d[32] + d[5] * d[28] * d[33] - d[4] * d[29] * d[33] -
                 d[2] * d[31] * d[33] + d[1] * d[32] * d[33] - d[5] * d[27] * d[34] + d[3] * d[29] * d[34] +
                 d[2] * d[30] * d[34] - d[0] * d[32] * d[34] + d[4] * d[27] * d[35] - d[3] * d[28] * d[35] -
                 d[1] * d[30] * d[35] + d[0] * d[31] * d[35];
    coeffs[17] = -d[17] * d[28] * d[30] + d[16] * d[29] * d[30] + d[17] * d[27] * d[31] - d[15] * d[29] * d[31] -
                 d[16] * d[27] * d[32] + d[15] * d[28] * d[32] + d[14] * d[28] * d[33] - d[13] * d[29] * d[33] -
                 d[11] * d[31] * d[33] + d[10] * d[32] * d[33] - d[14] * d[27] * d[34] + d[12] * d[29] * d[34] +
                 d[11] * d[30] * d[34] - d[9] * d[32] * d[34] + d[13] * d[27] * d[35] - d[12] * d[28] * d[35] -
                 d[10] * d[30] * d[35] + d[9] * d[31] * d[35];
    coeffs[18] = -d[26] * d[28] * d[30] + d[25] * d[29] * d[30] + d[26] * d[27] * d[31] - d[24] * d[29] * d[31] -
                 d[25] * d[27] * d[32] + d[24] * d[28] * d[32] + d[23] * d[28] * d[33] - d[22] * d[29] * d[33] -
                 d[20] * d[31] * d[33] + d[19] * d[32] * d[33] - d[23] * d[27] * d[34] + d[21] * d[29] * d[34] +
                 d[20] * d[30] * d[34] - d[18] * d[32] * d[34] + d[22] * d[27] * d[35] - d[21] * d[28] * d[35] -
                 d[19] * d[30] * d[35] + d[18] * d[31] * d[35];
    coeffs[19] = -d[29] * d[31] * d[33] + d[28] * d[32] * d[33] + d[29] * d[30] * d[34] - d[27] * d[32] * d[34] -
                 d[28] * d[30] * d[35] + d[27] * d[31] * d[35];
    coeffs[20] = -d[5] * d[7] * d[36] + d[4] * d[8] * d[36] + d[5] * d[6] * d[37] - d[3] * d[8] * d[37] -
                 d[4] * d[6] * d[38] + d[3] * d[7] * d[38] + d[2] * d[7] * d[39] - d[1] * d[8] * d[39] -
                 d[2] * d[6] * d[40] + d[0] * d[8] * d[40] + d[1] * d[6] * d[41] - d[0] * d[7] * d[41] -
                 d[2] * d[4] * d[42] + d[1] * d[5] * d[42] + d[2] * d[3] * d[43] - d[0] * d[5] * d[43] -
                 d[1] * d[3] * d[44] + d[0] * d[4] * d[44];
    coeffs[21] = d[8] * d[13] * d[36] - d[7] * d[14] * d[36] - d[5] * d[16] * d[36] + d[4] * d[17] * d[36] -
                 d[8] * d[12] * d[37] + d[6] * d[14] * d[37] + d[5] * d[15] * d[37] - d[3] * d[17] * d[37] +
                 d[7] * d[12] * d[38] - d[6] * d[13] * d[38] - d[4] * d[15] * d[38] + d[3] * d[16] * d[38] -
                 d[8] * d[10] * d[39] + d[7] * d[11] * d[39] + d[2] * d[16] * d[39] - d[1] * d[17] * d[39] +
                 d[8] * d[9] * d[40] - d[6] * d[11] * d[40] - d[2] * d[15] * d[40] + d[0] * d[17] * d[40] -
                 d[7] * d[9] * d[41] + d[6] * d[10] * d[41] + d[1] * d[15] * d[41] - d[0] * d[16] * d[41] +
                 d[5] * d[10] * d[42] - d[4] * d[11] * d[42] - d[2] * d[13] * d[42] + d[1] * d[14] * d[42] -
                 d[5] * d[9] * d[43] + d[3] * d[11] * d[43] + d[2] * d[12] * d[43] - d[0] * d[14] * d[43] +
                 d[4] * d[9] * d[44] - d[3] * d[10] * d[44] - d[1] * d[12] * d[44] + d[0] * d[13] * d[44];
    coeffs[22] = -d[14] * d[16] * d[36] + d[13] * d[17] * d[36] + d[14] * d[15] * d[37] - d[12] * d[17] * d[37] -
                 d[13] * d[15] * d[38] + d[12] * d[16] * d[38] + d[11] * d[16] * d[39] - d[10] * d[17] * d[39] -
                 d[11] * d[15] * d[40] + d[9] * d[17] * d[40] + d[10] * d[15] * d[41] - d[9] * d[16] * d[41] -
                 d[11] * d[13] * d[42] + d[10] * d[14] * d[42] + d[11] * d[12] * d[43] - d[9] * d[14] * d[43] -
                 d[10] * d[12] * d[44] + d[9] * d[13] * d[44];
    coeffs[23] = d[8] * d[22] * d[36] - d[7] * d[23] * d[36] - d[5] * d[25] * d[36] + d[4] * d[26] * d[36] -
                 d[8] * d[21] * d[37] + d[6] * d[23] * d[37] + d[5] * d[24] * d[37] - d[3] * d[26] * d[37] +
                 d[7] * d[21] * d[38] - d[6] * d[22] * d[38] - d[4] * d[24] * d[38] + d[3] * d[25] * d[38] -
                 d[8] * d[19] * d[39] + d[7] * d[20] * d[39] + d[2] * d[25] * d[39] - d[1] * d[26] * d[39] +
                 d[8] * d[18] * d[40] - d[6] * d[20] * d[40] - d[2] * d[24] * d[40] + d[0] * d[26] * d[40] -
                 d[7] * d[18] * d[41] + d[6] * d[19] * d[41] + d[1] * d[24] * d[41] - d[0] * d[25] * d[41] +
                 d[5] * d[19] * d[42] - d[4] * d[20] * d[42] - d[2] * d[22] * d[42] + d[1] * d[23] * d[42] -
                 d[5] * d[18] * d[43] + d[3] * d[20] * d[43] + d[2] * d[21] * d[43] - d[0] * d[23] * d[43] +
                 d[4] * d[18] * d[44] - d[3] * d[19] * d[44] - d[1] * d[21] * d[44] + d[0] * d[22] * d[44];
    coeffs[24] = d[17] * d[22] * d[36] - d[16] * d[23] * d[36] - d[14] * d[25] * d[36] + d[13] * d[26] * d[36] -
                 d[17] * d[21] * d[37] + d[15] * d[23] * d[37] + d[14] * d[24] * d[37] - d[12] * d[26] * d[37] +
                 d[16] * d[21] * d[38] - d[15] * d[22] * d[38] - d[13] * d[24] * d[38] + d[12] * d[25] * d[38] -
                 d[17] * d[19] * d[39] + d[16] * d[20] * d[39] + d[11] * d[25] * d[39] - d[10] * d[26] * d[39] +
                 d[17] * d[18] * d[40] - d[15] * d[20] * d[40] - d[11] * d[24] * d[40] + d[9] * d[26] * d[40] -
                 d[16] * d[18] * d[41] + d[15] * d[19] * d[41] + d[10] * d[24] * d[41] - d[9] * d[25] * d[41] +
                 d[14] * d[19] * d[42] - d[13] * d[20] * d[42] - d[11] * d[22] * d[42] + d[10] * d[23] * d[42] -
                 d[14] * d[18] * d[43] + d[12] * d[20] * d[43] + d[11] * d[21] * d[43] - d[9] * d[23] * d[43] +
                 d[13] * d[18] * d[44] - d[12] * d[19] * d[44] - d[10] * d[21] * d[44] + d[9] * d[22] * d[44];
    coeffs[25] = -d[23] * d[25] * d[36] + d[22] * d[26] * d[36] + d[23] * d[24] * d[37] - d[21] * d[26] * d[37] -
                 d[22] * d[24] * d[38] + d[21] * d[25] * d[38] + d[20] * d[25] * d[39] - d[19] * d[26] * d[39] -
                 d[20] * d[24] * d[40] + d[18] * d[26] * d[40] + d[19] * d[24] * d[41] - d[18] * d[25] * d[41] -
                 d[20] * d[22] * d[42] + d[19] * d[23] * d[42] + d[20] * d[21] * d[43] - d[18] * d[23] * d[43] -
                 d[19] * d[21] * d[44] + d[18] * d[22] * d[44];
    coeffs[26] = d[8] * d[31] * d[36] - d[7] * d[32] * d[36] - d[5] * d[34] * d[36] + d[4] * d[35] * d[36] -
                 d[8] * d[30] * d[37] + d[6] * d[32] * d[37] + d[5] * d[33] * d[37] - d[3] * d[35] * d[37] +
                 d[7] * d[30] * d[38] - d[6] * d[31] * d[38] - d[4] * d[33] * d[38] + d[3] * d[34] * d[38] -
                 d[8] * d[28] * d[39] + d[7] * d[29] * d[39] + d[2] * d[34] * d[39] - d[1] * d[35] * d[39] +
                 d[8] * d[27] * d[40] - d[6] * d[29] * d[40] - d[2] * d[33] * d[40] + d[0] * d[35] * d[40] -
                 d[7] * d[27] * d[41] + d[6] * d[28] * d[41] + d[1] * d[33] * d[41] - d[0] * d[34] * d[41] +
                 d[5] * d[28] * d[42] - d[4] * d[29] * d[42] - d[2] * d[31] * d[42] + d[1] * d[32] * d[42] -
                 d[5] * d[27] * d[43] + d[3] * d[29] * d[43] + d[2] * d[30] * d[43] - d[0] * d[32] * d[43] +
                 d[4] * d[27] * d[44] - d[3] * d[28] * d[44] - d[1] * d[30] * d[44] + d[0] * d[31] * d[44];
    coeffs[27] = d[17] * d[31] * d[36] - d[16] * d[32] * d[36] - d[14] * d[34] * d[36] + d[13] * d[35] * d[36] -
                 d[17] * d[30] * d[37] + d[15] * d[32] * d[37] + d[14] * d[33] * d[37] - d[12] * d[35] * d[37] +
                 d[16] * d[30] * d[38] - d[15] * d[31] * d[38] - d[13] * d[33] * d[38] + d[12] * d[34] * d[38] -
                 d[17] * d[28] * d[39] + d[16] * d[29] * d[39] + d[11] * d[34] * d[39] - d[10] * d[35] * d[39] +
                 d[17] * d[27] * d[40] - d[15] * d[29] * d[40] - d[11] * d[33] * d[40] + d[9] * d[35] * d[40] -
                 d[16] * d[27] * d[41] + d[15] * d[28] * d[41] + d[10] * d[33] * d[41] - d[9] * d[34] * d[41] +
                 d[14] * d[28] * d[42] - d[13] * d[29] * d[42] - d[11] * d[31] * d[42] + d[10] * d[32] * d[42] -
                 d[14] * d[27] * d[43] + d[12] * d[29] * d[43] + d[11] * d[30] * d[43] - d[9] * d[32] * d[43] +
                 d[13] * d[27] * d[44] - d[12] * d[28] * d[44] - d[10] * d[30] * d[44] + d[9] * d[31] * d[44];
    coeffs[28] = d[26] * d[31] * d[36] - d[25] * d[32] * d[36] - d[23] * d[34] * d[36] + d[22] * d[35] * d[36] -
                 d[26] * d[30] * d[37] + d[24] * d[32] * d[37] + d[23] * d[33] * d[37] - d[21] * d[35] * d[37] +
                 d[25] * d[30] * d[38] - d[24] * d[31] * d[38] - d[22] * d[33] * d[38] + d[21] * d[34] * d[38] -
                 d[26] * d[28] * d[39] + d[25] * d[29] * d[39] + d[20] * d[34] * d[39] - d[19] * d[35] * d[39] +
                 d[26] * d[27] * d[40] - d[24] * d[29] * d[40] - d[20] * d[33] * d[40] + d[18] * d[35] * d[40] -
                 d[25] * d[27] * d[41] + d[24] * d[28] * d[41] + d[19] * d[33] * d[41] - d[18] * d[34] * d[41] +
                 d[23] * d[28] * d[42] - d[22] * d[29] * d[42] - d[20] * d[31] * d[42] + d[19] * d[32] * d[42] -
                 d[23] * d[27] * d[43] + d[21] * d[29] * d[43] + d[20] * d[30] * d[43] - d[18] * d[32] * d[43] +
                 d[22] * d[27] * d[44] - d[21] * d[28] * d[44] - d[19] * d[30] * d[44] + d[18] * d[31] * d[44];
    coeffs[29] = -d[32] * d[34] * d[36] + d[31] * d[35] * d[36] + d[32] * d[33] * d[37] - d[30] * d[35] * d[37] -
                 d[31] * d[33] * d[38] + d[30] * d[34] * d[38] + d[29] * d[34] * d[39] - d[28] * d[35] * d[39] -
                 d[29] * d[33] * d[40] + d[27] * d[35] * d[40] + d[28] * d[33] * d[41] - d[27] * d[34] * d[41] -
                 d[29] * d[31] * d[42] + d[28] * d[32] * d[42] + d[29] * d[30] * d[43] - d[27] * d[32] * d[43] -
                 d[28] * d[30] * d[44] + d[27] * d[31] * d[44];
    coeffs[30] = -d[8] * d[37] * d[39] + d[7] * d[38] * d[39] + d[8] * d[36] * d[40] - d[6] * d[38] * d[40] -
                 d[7] * d[36] * d[41] + d[6] * d[37] * d[41] + d[5] * d[37] * d[42] - d[4] * d[38] * d[42] -
                 d[2] * d[40] * d[42] + d[1] * d[41] * d[42] - d[5] * d[36] * d[43] + d[3] * d[38] * d[43] +
                 d[2] * d[39] * d[43] - d[0] * d[41] * d[43] + d[4] * d[36] * d[44] - d[3] * d[37] * d[44] -
                 d[1] * d[39] * d[44] + d[0] * d[40] * d[44];
    coeffs[31] = -d[17] * d[37] * d[39] + d[16] * d[38] * d[39] + d[17] * d[36] * d[40] - d[15] * d[38] * d[40] -
                 d[16] * d[36] * d[41] + d[15] * d[37] * d[41] + d[14] * d[37] * d[42] - d[13] * d[38] * d[42] -
                 d[11] * d[40] * d[42] + d[10] * d[41] * d[42] - d[14] * d[36] * d[43] + d[12] * d[38] * d[43] +
                 d[11] * d[39] * d[43] - d[9] * d[41] * d[43] + d[13] * d[36] * d[44] - d[12] * d[37] * d[44] -
                 d[10] * d[39] * d[44] + d[9] * d[40] * d[44];
    coeffs[32] = -d[26] * d[37] * d[39] + d[25] * d[38] * d[39] + d[26] * d[36] * d[40] - d[24] * d[38] * d[40] -
                 d[25] * d[36] * d[41] + d[24] * d[37] * d[41] + d[23] * d[37] * d[42] - d[22] * d[38] * d[42] -
                 d[20] * d[40] * d[42] + d[19] * d[41] * d[42] - d[23] * d[36] * d[43] + d[21] * d[38] * d[43] +
                 d[20] * d[39] * d[43] - d[18] * d[41] * d[43] + d[22] * d[36] * d[44] - d[21] * d[37] * d[44] -
                 d[19] * d[39] * d[44] + d[18] * d[40] * d[44];
    coeffs[33] = -d[35] * d[37] * d[39] + d[34] * d[38] * d[39] + d[35] * d[36] * d[40] - d[33] * d[38] * d[40] -
                 d[34] * d[36] * d[41] + d[33] * d[37] * d[41] + d[32] * d[37] * d[42] - d[31] * d[38] * d[42] -
                 d[29] * d[40] * d[42] + d[28] * d[41] * d[42] - d[32] * d[36] * d[43] + d[30] * d[38] * d[43] +
                 d[29] * d[39] * d[43] - d[27] * d[41] * d[43] + d[31] * d[36] * d[44] - d[30] * d[37] * d[44] -
                 d[28] * d[39] * d[44] + d[27] * d[40] * d[44];
    coeffs[34] = -d[38] * d[40] * d[42] + d[37] * d[41] * d[42] + d[38] * d[39] * d[43] - d[36] * d[41] * d[43] -
                 d[37] * d[39] * d[44] + d[36] * d[40] * d[44];
    coeffs[35] = -2 * std::pow(d[2], 2) * d[4] + 2 * d[1] * d[2] * d[5] + 2 * d[2] * d[3] * d[5] -
                 2 * d[0] * std::pow(d[5], 2) - 4 * d[2] * d[4] * d[6] + 2 * d[1] * d[5] * d[6] +
                 2 * d[3] * d[5] * d[6] - 2 * d[4] * std::pow(d[6], 2) + 2 * d[1] * d[2] * d[7] +
                 2 * d[2] * d[3] * d[7] - 4 * d[0] * d[5] * d[7] + 2 * d[1] * d[6] * d[7] + 2 * d[3] * d[6] * d[7] -
                 2 * d[0] * std::pow(d[7], 2) - 2 * std::pow(d[1], 2) * d[8] - 4 * d[1] * d[3] * d[8] -
                 2 * std::pow(d[3], 2) * d[8] + 8 * d[0] * d[4] * d[8];
    coeffs[36] = -2 * std::pow(d[5], 2) * d[9] - 4 * d[5] * d[7] * d[9] - 2 * std::pow(d[7], 2) * d[9] +
                 8 * d[4] * d[8] * d[9] + 2 * d[2] * d[5] * d[10] + 2 * d[5] * d[6] * d[10] + 2 * d[2] * d[7] * d[10] +
                 2 * d[6] * d[7] * d[10] - 4 * d[1] * d[8] * d[10] - 4 * d[3] * d[8] * d[10] - 4 * d[2] * d[4] * d[11] +
                 2 * d[1] * d[5] * d[11] + 2 * d[3] * d[5] * d[11] - 4 * d[4] * d[6] * d[11] + 2 * d[1] * d[7] * d[11] +
                 2 * d[3] * d[7] * d[11] + 2 * d[2] * d[5] * d[12] + 2 * d[5] * d[6] * d[12] + 2 * d[2] * d[7] * d[12] +
                 2 * d[6] * d[7] * d[12] - 4 * d[1] * d[8] * d[12] - 4 * d[3] * d[8] * d[12] -
                 2 * std::pow(d[2], 2) * d[13] - 4 * d[2] * d[6] * d[13] - 2 * std::pow(d[6], 2) * d[13] +
                 8 * d[0] * d[8] * d[13] + 2 * d[1] * d[2] * d[14] + 2 * d[2] * d[3] * d[14] - 4 * d[0] * d[5] * d[14] +
                 2 * d[1] * d[6] * d[14] + 2 * d[3] * d[6] * d[14] - 4 * d[0] * d[7] * d[14] - 4 * d[2] * d[4] * d[15] +
                 2 * d[1] * d[5] * d[15] + 2 * d[3] * d[5] * d[15] - 4 * d[4] * d[6] * d[15] + 2 * d[1] * d[7] * d[15] +
                 2 * d[3] * d[7] * d[15] + 2 * d[1] * d[2] * d[16] + 2 * d[2] * d[3] * d[16] - 4 * d[0] * d[5] * d[16] +
                 2 * d[1] * d[6] * d[16] + 2 * d[3] * d[6] * d[16] - 4 * d[0] * d[7] * d[16] -
                 2 * std::pow(d[1], 2) * d[17] - 4 * d[1] * d[3] * d[17] - 2 * std::pow(d[3], 2) * d[17] +
                 8 * d[0] * d[4] * d[17];
    coeffs[37] =
        -2 * d[8] * std::pow(d[10], 2) + 2 * d[5] * d[10] * d[11] + 2 * d[7] * d[10] * d[11] -
        2 * d[4] * std::pow(d[11], 2) - 4 * d[8] * d[10] * d[12] + 2 * d[5] * d[11] * d[12] + 2 * d[7] * d[11] * d[12] -
        2 * d[8] * std::pow(d[12], 2) + 8 * d[8] * d[9] * d[13] - 4 * d[2] * d[11] * d[13] - 4 * d[6] * d[11] * d[13] -
        4 * d[5] * d[9] * d[14] - 4 * d[7] * d[9] * d[14] + 2 * d[2] * d[10] * d[14] + 2 * d[6] * d[10] * d[14] +
        2 * d[1] * d[11] * d[14] + 2 * d[3] * d[11] * d[14] + 2 * d[2] * d[12] * d[14] + 2 * d[6] * d[12] * d[14] -
        2 * d[0] * std::pow(d[14], 2) + 2 * d[5] * d[10] * d[15] + 2 * d[7] * d[10] * d[15] - 4 * d[4] * d[11] * d[15] +
        2 * d[5] * d[12] * d[15] + 2 * d[7] * d[12] * d[15] - 4 * d[2] * d[13] * d[15] - 4 * d[6] * d[13] * d[15] +
        2 * d[1] * d[14] * d[15] + 2 * d[3] * d[14] * d[15] - 2 * d[4] * std::pow(d[15], 2) - 4 * d[5] * d[9] * d[16] -
        4 * d[7] * d[9] * d[16] + 2 * d[2] * d[10] * d[16] + 2 * d[6] * d[10] * d[16] + 2 * d[1] * d[11] * d[16] +
        2 * d[3] * d[11] * d[16] + 2 * d[2] * d[12] * d[16] + 2 * d[6] * d[12] * d[16] - 4 * d[0] * d[14] * d[16] +
        2 * d[1] * d[15] * d[16] + 2 * d[3] * d[15] * d[16] - 2 * d[0] * std::pow(d[16], 2) + 8 * d[4] * d[9] * d[17] -
        4 * d[1] * d[10] * d[17] - 4 * d[3] * d[10] * d[17] - 4 * d[1] * d[12] * d[17] - 4 * d[3] * d[12] * d[17] +
        8 * d[0] * d[13] * d[17];
    coeffs[38] = -2 * std::pow(d[11], 2) * d[13] + 2 * d[10] * d[11] * d[14] + 2 * d[11] * d[12] * d[14] -
                 2 * d[9] * std::pow(d[14], 2) - 4 * d[11] * d[13] * d[15] + 2 * d[10] * d[14] * d[15] +
                 2 * d[12] * d[14] * d[15] - 2 * d[13] * std::pow(d[15], 2) + 2 * d[10] * d[11] * d[16] +
                 2 * d[11] * d[12] * d[16] - 4 * d[9] * d[14] * d[16] + 2 * d[10] * d[15] * d[16] +
                 2 * d[12] * d[15] * d[16] - 2 * d[9] * std::pow(d[16], 2) - 2 * std::pow(d[10], 2) * d[17] -
                 4 * d[10] * d[12] * d[17] - 2 * std::pow(d[12], 2) * d[17] + 8 * d[9] * d[13] * d[17];
    coeffs[39] = -2 * std::pow(d[5], 2) * d[18] - 4 * d[5] * d[7] * d[18] - 2 * std::pow(d[7], 2) * d[18] +
                 8 * d[4] * d[8] * d[18] + 2 * d[2] * d[5] * d[19] + 2 * d[5] * d[6] * d[19] + 2 * d[2] * d[7] * d[19] +
                 2 * d[6] * d[7] * d[19] - 4 * d[1] * d[8] * d[19] - 4 * d[3] * d[8] * d[19] - 4 * d[2] * d[4] * d[20] +
                 2 * d[1] * d[5] * d[20] + 2 * d[3] * d[5] * d[20] - 4 * d[4] * d[6] * d[20] + 2 * d[1] * d[7] * d[20] +
                 2 * d[3] * d[7] * d[20] + 2 * d[2] * d[5] * d[21] + 2 * d[5] * d[6] * d[21] + 2 * d[2] * d[7] * d[21] +
                 2 * d[6] * d[7] * d[21] - 4 * d[1] * d[8] * d[21] - 4 * d[3] * d[8] * d[21] -
                 2 * std::pow(d[2], 2) * d[22] - 4 * d[2] * d[6] * d[22] - 2 * std::pow(d[6], 2) * d[22] +
                 8 * d[0] * d[8] * d[22] + 2 * d[1] * d[2] * d[23] + 2 * d[2] * d[3] * d[23] - 4 * d[0] * d[5] * d[23] +
                 2 * d[1] * d[6] * d[23] + 2 * d[3] * d[6] * d[23] - 4 * d[0] * d[7] * d[23] - 4 * d[2] * d[4] * d[24] +
                 2 * d[1] * d[5] * d[24] + 2 * d[3] * d[5] * d[24] - 4 * d[4] * d[6] * d[24] + 2 * d[1] * d[7] * d[24] +
                 2 * d[3] * d[7] * d[24] + 2 * d[1] * d[2] * d[25] + 2 * d[2] * d[3] * d[25] - 4 * d[0] * d[5] * d[25] +
                 2 * d[1] * d[6] * d[25] + 2 * d[3] * d[6] * d[25] - 4 * d[0] * d[7] * d[25] -
                 2 * std::pow(d[1], 2) * d[26] - 4 * d[1] * d[3] * d[26] - 2 * std::pow(d[3], 2) * d[26] +
                 8 * d[0] * d[4] * d[26];
    coeffs[40] =
        8 * d[8] * d[13] * d[18] - 4 * d[5] * d[14] * d[18] - 4 * d[7] * d[14] * d[18] - 4 * d[5] * d[16] * d[18] -
        4 * d[7] * d[16] * d[18] + 8 * d[4] * d[17] * d[18] - 4 * d[8] * d[10] * d[19] + 2 * d[5] * d[11] * d[19] +
        2 * d[7] * d[11] * d[19] - 4 * d[8] * d[12] * d[19] + 2 * d[2] * d[14] * d[19] + 2 * d[6] * d[14] * d[19] +
        2 * d[5] * d[15] * d[19] + 2 * d[7] * d[15] * d[19] + 2 * d[2] * d[16] * d[19] + 2 * d[6] * d[16] * d[19] -
        4 * d[1] * d[17] * d[19] - 4 * d[3] * d[17] * d[19] + 2 * d[5] * d[10] * d[20] + 2 * d[7] * d[10] * d[20] -
        4 * d[4] * d[11] * d[20] + 2 * d[5] * d[12] * d[20] + 2 * d[7] * d[12] * d[20] - 4 * d[2] * d[13] * d[20] -
        4 * d[6] * d[13] * d[20] + 2 * d[1] * d[14] * d[20] + 2 * d[3] * d[14] * d[20] - 4 * d[4] * d[15] * d[20] +
        2 * d[1] * d[16] * d[20] + 2 * d[3] * d[16] * d[20] - 4 * d[8] * d[10] * d[21] + 2 * d[5] * d[11] * d[21] +
        2 * d[7] * d[11] * d[21] - 4 * d[8] * d[12] * d[21] + 2 * d[2] * d[14] * d[21] + 2 * d[6] * d[14] * d[21] +
        2 * d[5] * d[15] * d[21] + 2 * d[7] * d[15] * d[21] + 2 * d[2] * d[16] * d[21] + 2 * d[6] * d[16] * d[21] -
        4 * d[1] * d[17] * d[21] - 4 * d[3] * d[17] * d[21] + 8 * d[8] * d[9] * d[22] - 4 * d[2] * d[11] * d[22] -
        4 * d[6] * d[11] * d[22] - 4 * d[2] * d[15] * d[22] - 4 * d[6] * d[15] * d[22] + 8 * d[0] * d[17] * d[22] -
        4 * d[5] * d[9] * d[23] - 4 * d[7] * d[9] * d[23] + 2 * d[2] * d[10] * d[23] + 2 * d[6] * d[10] * d[23] +
        2 * d[1] * d[11] * d[23] + 2 * d[3] * d[11] * d[23] + 2 * d[2] * d[12] * d[23] + 2 * d[6] * d[12] * d[23] -
        4 * d[0] * d[14] * d[23] + 2 * d[1] * d[15] * d[23] + 2 * d[3] * d[15] * d[23] - 4 * d[0] * d[16] * d[23] +
        2 * d[5] * d[10] * d[24] + 2 * d[7] * d[10] * d[24] - 4 * d[4] * d[11] * d[24] + 2 * d[5] * d[12] * d[24] +
        2 * d[7] * d[12] * d[24] - 4 * d[2] * d[13] * d[24] - 4 * d[6] * d[13] * d[24] + 2 * d[1] * d[14] * d[24] +
        2 * d[3] * d[14] * d[24] - 4 * d[4] * d[15] * d[24] + 2 * d[1] * d[16] * d[24] + 2 * d[3] * d[16] * d[24] -
        4 * d[5] * d[9] * d[25] - 4 * d[7] * d[9] * d[25] + 2 * d[2] * d[10] * d[25] + 2 * d[6] * d[10] * d[25] +
        2 * d[1] * d[11] * d[25] + 2 * d[3] * d[11] * d[25] + 2 * d[2] * d[12] * d[25] + 2 * d[6] * d[12] * d[25] -
        4 * d[0] * d[14] * d[25] + 2 * d[1] * d[15] * d[25] + 2 * d[3] * d[15] * d[25] - 4 * d[0] * d[16] * d[25] +
        8 * d[4] * d[9] * d[26] - 4 * d[1] * d[10] * d[26] - 4 * d[3] * d[10] * d[26] - 4 * d[1] * d[12] * d[26] -
        4 * d[3] * d[12] * d[26] + 8 * d[0] * d[13] * d[26];
    coeffs[41] =
        -2 * std::pow(d[14], 2) * d[18] - 4 * d[14] * d[16] * d[18] - 2 * std::pow(d[16], 2) * d[18] +
        8 * d[13] * d[17] * d[18] + 2 * d[11] * d[14] * d[19] + 2 * d[14] * d[15] * d[19] + 2 * d[11] * d[16] * d[19] +
        2 * d[15] * d[16] * d[19] - 4 * d[10] * d[17] * d[19] - 4 * d[12] * d[17] * d[19] - 4 * d[11] * d[13] * d[20] +
        2 * d[10] * d[14] * d[20] + 2 * d[12] * d[14] * d[20] - 4 * d[13] * d[15] * d[20] + 2 * d[10] * d[16] * d[20] +
        2 * d[12] * d[16] * d[20] + 2 * d[11] * d[14] * d[21] + 2 * d[14] * d[15] * d[21] + 2 * d[11] * d[16] * d[21] +
        2 * d[15] * d[16] * d[21] - 4 * d[10] * d[17] * d[21] - 4 * d[12] * d[17] * d[21] -
        2 * std::pow(d[11], 2) * d[22] - 4 * d[11] * d[15] * d[22] - 2 * std::pow(d[15], 2) * d[22] +
        8 * d[9] * d[17] * d[22] + 2 * d[10] * d[11] * d[23] + 2 * d[11] * d[12] * d[23] - 4 * d[9] * d[14] * d[23] +
        2 * d[10] * d[15] * d[23] + 2 * d[12] * d[15] * d[23] - 4 * d[9] * d[16] * d[23] - 4 * d[11] * d[13] * d[24] +
        2 * d[10] * d[14] * d[24] + 2 * d[12] * d[14] * d[24] - 4 * d[13] * d[15] * d[24] + 2 * d[10] * d[16] * d[24] +
        2 * d[12] * d[16] * d[24] + 2 * d[10] * d[11] * d[25] + 2 * d[11] * d[12] * d[25] - 4 * d[9] * d[14] * d[25] +
        2 * d[10] * d[15] * d[25] + 2 * d[12] * d[15] * d[25] - 4 * d[9] * d[16] * d[25] -
        2 * std::pow(d[10], 2) * d[26] - 4 * d[10] * d[12] * d[26] - 2 * std::pow(d[12], 2) * d[26] +
        8 * d[9] * d[13] * d[26];
    coeffs[42] =
        -2 * d[8] * std::pow(d[19], 2) + 2 * d[5] * d[19] * d[20] + 2 * d[7] * d[19] * d[20] -
        2 * d[4] * std::pow(d[20], 2) - 4 * d[8] * d[19] * d[21] + 2 * d[5] * d[20] * d[21] + 2 * d[7] * d[20] * d[21] -
        2 * d[8] * std::pow(d[21], 2) + 8 * d[8] * d[18] * d[22] - 4 * d[2] * d[20] * d[22] - 4 * d[6] * d[20] * d[22] -
        4 * d[5] * d[18] * d[23] - 4 * d[7] * d[18] * d[23] + 2 * d[2] * d[19] * d[23] + 2 * d[6] * d[19] * d[23] +
        2 * d[1] * d[20] * d[23] + 2 * d[3] * d[20] * d[23] + 2 * d[2] * d[21] * d[23] + 2 * d[6] * d[21] * d[23] -
        2 * d[0] * std::pow(d[23], 2) + 2 * d[5] * d[19] * d[24] + 2 * d[7] * d[19] * d[24] - 4 * d[4] * d[20] * d[24] +
        2 * d[5] * d[21] * d[24] + 2 * d[7] * d[21] * d[24] - 4 * d[2] * d[22] * d[24] - 4 * d[6] * d[22] * d[24] +
        2 * d[1] * d[23] * d[24] + 2 * d[3] * d[23] * d[24] - 2 * d[4] * std::pow(d[24], 2) - 4 * d[5] * d[18] * d[25] -
        4 * d[7] * d[18] * d[25] + 2 * d[2] * d[19] * d[25] + 2 * d[6] * d[19] * d[25] + 2 * d[1] * d[20] * d[25] +
        2 * d[3] * d[20] * d[25] + 2 * d[2] * d[21] * d[25] + 2 * d[6] * d[21] * d[25] - 4 * d[0] * d[23] * d[25] +
        2 * d[1] * d[24] * d[25] + 2 * d[3] * d[24] * d[25] - 2 * d[0] * std::pow(d[25], 2) + 8 * d[4] * d[18] * d[26] -
        4 * d[1] * d[19] * d[26] - 4 * d[3] * d[19] * d[26] - 4 * d[1] * d[21] * d[26] - 4 * d[3] * d[21] * d[26] +
        8 * d[0] * d[22] * d[26];
    coeffs[43] =
        -2 * d[17] * std::pow(d[19], 2) + 2 * d[14] * d[19] * d[20] + 2 * d[16] * d[19] * d[20] -
        2 * d[13] * std::pow(d[20], 2) - 4 * d[17] * d[19] * d[21] + 2 * d[14] * d[20] * d[21] +
        2 * d[16] * d[20] * d[21] - 2 * d[17] * std::pow(d[21], 2) + 8 * d[17] * d[18] * d[22] -
        4 * d[11] * d[20] * d[22] - 4 * d[15] * d[20] * d[22] - 4 * d[14] * d[18] * d[23] - 4 * d[16] * d[18] * d[23] +
        2 * d[11] * d[19] * d[23] + 2 * d[15] * d[19] * d[23] + 2 * d[10] * d[20] * d[23] + 2 * d[12] * d[20] * d[23] +
        2 * d[11] * d[21] * d[23] + 2 * d[15] * d[21] * d[23] - 2 * d[9] * std::pow(d[23], 2) +
        2 * d[14] * d[19] * d[24] + 2 * d[16] * d[19] * d[24] - 4 * d[13] * d[20] * d[24] + 2 * d[14] * d[21] * d[24] +
        2 * d[16] * d[21] * d[24] - 4 * d[11] * d[22] * d[24] - 4 * d[15] * d[22] * d[24] + 2 * d[10] * d[23] * d[24] +
        2 * d[12] * d[23] * d[24] - 2 * d[13] * std::pow(d[24], 2) - 4 * d[14] * d[18] * d[25] -
        4 * d[16] * d[18] * d[25] + 2 * d[11] * d[19] * d[25] + 2 * d[15] * d[19] * d[25] + 2 * d[10] * d[20] * d[25] +
        2 * d[12] * d[20] * d[25] + 2 * d[11] * d[21] * d[25] + 2 * d[15] * d[21] * d[25] - 4 * d[9] * d[23] * d[25] +
        2 * d[10] * d[24] * d[25] + 2 * d[12] * d[24] * d[25] - 2 * d[9] * std::pow(d[25], 2) +
        8 * d[13] * d[18] * d[26] - 4 * d[10] * d[19] * d[26] - 4 * d[12] * d[19] * d[26] - 4 * d[10] * d[21] * d[26] -
        4 * d[12] * d[21] * d[26] + 8 * d[9] * d[22] * d[26];
    coeffs[44] = -2 * std::pow(d[20], 2) * d[22] + 2 * d[19] * d[20] * d[23] + 2 * d[20] * d[21] * d[23] -
                 2 * d[18] * std::pow(d[23], 2) - 4 * d[20] * d[22] * d[24] + 2 * d[19] * d[23] * d[24] +
                 2 * d[21] * d[23] * d[24] - 2 * d[22] * std::pow(d[24], 2) + 2 * d[19] * d[20] * d[25] +
                 2 * d[20] * d[21] * d[25] - 4 * d[18] * d[23] * d[25] + 2 * d[19] * d[24] * d[25] +
                 2 * d[21] * d[24] * d[25] - 2 * d[18] * std::pow(d[25], 2) - 2 * std::pow(d[19], 2) * d[26] -
                 4 * d[19] * d[21] * d[26] - 2 * std::pow(d[21], 2) * d[26] + 8 * d[18] * d[22] * d[26];
    coeffs[45] = -2 * std::pow(d[5], 2) * d[27] - 4 * d[5] * d[7] * d[27] - 2 * std::pow(d[7], 2) * d[27] +
                 8 * d[4] * d[8] * d[27] + 2 * d[2] * d[5] * d[28] + 2 * d[5] * d[6] * d[28] + 2 * d[2] * d[7] * d[28] +
                 2 * d[6] * d[7] * d[28] - 4 * d[1] * d[8] * d[28] - 4 * d[3] * d[8] * d[28] - 4 * d[2] * d[4] * d[29] +
                 2 * d[1] * d[5] * d[29] + 2 * d[3] * d[5] * d[29] - 4 * d[4] * d[6] * d[29] + 2 * d[1] * d[7] * d[29] +
                 2 * d[3] * d[7] * d[29] + 2 * d[2] * d[5] * d[30] + 2 * d[5] * d[6] * d[30] + 2 * d[2] * d[7] * d[30] +
                 2 * d[6] * d[7] * d[30] - 4 * d[1] * d[8] * d[30] - 4 * d[3] * d[8] * d[30] -
                 2 * std::pow(d[2], 2) * d[31] - 4 * d[2] * d[6] * d[31] - 2 * std::pow(d[6], 2) * d[31] +
                 8 * d[0] * d[8] * d[31] + 2 * d[1] * d[2] * d[32] + 2 * d[2] * d[3] * d[32] - 4 * d[0] * d[5] * d[32] +
                 2 * d[1] * d[6] * d[32] + 2 * d[3] * d[6] * d[32] - 4 * d[0] * d[7] * d[32] - 4 * d[2] * d[4] * d[33] +
                 2 * d[1] * d[5] * d[33] + 2 * d[3] * d[5] * d[33] - 4 * d[4] * d[6] * d[33] + 2 * d[1] * d[7] * d[33] +
                 2 * d[3] * d[7] * d[33] + 2 * d[1] * d[2] * d[34] + 2 * d[2] * d[3] * d[34] - 4 * d[0] * d[5] * d[34] +
                 2 * d[1] * d[6] * d[34] + 2 * d[3] * d[6] * d[34] - 4 * d[0] * d[7] * d[34] -
                 2 * std::pow(d[1], 2) * d[35] - 4 * d[1] * d[3] * d[35] - 2 * std::pow(d[3], 2) * d[35] +
                 8 * d[0] * d[4] * d[35];
    coeffs[46] =
        8 * d[8] * d[13] * d[27] - 4 * d[5] * d[14] * d[27] - 4 * d[7] * d[14] * d[27] - 4 * d[5] * d[16] * d[27] -
        4 * d[7] * d[16] * d[27] + 8 * d[4] * d[17] * d[27] - 4 * d[8] * d[10] * d[28] + 2 * d[5] * d[11] * d[28] +
        2 * d[7] * d[11] * d[28] - 4 * d[8] * d[12] * d[28] + 2 * d[2] * d[14] * d[28] + 2 * d[6] * d[14] * d[28] +
        2 * d[5] * d[15] * d[28] + 2 * d[7] * d[15] * d[28] + 2 * d[2] * d[16] * d[28] + 2 * d[6] * d[16] * d[28] -
        4 * d[1] * d[17] * d[28] - 4 * d[3] * d[17] * d[28] + 2 * d[5] * d[10] * d[29] + 2 * d[7] * d[10] * d[29] -
        4 * d[4] * d[11] * d[29] + 2 * d[5] * d[12] * d[29] + 2 * d[7] * d[12] * d[29] - 4 * d[2] * d[13] * d[29] -
        4 * d[6] * d[13] * d[29] + 2 * d[1] * d[14] * d[29] + 2 * d[3] * d[14] * d[29] - 4 * d[4] * d[15] * d[29] +
        2 * d[1] * d[16] * d[29] + 2 * d[3] * d[16] * d[29] - 4 * d[8] * d[10] * d[30] + 2 * d[5] * d[11] * d[30] +
        2 * d[7] * d[11] * d[30] - 4 * d[8] * d[12] * d[30] + 2 * d[2] * d[14] * d[30] + 2 * d[6] * d[14] * d[30] +
        2 * d[5] * d[15] * d[30] + 2 * d[7] * d[15] * d[30] + 2 * d[2] * d[16] * d[30] + 2 * d[6] * d[16] * d[30] -
        4 * d[1] * d[17] * d[30] - 4 * d[3] * d[17] * d[30] + 8 * d[8] * d[9] * d[31] - 4 * d[2] * d[11] * d[31] -
        4 * d[6] * d[11] * d[31] - 4 * d[2] * d[15] * d[31] - 4 * d[6] * d[15] * d[31] + 8 * d[0] * d[17] * d[31] -
        4 * d[5] * d[9] * d[32] - 4 * d[7] * d[9] * d[32] + 2 * d[2] * d[10] * d[32] + 2 * d[6] * d[10] * d[32] +
        2 * d[1] * d[11] * d[32] + 2 * d[3] * d[11] * d[32] + 2 * d[2] * d[12] * d[32] + 2 * d[6] * d[12] * d[32] -
        4 * d[0] * d[14] * d[32] + 2 * d[1] * d[15] * d[32] + 2 * d[3] * d[15] * d[32] - 4 * d[0] * d[16] * d[32] +
        2 * d[5] * d[10] * d[33] + 2 * d[7] * d[10] * d[33] - 4 * d[4] * d[11] * d[33] + 2 * d[5] * d[12] * d[33] +
        2 * d[7] * d[12] * d[33] - 4 * d[2] * d[13] * d[33] - 4 * d[6] * d[13] * d[33] + 2 * d[1] * d[14] * d[33] +
        2 * d[3] * d[14] * d[33] - 4 * d[4] * d[15] * d[33] + 2 * d[1] * d[16] * d[33] + 2 * d[3] * d[16] * d[33] -
        4 * d[5] * d[9] * d[34] - 4 * d[7] * d[9] * d[34] + 2 * d[2] * d[10] * d[34] + 2 * d[6] * d[10] * d[34] +
        2 * d[1] * d[11] * d[34] + 2 * d[3] * d[11] * d[34] + 2 * d[2] * d[12] * d[34] + 2 * d[6] * d[12] * d[34] -
        4 * d[0] * d[14] * d[34] + 2 * d[1] * d[15] * d[34] + 2 * d[3] * d[15] * d[34] - 4 * d[0] * d[16] * d[34] +
        8 * d[4] * d[9] * d[35] - 4 * d[1] * d[10] * d[35] - 4 * d[3] * d[10] * d[35] - 4 * d[1] * d[12] * d[35] -
        4 * d[3] * d[12] * d[35] + 8 * d[0] * d[13] * d[35];
    coeffs[47] =
        -2 * std::pow(d[14], 2) * d[27] - 4 * d[14] * d[16] * d[27] - 2 * std::pow(d[16], 2) * d[27] +
        8 * d[13] * d[17] * d[27] + 2 * d[11] * d[14] * d[28] + 2 * d[14] * d[15] * d[28] + 2 * d[11] * d[16] * d[28] +
        2 * d[15] * d[16] * d[28] - 4 * d[10] * d[17] * d[28] - 4 * d[12] * d[17] * d[28] - 4 * d[11] * d[13] * d[29] +
        2 * d[10] * d[14] * d[29] + 2 * d[12] * d[14] * d[29] - 4 * d[13] * d[15] * d[29] + 2 * d[10] * d[16] * d[29] +
        2 * d[12] * d[16] * d[29] + 2 * d[11] * d[14] * d[30] + 2 * d[14] * d[15] * d[30] + 2 * d[11] * d[16] * d[30] +
        2 * d[15] * d[16] * d[30] - 4 * d[10] * d[17] * d[30] - 4 * d[12] * d[17] * d[30] -
        2 * std::pow(d[11], 2) * d[31] - 4 * d[11] * d[15] * d[31] - 2 * std::pow(d[15], 2) * d[31] +
        8 * d[9] * d[17] * d[31] + 2 * d[10] * d[11] * d[32] + 2 * d[11] * d[12] * d[32] - 4 * d[9] * d[14] * d[32] +
        2 * d[10] * d[15] * d[32] + 2 * d[12] * d[15] * d[32] - 4 * d[9] * d[16] * d[32] - 4 * d[11] * d[13] * d[33] +
        2 * d[10] * d[14] * d[33] + 2 * d[12] * d[14] * d[33] - 4 * d[13] * d[15] * d[33] + 2 * d[10] * d[16] * d[33] +
        2 * d[12] * d[16] * d[33] + 2 * d[10] * d[11] * d[34] + 2 * d[11] * d[12] * d[34] - 4 * d[9] * d[14] * d[34] +
        2 * d[10] * d[15] * d[34] + 2 * d[12] * d[15] * d[34] - 4 * d[9] * d[16] * d[34] -
        2 * std::pow(d[10], 2) * d[35] - 4 * d[10] * d[12] * d[35] - 2 * std::pow(d[12], 2) * d[35] +
        8 * d[9] * d[13] * d[35];
    coeffs[48] =
        8 * d[8] * d[22] * d[27] - 4 * d[5] * d[23] * d[27] - 4 * d[7] * d[23] * d[27] - 4 * d[5] * d[25] * d[27] -
        4 * d[7] * d[25] * d[27] + 8 * d[4] * d[26] * d[27] - 4 * d[8] * d[19] * d[28] + 2 * d[5] * d[20] * d[28] +
        2 * d[7] * d[20] * d[28] - 4 * d[8] * d[21] * d[28] + 2 * d[2] * d[23] * d[28] + 2 * d[6] * d[23] * d[28] +
        2 * d[5] * d[24] * d[28] + 2 * d[7] * d[24] * d[28] + 2 * d[2] * d[25] * d[28] + 2 * d[6] * d[25] * d[28] -
        4 * d[1] * d[26] * d[28] - 4 * d[3] * d[26] * d[28] + 2 * d[5] * d[19] * d[29] + 2 * d[7] * d[19] * d[29] -
        4 * d[4] * d[20] * d[29] + 2 * d[5] * d[21] * d[29] + 2 * d[7] * d[21] * d[29] - 4 * d[2] * d[22] * d[29] -
        4 * d[6] * d[22] * d[29] + 2 * d[1] * d[23] * d[29] + 2 * d[3] * d[23] * d[29] - 4 * d[4] * d[24] * d[29] +
        2 * d[1] * d[25] * d[29] + 2 * d[3] * d[25] * d[29] - 4 * d[8] * d[19] * d[30] + 2 * d[5] * d[20] * d[30] +
        2 * d[7] * d[20] * d[30] - 4 * d[8] * d[21] * d[30] + 2 * d[2] * d[23] * d[30] + 2 * d[6] * d[23] * d[30] +
        2 * d[5] * d[24] * d[30] + 2 * d[7] * d[24] * d[30] + 2 * d[2] * d[25] * d[30] + 2 * d[6] * d[25] * d[30] -
        4 * d[1] * d[26] * d[30] - 4 * d[3] * d[26] * d[30] + 8 * d[8] * d[18] * d[31] - 4 * d[2] * d[20] * d[31] -
        4 * d[6] * d[20] * d[31] - 4 * d[2] * d[24] * d[31] - 4 * d[6] * d[24] * d[31] + 8 * d[0] * d[26] * d[31] -
        4 * d[5] * d[18] * d[32] - 4 * d[7] * d[18] * d[32] + 2 * d[2] * d[19] * d[32] + 2 * d[6] * d[19] * d[32] +
        2 * d[1] * d[20] * d[32] + 2 * d[3] * d[20] * d[32] + 2 * d[2] * d[21] * d[32] + 2 * d[6] * d[21] * d[32] -
        4 * d[0] * d[23] * d[32] + 2 * d[1] * d[24] * d[32] + 2 * d[3] * d[24] * d[32] - 4 * d[0] * d[25] * d[32] +
        2 * d[5] * d[19] * d[33] + 2 * d[7] * d[19] * d[33] - 4 * d[4] * d[20] * d[33] + 2 * d[5] * d[21] * d[33] +
        2 * d[7] * d[21] * d[33] - 4 * d[2] * d[22] * d[33] - 4 * d[6] * d[22] * d[33] + 2 * d[1] * d[23] * d[33] +
        2 * d[3] * d[23] * d[33] - 4 * d[4] * d[24] * d[33] + 2 * d[1] * d[25] * d[33] + 2 * d[3] * d[25] * d[33] -
        4 * d[5] * d[18] * d[34] - 4 * d[7] * d[18] * d[34] + 2 * d[2] * d[19] * d[34] + 2 * d[6] * d[19] * d[34] +
        2 * d[1] * d[20] * d[34] + 2 * d[3] * d[20] * d[34] + 2 * d[2] * d[21] * d[34] + 2 * d[6] * d[21] * d[34] -
        4 * d[0] * d[23] * d[34] + 2 * d[1] * d[24] * d[34] + 2 * d[3] * d[24] * d[34] - 4 * d[0] * d[25] * d[34] +
        8 * d[4] * d[18] * d[35] - 4 * d[1] * d[19] * d[35] - 4 * d[3] * d[19] * d[35] - 4 * d[1] * d[21] * d[35] -
        4 * d[3] * d[21] * d[35] + 8 * d[0] * d[22] * d[35];
    coeffs[49] =
        8 * d[17] * d[22] * d[27] - 4 * d[14] * d[23] * d[27] - 4 * d[16] * d[23] * d[27] - 4 * d[14] * d[25] * d[27] -
        4 * d[16] * d[25] * d[27] + 8 * d[13] * d[26] * d[27] - 4 * d[17] * d[19] * d[28] + 2 * d[14] * d[20] * d[28] +
        2 * d[16] * d[20] * d[28] - 4 * d[17] * d[21] * d[28] + 2 * d[11] * d[23] * d[28] + 2 * d[15] * d[23] * d[28] +
        2 * d[14] * d[24] * d[28] + 2 * d[16] * d[24] * d[28] + 2 * d[11] * d[25] * d[28] + 2 * d[15] * d[25] * d[28] -
        4 * d[10] * d[26] * d[28] - 4 * d[12] * d[26] * d[28] + 2 * d[14] * d[19] * d[29] + 2 * d[16] * d[19] * d[29] -
        4 * d[13] * d[20] * d[29] + 2 * d[14] * d[21] * d[29] + 2 * d[16] * d[21] * d[29] - 4 * d[11] * d[22] * d[29] -
        4 * d[15] * d[22] * d[29] + 2 * d[10] * d[23] * d[29] + 2 * d[12] * d[23] * d[29] - 4 * d[13] * d[24] * d[29] +
        2 * d[10] * d[25] * d[29] + 2 * d[12] * d[25] * d[29] - 4 * d[17] * d[19] * d[30] + 2 * d[14] * d[20] * d[30] +
        2 * d[16] * d[20] * d[30] - 4 * d[17] * d[21] * d[30] + 2 * d[11] * d[23] * d[30] + 2 * d[15] * d[23] * d[30] +
        2 * d[14] * d[24] * d[30] + 2 * d[16] * d[24] * d[30] + 2 * d[11] * d[25] * d[30] + 2 * d[15] * d[25] * d[30] -
        4 * d[10] * d[26] * d[30] - 4 * d[12] * d[26] * d[30] + 8 * d[17] * d[18] * d[31] - 4 * d[11] * d[20] * d[31] -
        4 * d[15] * d[20] * d[31] - 4 * d[11] * d[24] * d[31] - 4 * d[15] * d[24] * d[31] + 8 * d[9] * d[26] * d[31] -
        4 * d[14] * d[18] * d[32] - 4 * d[16] * d[18] * d[32] + 2 * d[11] * d[19] * d[32] + 2 * d[15] * d[19] * d[32] +
        2 * d[10] * d[20] * d[32] + 2 * d[12] * d[20] * d[32] + 2 * d[11] * d[21] * d[32] + 2 * d[15] * d[21] * d[32] -
        4 * d[9] * d[23] * d[32] + 2 * d[10] * d[24] * d[32] + 2 * d[12] * d[24] * d[32] - 4 * d[9] * d[25] * d[32] +
        2 * d[14] * d[19] * d[33] + 2 * d[16] * d[19] * d[33] - 4 * d[13] * d[20] * d[33] + 2 * d[14] * d[21] * d[33] +
        2 * d[16] * d[21] * d[33] - 4 * d[11] * d[22] * d[33] - 4 * d[15] * d[22] * d[33] + 2 * d[10] * d[23] * d[33] +
        2 * d[12] * d[23] * d[33] - 4 * d[13] * d[24] * d[33] + 2 * d[10] * d[25] * d[33] + 2 * d[12] * d[25] * d[33] -
        4 * d[14] * d[18] * d[34] - 4 * d[16] * d[18] * d[34] + 2 * d[11] * d[19] * d[34] + 2 * d[15] * d[19] * d[34] +
        2 * d[10] * d[20] * d[34] + 2 * d[12] * d[20] * d[34] + 2 * d[11] * d[21] * d[34] + 2 * d[15] * d[21] * d[34] -
        4 * d[9] * d[23] * d[34] + 2 * d[10] * d[24] * d[34] + 2 * d[12] * d[24] * d[34] - 4 * d[9] * d[25] * d[34] +
        8 * d[13] * d[18] * d[35] - 4 * d[10] * d[19] * d[35] - 4 * d[12] * d[19] * d[35] - 4 * d[10] * d[21] * d[35] -
        4 * d[12] * d[21] * d[35] + 8 * d[9] * d[22] * d[35];
    coeffs[50] =
        -2 * std::pow(d[23], 2) * d[27] - 4 * d[23] * d[25] * d[27] - 2 * std::pow(d[25], 2) * d[27] +
        8 * d[22] * d[26] * d[27] + 2 * d[20] * d[23] * d[28] + 2 * d[23] * d[24] * d[28] + 2 * d[20] * d[25] * d[28] +
        2 * d[24] * d[25] * d[28] - 4 * d[19] * d[26] * d[28] - 4 * d[21] * d[26] * d[28] - 4 * d[20] * d[22] * d[29] +
        2 * d[19] * d[23] * d[29] + 2 * d[21] * d[23] * d[29] - 4 * d[22] * d[24] * d[29] + 2 * d[19] * d[25] * d[29] +
        2 * d[21] * d[25] * d[29] + 2 * d[20] * d[23] * d[30] + 2 * d[23] * d[24] * d[30] + 2 * d[20] * d[25] * d[30] +
        2 * d[24] * d[25] * d[30] - 4 * d[19] * d[26] * d[30] - 4 * d[21] * d[26] * d[30] -
        2 * std::pow(d[20], 2) * d[31] - 4 * d[20] * d[24] * d[31] - 2 * std::pow(d[24], 2) * d[31] +
        8 * d[18] * d[26] * d[31] + 2 * d[19] * d[20] * d[32] + 2 * d[20] * d[21] * d[32] - 4 * d[18] * d[23] * d[32] +
        2 * d[19] * d[24] * d[32] + 2 * d[21] * d[24] * d[32] - 4 * d[18] * d[25] * d[32] - 4 * d[20] * d[22] * d[33] +
        2 * d[19] * d[23] * d[33] + 2 * d[21] * d[23] * d[33] - 4 * d[22] * d[24] * d[33] + 2 * d[19] * d[25] * d[33] +
        2 * d[21] * d[25] * d[33] + 2 * d[19] * d[20] * d[34] + 2 * d[20] * d[21] * d[34] - 4 * d[18] * d[23] * d[34] +
        2 * d[19] * d[24] * d[34] + 2 * d[21] * d[24] * d[34] - 4 * d[18] * d[25] * d[34] -
        2 * std::pow(d[19], 2) * d[35] - 4 * d[19] * d[21] * d[35] - 2 * std::pow(d[21], 2) * d[35] +
        8 * d[18] * d[22] * d[35];
    coeffs[51] =
        -2 * d[8] * std::pow(d[28], 2) + 2 * d[5] * d[28] * d[29] + 2 * d[7] * d[28] * d[29] -
        2 * d[4] * std::pow(d[29], 2) - 4 * d[8] * d[28] * d[30] + 2 * d[5] * d[29] * d[30] + 2 * d[7] * d[29] * d[30] -
        2 * d[8] * std::pow(d[30], 2) + 8 * d[8] * d[27] * d[31] - 4 * d[2] * d[29] * d[31] - 4 * d[6] * d[29] * d[31] -
        4 * d[5] * d[27] * d[32] - 4 * d[7] * d[27] * d[32] + 2 * d[2] * d[28] * d[32] + 2 * d[6] * d[28] * d[32] +
        2 * d[1] * d[29] * d[32] + 2 * d[3] * d[29] * d[32] + 2 * d[2] * d[30] * d[32] + 2 * d[6] * d[30] * d[32] -
        2 * d[0] * std::pow(d[32], 2) + 2 * d[5] * d[28] * d[33] + 2 * d[7] * d[28] * d[33] - 4 * d[4] * d[29] * d[33] +
        2 * d[5] * d[30] * d[33] + 2 * d[7] * d[30] * d[33] - 4 * d[2] * d[31] * d[33] - 4 * d[6] * d[31] * d[33] +
        2 * d[1] * d[32] * d[33] + 2 * d[3] * d[32] * d[33] - 2 * d[4] * std::pow(d[33], 2) - 4 * d[5] * d[27] * d[34] -
        4 * d[7] * d[27] * d[34] + 2 * d[2] * d[28] * d[34] + 2 * d[6] * d[28] * d[34] + 2 * d[1] * d[29] * d[34] +
        2 * d[3] * d[29] * d[34] + 2 * d[2] * d[30] * d[34] + 2 * d[6] * d[30] * d[34] - 4 * d[0] * d[32] * d[34] +
        2 * d[1] * d[33] * d[34] + 2 * d[3] * d[33] * d[34] - 2 * d[0] * std::pow(d[34], 2) + 8 * d[4] * d[27] * d[35] -
        4 * d[1] * d[28] * d[35] - 4 * d[3] * d[28] * d[35] - 4 * d[1] * d[30] * d[35] - 4 * d[3] * d[30] * d[35] +
        8 * d[0] * d[31] * d[35];
    coeffs[52] =
        -2 * d[17] * std::pow(d[28], 2) + 2 * d[14] * d[28] * d[29] + 2 * d[16] * d[28] * d[29] -
        2 * d[13] * std::pow(d[29], 2) - 4 * d[17] * d[28] * d[30] + 2 * d[14] * d[29] * d[30] +
        2 * d[16] * d[29] * d[30] - 2 * d[17] * std::pow(d[30], 2) + 8 * d[17] * d[27] * d[31] -
        4 * d[11] * d[29] * d[31] - 4 * d[15] * d[29] * d[31] - 4 * d[14] * d[27] * d[32] - 4 * d[16] * d[27] * d[32] +
        2 * d[11] * d[28] * d[32] + 2 * d[15] * d[28] * d[32] + 2 * d[10] * d[29] * d[32] + 2 * d[12] * d[29] * d[32] +
        2 * d[11] * d[30] * d[32] + 2 * d[15] * d[30] * d[32] - 2 * d[9] * std::pow(d[32], 2) +
        2 * d[14] * d[28] * d[33] + 2 * d[16] * d[28] * d[33] - 4 * d[13] * d[29] * d[33] + 2 * d[14] * d[30] * d[33] +
        2 * d[16] * d[30] * d[33] - 4 * d[11] * d[31] * d[33] - 4 * d[15] * d[31] * d[33] + 2 * d[10] * d[32] * d[33] +
        2 * d[12] * d[32] * d[33] - 2 * d[13] * std::pow(d[33], 2) - 4 * d[14] * d[27] * d[34] -
        4 * d[16] * d[27] * d[34] + 2 * d[11] * d[28] * d[34] + 2 * d[15] * d[28] * d[34] + 2 * d[10] * d[29] * d[34] +
        2 * d[12] * d[29] * d[34] + 2 * d[11] * d[30] * d[34] + 2 * d[15] * d[30] * d[34] - 4 * d[9] * d[32] * d[34] +
        2 * d[10] * d[33] * d[34] + 2 * d[12] * d[33] * d[34] - 2 * d[9] * std::pow(d[34], 2) +
        8 * d[13] * d[27] * d[35] - 4 * d[10] * d[28] * d[35] - 4 * d[12] * d[28] * d[35] - 4 * d[10] * d[30] * d[35] -
        4 * d[12] * d[30] * d[35] + 8 * d[9] * d[31] * d[35];
    coeffs[53] =
        -2 * d[26] * std::pow(d[28], 2) + 2 * d[23] * d[28] * d[29] + 2 * d[25] * d[28] * d[29] -
        2 * d[22] * std::pow(d[29], 2) - 4 * d[26] * d[28] * d[30] + 2 * d[23] * d[29] * d[30] +
        2 * d[25] * d[29] * d[30] - 2 * d[26] * std::pow(d[30], 2) + 8 * d[26] * d[27] * d[31] -
        4 * d[20] * d[29] * d[31] - 4 * d[24] * d[29] * d[31] - 4 * d[23] * d[27] * d[32] - 4 * d[25] * d[27] * d[32] +
        2 * d[20] * d[28] * d[32] + 2 * d[24] * d[28] * d[32] + 2 * d[19] * d[29] * d[32] + 2 * d[21] * d[29] * d[32] +
        2 * d[20] * d[30] * d[32] + 2 * d[24] * d[30] * d[32] - 2 * d[18] * std::pow(d[32], 2) +
        2 * d[23] * d[28] * d[33] + 2 * d[25] * d[28] * d[33] - 4 * d[22] * d[29] * d[33] + 2 * d[23] * d[30] * d[33] +
        2 * d[25] * d[30] * d[33] - 4 * d[20] * d[31] * d[33] - 4 * d[24] * d[31] * d[33] + 2 * d[19] * d[32] * d[33] +
        2 * d[21] * d[32] * d[33] - 2 * d[22] * std::pow(d[33], 2) - 4 * d[23] * d[27] * d[34] -
        4 * d[25] * d[27] * d[34] + 2 * d[20] * d[28] * d[34] + 2 * d[24] * d[28] * d[34] + 2 * d[19] * d[29] * d[34] +
        2 * d[21] * d[29] * d[34] + 2 * d[20] * d[30] * d[34] + 2 * d[24] * d[30] * d[34] - 4 * d[18] * d[32] * d[34] +
        2 * d[19] * d[33] * d[34] + 2 * d[21] * d[33] * d[34] - 2 * d[18] * std::pow(d[34], 2) +
        8 * d[22] * d[27] * d[35] - 4 * d[19] * d[28] * d[35] - 4 * d[21] * d[28] * d[35] - 4 * d[19] * d[30] * d[35] -
        4 * d[21] * d[30] * d[35] + 8 * d[18] * d[31] * d[35];
    coeffs[54] = -2 * std::pow(d[29], 2) * d[31] + 2 * d[28] * d[29] * d[32] + 2 * d[29] * d[30] * d[32] -
                 2 * d[27] * std::pow(d[32], 2) - 4 * d[29] * d[31] * d[33] + 2 * d[28] * d[32] * d[33] +
                 2 * d[30] * d[32] * d[33] - 2 * d[31] * std::pow(d[33], 2) + 2 * d[28] * d[29] * d[34] +
                 2 * d[29] * d[30] * d[34] - 4 * d[27] * d[32] * d[34] + 2 * d[28] * d[33] * d[34] +
                 2 * d[30] * d[33] * d[34] - 2 * d[27] * std::pow(d[34], 2) - 2 * std::pow(d[28], 2) * d[35] -
                 4 * d[28] * d[30] * d[35] - 2 * std::pow(d[30], 2) * d[35] + 8 * d[27] * d[31] * d[35];
    coeffs[55] = -2 * std::pow(d[5], 2) * d[36] - 4 * d[5] * d[7] * d[36] - 2 * std::pow(d[7], 2) * d[36] +
                 8 * d[4] * d[8] * d[36] + 2 * d[2] * d[5] * d[37] + 2 * d[5] * d[6] * d[37] + 2 * d[2] * d[7] * d[37] +
                 2 * d[6] * d[7] * d[37] - 4 * d[1] * d[8] * d[37] - 4 * d[3] * d[8] * d[37] - 4 * d[2] * d[4] * d[38] +
                 2 * d[1] * d[5] * d[38] + 2 * d[3] * d[5] * d[38] - 4 * d[4] * d[6] * d[38] + 2 * d[1] * d[7] * d[38] +
                 2 * d[3] * d[7] * d[38] + 2 * d[2] * d[5] * d[39] + 2 * d[5] * d[6] * d[39] + 2 * d[2] * d[7] * d[39] +
                 2 * d[6] * d[7] * d[39] - 4 * d[1] * d[8] * d[39] - 4 * d[3] * d[8] * d[39] -
                 2 * std::pow(d[2], 2) * d[40] - 4 * d[2] * d[6] * d[40] - 2 * std::pow(d[6], 2) * d[40] +
                 8 * d[0] * d[8] * d[40] + 2 * d[1] * d[2] * d[41] + 2 * d[2] * d[3] * d[41] - 4 * d[0] * d[5] * d[41] +
                 2 * d[1] * d[6] * d[41] + 2 * d[3] * d[6] * d[41] - 4 * d[0] * d[7] * d[41] - 4 * d[2] * d[4] * d[42] +
                 2 * d[1] * d[5] * d[42] + 2 * d[3] * d[5] * d[42] - 4 * d[4] * d[6] * d[42] + 2 * d[1] * d[7] * d[42] +
                 2 * d[3] * d[7] * d[42] + 2 * d[1] * d[2] * d[43] + 2 * d[2] * d[3] * d[43] - 4 * d[0] * d[5] * d[43] +
                 2 * d[1] * d[6] * d[43] + 2 * d[3] * d[6] * d[43] - 4 * d[0] * d[7] * d[43] -
                 2 * std::pow(d[1], 2) * d[44] - 4 * d[1] * d[3] * d[44] - 2 * std::pow(d[3], 2) * d[44] +
                 8 * d[0] * d[4] * d[44];
    coeffs[56] =
        8 * d[8] * d[13] * d[36] - 4 * d[5] * d[14] * d[36] - 4 * d[7] * d[14] * d[36] - 4 * d[5] * d[16] * d[36] -
        4 * d[7] * d[16] * d[36] + 8 * d[4] * d[17] * d[36] - 4 * d[8] * d[10] * d[37] + 2 * d[5] * d[11] * d[37] +
        2 * d[7] * d[11] * d[37] - 4 * d[8] * d[12] * d[37] + 2 * d[2] * d[14] * d[37] + 2 * d[6] * d[14] * d[37] +
        2 * d[5] * d[15] * d[37] + 2 * d[7] * d[15] * d[37] + 2 * d[2] * d[16] * d[37] + 2 * d[6] * d[16] * d[37] -
        4 * d[1] * d[17] * d[37] - 4 * d[3] * d[17] * d[37] + 2 * d[5] * d[10] * d[38] + 2 * d[7] * d[10] * d[38] -
        4 * d[4] * d[11] * d[38] + 2 * d[5] * d[12] * d[38] + 2 * d[7] * d[12] * d[38] - 4 * d[2] * d[13] * d[38] -
        4 * d[6] * d[13] * d[38] + 2 * d[1] * d[14] * d[38] + 2 * d[3] * d[14] * d[38] - 4 * d[4] * d[15] * d[38] +
        2 * d[1] * d[16] * d[38] + 2 * d[3] * d[16] * d[38] - 4 * d[8] * d[10] * d[39] + 2 * d[5] * d[11] * d[39] +
        2 * d[7] * d[11] * d[39] - 4 * d[8] * d[12] * d[39] + 2 * d[2] * d[14] * d[39] + 2 * d[6] * d[14] * d[39] +
        2 * d[5] * d[15] * d[39] + 2 * d[7] * d[15] * d[39] + 2 * d[2] * d[16] * d[39] + 2 * d[6] * d[16] * d[39] -
        4 * d[1] * d[17] * d[39] - 4 * d[3] * d[17] * d[39] + 8 * d[8] * d[9] * d[40] - 4 * d[2] * d[11] * d[40] -
        4 * d[6] * d[11] * d[40] - 4 * d[2] * d[15] * d[40] - 4 * d[6] * d[15] * d[40] + 8 * d[0] * d[17] * d[40] -
        4 * d[5] * d[9] * d[41] - 4 * d[7] * d[9] * d[41] + 2 * d[2] * d[10] * d[41] + 2 * d[6] * d[10] * d[41] +
        2 * d[1] * d[11] * d[41] + 2 * d[3] * d[11] * d[41] + 2 * d[2] * d[12] * d[41] + 2 * d[6] * d[12] * d[41] -
        4 * d[0] * d[14] * d[41] + 2 * d[1] * d[15] * d[41] + 2 * d[3] * d[15] * d[41] - 4 * d[0] * d[16] * d[41] +
        2 * d[5] * d[10] * d[42] + 2 * d[7] * d[10] * d[42] - 4 * d[4] * d[11] * d[42] + 2 * d[5] * d[12] * d[42] +
        2 * d[7] * d[12] * d[42] - 4 * d[2] * d[13] * d[42] - 4 * d[6] * d[13] * d[42] + 2 * d[1] * d[14] * d[42] +
        2 * d[3] * d[14] * d[42] - 4 * d[4] * d[15] * d[42] + 2 * d[1] * d[16] * d[42] + 2 * d[3] * d[16] * d[42] -
        4 * d[5] * d[9] * d[43] - 4 * d[7] * d[9] * d[43] + 2 * d[2] * d[10] * d[43] + 2 * d[6] * d[10] * d[43] +
        2 * d[1] * d[11] * d[43] + 2 * d[3] * d[11] * d[43] + 2 * d[2] * d[12] * d[43] + 2 * d[6] * d[12] * d[43] -
        4 * d[0] * d[14] * d[43] + 2 * d[1] * d[15] * d[43] + 2 * d[3] * d[15] * d[43] - 4 * d[0] * d[16] * d[43] +
        8 * d[4] * d[9] * d[44] - 4 * d[1] * d[10] * d[44] - 4 * d[3] * d[10] * d[44] - 4 * d[1] * d[12] * d[44] -
        4 * d[3] * d[12] * d[44] + 8 * d[0] * d[13] * d[44];
    coeffs[57] =
        -2 * std::pow(d[14], 2) * d[36] - 4 * d[14] * d[16] * d[36] - 2 * std::pow(d[16], 2) * d[36] +
        8 * d[13] * d[17] * d[36] + 2 * d[11] * d[14] * d[37] + 2 * d[14] * d[15] * d[37] + 2 * d[11] * d[16] * d[37] +
        2 * d[15] * d[16] * d[37] - 4 * d[10] * d[17] * d[37] - 4 * d[12] * d[17] * d[37] - 4 * d[11] * d[13] * d[38] +
        2 * d[10] * d[14] * d[38] + 2 * d[12] * d[14] * d[38] - 4 * d[13] * d[15] * d[38] + 2 * d[10] * d[16] * d[38] +
        2 * d[12] * d[16] * d[38] + 2 * d[11] * d[14] * d[39] + 2 * d[14] * d[15] * d[39] + 2 * d[11] * d[16] * d[39] +
        2 * d[15] * d[16] * d[39] - 4 * d[10] * d[17] * d[39] - 4 * d[12] * d[17] * d[39] -
        2 * std::pow(d[11], 2) * d[40] - 4 * d[11] * d[15] * d[40] - 2 * std::pow(d[15], 2) * d[40] +
        8 * d[9] * d[17] * d[40] + 2 * d[10] * d[11] * d[41] + 2 * d[11] * d[12] * d[41] - 4 * d[9] * d[14] * d[41] +
        2 * d[10] * d[15] * d[41] + 2 * d[12] * d[15] * d[41] - 4 * d[9] * d[16] * d[41] - 4 * d[11] * d[13] * d[42] +
        2 * d[10] * d[14] * d[42] + 2 * d[12] * d[14] * d[42] - 4 * d[13] * d[15] * d[42] + 2 * d[10] * d[16] * d[42] +
        2 * d[12] * d[16] * d[42] + 2 * d[10] * d[11] * d[43] + 2 * d[11] * d[12] * d[43] - 4 * d[9] * d[14] * d[43] +
        2 * d[10] * d[15] * d[43] + 2 * d[12] * d[15] * d[43] - 4 * d[9] * d[16] * d[43] -
        2 * std::pow(d[10], 2) * d[44] - 4 * d[10] * d[12] * d[44] - 2 * std::pow(d[12], 2) * d[44] +
        8 * d[9] * d[13] * d[44];
    coeffs[58] =
        8 * d[8] * d[22] * d[36] - 4 * d[5] * d[23] * d[36] - 4 * d[7] * d[23] * d[36] - 4 * d[5] * d[25] * d[36] -
        4 * d[7] * d[25] * d[36] + 8 * d[4] * d[26] * d[36] - 4 * d[8] * d[19] * d[37] + 2 * d[5] * d[20] * d[37] +
        2 * d[7] * d[20] * d[37] - 4 * d[8] * d[21] * d[37] + 2 * d[2] * d[23] * d[37] + 2 * d[6] * d[23] * d[37] +
        2 * d[5] * d[24] * d[37] + 2 * d[7] * d[24] * d[37] + 2 * d[2] * d[25] * d[37] + 2 * d[6] * d[25] * d[37] -
        4 * d[1] * d[26] * d[37] - 4 * d[3] * d[26] * d[37] + 2 * d[5] * d[19] * d[38] + 2 * d[7] * d[19] * d[38] -
        4 * d[4] * d[20] * d[38] + 2 * d[5] * d[21] * d[38] + 2 * d[7] * d[21] * d[38] - 4 * d[2] * d[22] * d[38] -
        4 * d[6] * d[22] * d[38] + 2 * d[1] * d[23] * d[38] + 2 * d[3] * d[23] * d[38] - 4 * d[4] * d[24] * d[38] +
        2 * d[1] * d[25] * d[38] + 2 * d[3] * d[25] * d[38] - 4 * d[8] * d[19] * d[39] + 2 * d[5] * d[20] * d[39] +
        2 * d[7] * d[20] * d[39] - 4 * d[8] * d[21] * d[39] + 2 * d[2] * d[23] * d[39] + 2 * d[6] * d[23] * d[39] +
        2 * d[5] * d[24] * d[39] + 2 * d[7] * d[24] * d[39] + 2 * d[2] * d[25] * d[39] + 2 * d[6] * d[25] * d[39] -
        4 * d[1] * d[26] * d[39] - 4 * d[3] * d[26] * d[39] + 8 * d[8] * d[18] * d[40] - 4 * d[2] * d[20] * d[40] -
        4 * d[6] * d[20] * d[40] - 4 * d[2] * d[24] * d[40] - 4 * d[6] * d[24] * d[40] + 8 * d[0] * d[26] * d[40] -
        4 * d[5] * d[18] * d[41] - 4 * d[7] * d[18] * d[41] + 2 * d[2] * d[19] * d[41] + 2 * d[6] * d[19] * d[41] +
        2 * d[1] * d[20] * d[41] + 2 * d[3] * d[20] * d[41] + 2 * d[2] * d[21] * d[41] + 2 * d[6] * d[21] * d[41] -
        4 * d[0] * d[23] * d[41] + 2 * d[1] * d[24] * d[41] + 2 * d[3] * d[24] * d[41] - 4 * d[0] * d[25] * d[41] +
        2 * d[5] * d[19] * d[42] + 2 * d[7] * d[19] * d[42] - 4 * d[4] * d[20] * d[42] + 2 * d[5] * d[21] * d[42] +
        2 * d[7] * d[21] * d[42] - 4 * d[2] * d[22] * d[42] - 4 * d[6] * d[22] * d[42] + 2 * d[1] * d[23] * d[42] +
        2 * d[3] * d[23] * d[42] - 4 * d[4] * d[24] * d[42] + 2 * d[1] * d[25] * d[42] + 2 * d[3] * d[25] * d[42] -
        4 * d[5] * d[18] * d[43] - 4 * d[7] * d[18] * d[43] + 2 * d[2] * d[19] * d[43] + 2 * d[6] * d[19] * d[43] +
        2 * d[1] * d[20] * d[43] + 2 * d[3] * d[20] * d[43] + 2 * d[2] * d[21] * d[43] + 2 * d[6] * d[21] * d[43] -
        4 * d[0] * d[23] * d[43] + 2 * d[1] * d[24] * d[43] + 2 * d[3] * d[24] * d[43] - 4 * d[0] * d[25] * d[43] +
        8 * d[4] * d[18] * d[44] - 4 * d[1] * d[19] * d[44] - 4 * d[3] * d[19] * d[44] - 4 * d[1] * d[21] * d[44] -
        4 * d[3] * d[21] * d[44] + 8 * d[0] * d[22] * d[44];
    coeffs[59] =
        8 * d[17] * d[22] * d[36] - 4 * d[14] * d[23] * d[36] - 4 * d[16] * d[23] * d[36] - 4 * d[14] * d[25] * d[36] -
        4 * d[16] * d[25] * d[36] + 8 * d[13] * d[26] * d[36] - 4 * d[17] * d[19] * d[37] + 2 * d[14] * d[20] * d[37] +
        2 * d[16] * d[20] * d[37] - 4 * d[17] * d[21] * d[37] + 2 * d[11] * d[23] * d[37] + 2 * d[15] * d[23] * d[37] +
        2 * d[14] * d[24] * d[37] + 2 * d[16] * d[24] * d[37] + 2 * d[11] * d[25] * d[37] + 2 * d[15] * d[25] * d[37] -
        4 * d[10] * d[26] * d[37] - 4 * d[12] * d[26] * d[37] + 2 * d[14] * d[19] * d[38] + 2 * d[16] * d[19] * d[38] -
        4 * d[13] * d[20] * d[38] + 2 * d[14] * d[21] * d[38] + 2 * d[16] * d[21] * d[38] - 4 * d[11] * d[22] * d[38] -
        4 * d[15] * d[22] * d[38] + 2 * d[10] * d[23] * d[38] + 2 * d[12] * d[23] * d[38] - 4 * d[13] * d[24] * d[38] +
        2 * d[10] * d[25] * d[38] + 2 * d[12] * d[25] * d[38] - 4 * d[17] * d[19] * d[39] + 2 * d[14] * d[20] * d[39] +
        2 * d[16] * d[20] * d[39] - 4 * d[17] * d[21] * d[39] + 2 * d[11] * d[23] * d[39] + 2 * d[15] * d[23] * d[39] +
        2 * d[14] * d[24] * d[39] + 2 * d[16] * d[24] * d[39] + 2 * d[11] * d[25] * d[39] + 2 * d[15] * d[25] * d[39] -
        4 * d[10] * d[26] * d[39] - 4 * d[12] * d[26] * d[39] + 8 * d[17] * d[18] * d[40] - 4 * d[11] * d[20] * d[40] -
        4 * d[15] * d[20] * d[40] - 4 * d[11] * d[24] * d[40] - 4 * d[15] * d[24] * d[40] + 8 * d[9] * d[26] * d[40] -
        4 * d[14] * d[18] * d[41] - 4 * d[16] * d[18] * d[41] + 2 * d[11] * d[19] * d[41] + 2 * d[15] * d[19] * d[41] +
        2 * d[10] * d[20] * d[41] + 2 * d[12] * d[20] * d[41] + 2 * d[11] * d[21] * d[41] + 2 * d[15] * d[21] * d[41] -
        4 * d[9] * d[23] * d[41] + 2 * d[10] * d[24] * d[41] + 2 * d[12] * d[24] * d[41] - 4 * d[9] * d[25] * d[41] +
        2 * d[14] * d[19] * d[42] + 2 * d[16] * d[19] * d[42] - 4 * d[13] * d[20] * d[42] + 2 * d[14] * d[21] * d[42] +
        2 * d[16] * d[21] * d[42] - 4 * d[11] * d[22] * d[42] - 4 * d[15] * d[22] * d[42] + 2 * d[10] * d[23] * d[42] +
        2 * d[12] * d[23] * d[42] - 4 * d[13] * d[24] * d[42] + 2 * d[10] * d[25] * d[42] + 2 * d[12] * d[25] * d[42] -
        4 * d[14] * d[18] * d[43] - 4 * d[16] * d[18] * d[43] + 2 * d[11] * d[19] * d[43] + 2 * d[15] * d[19] * d[43] +
        2 * d[10] * d[20] * d[43] + 2 * d[12] * d[20] * d[43] + 2 * d[11] * d[21] * d[43] + 2 * d[15] * d[21] * d[43] -
        4 * d[9] * d[23] * d[43] + 2 * d[10] * d[24] * d[43] + 2 * d[12] * d[24] * d[43] - 4 * d[9] * d[25] * d[43] +
        8 * d[13] * d[18] * d[44] - 4 * d[10] * d[19] * d[44] - 4 * d[12] * d[19] * d[44] - 4 * d[10] * d[21] * d[44] -
        4 * d[12] * d[21] * d[44] + 8 * d[9] * d[22] * d[44];
    coeffs[60] =
        -2 * std::pow(d[23], 2) * d[36] - 4 * d[23] * d[25] * d[36] - 2 * std::pow(d[25], 2) * d[36] +
        8 * d[22] * d[26] * d[36] + 2 * d[20] * d[23] * d[37] + 2 * d[23] * d[24] * d[37] + 2 * d[20] * d[25] * d[37] +
        2 * d[24] * d[25] * d[37] - 4 * d[19] * d[26] * d[37] - 4 * d[21] * d[26] * d[37] - 4 * d[20] * d[22] * d[38] +
        2 * d[19] * d[23] * d[38] + 2 * d[21] * d[23] * d[38] - 4 * d[22] * d[24] * d[38] + 2 * d[19] * d[25] * d[38] +
        2 * d[21] * d[25] * d[38] + 2 * d[20] * d[23] * d[39] + 2 * d[23] * d[24] * d[39] + 2 * d[20] * d[25] * d[39] +
        2 * d[24] * d[25] * d[39] - 4 * d[19] * d[26] * d[39] - 4 * d[21] * d[26] * d[39] -
        2 * std::pow(d[20], 2) * d[40] - 4 * d[20] * d[24] * d[40] - 2 * std::pow(d[24], 2) * d[40] +
        8 * d[18] * d[26] * d[40] + 2 * d[19] * d[20] * d[41] + 2 * d[20] * d[21] * d[41] - 4 * d[18] * d[23] * d[41] +
        2 * d[19] * d[24] * d[41] + 2 * d[21] * d[24] * d[41] - 4 * d[18] * d[25] * d[41] - 4 * d[20] * d[22] * d[42] +
        2 * d[19] * d[23] * d[42] + 2 * d[21] * d[23] * d[42] - 4 * d[22] * d[24] * d[42] + 2 * d[19] * d[25] * d[42] +
        2 * d[21] * d[25] * d[42] + 2 * d[19] * d[20] * d[43] + 2 * d[20] * d[21] * d[43] - 4 * d[18] * d[23] * d[43] +
        2 * d[19] * d[24] * d[43] + 2 * d[21] * d[24] * d[43] - 4 * d[18] * d[25] * d[43] -
        2 * std::pow(d[19], 2) * d[44] - 4 * d[19] * d[21] * d[44] - 2 * std::pow(d[21], 2) * d[44] +
        8 * d[18] * d[22] * d[44];
    coeffs[61] =
        8 * d[8] * d[31] * d[36] - 4 * d[5] * d[32] * d[36] - 4 * d[7] * d[32] * d[36] - 4 * d[5] * d[34] * d[36] -
        4 * d[7] * d[34] * d[36] + 8 * d[4] * d[35] * d[36] - 4 * d[8] * d[28] * d[37] + 2 * d[5] * d[29] * d[37] +
        2 * d[7] * d[29] * d[37] - 4 * d[8] * d[30] * d[37] + 2 * d[2] * d[32] * d[37] + 2 * d[6] * d[32] * d[37] +
        2 * d[5] * d[33] * d[37] + 2 * d[7] * d[33] * d[37] + 2 * d[2] * d[34] * d[37] + 2 * d[6] * d[34] * d[37] -
        4 * d[1] * d[35] * d[37] - 4 * d[3] * d[35] * d[37] + 2 * d[5] * d[28] * d[38] + 2 * d[7] * d[28] * d[38] -
        4 * d[4] * d[29] * d[38] + 2 * d[5] * d[30] * d[38] + 2 * d[7] * d[30] * d[38] - 4 * d[2] * d[31] * d[38] -
        4 * d[6] * d[31] * d[38] + 2 * d[1] * d[32] * d[38] + 2 * d[3] * d[32] * d[38] - 4 * d[4] * d[33] * d[38] +
        2 * d[1] * d[34] * d[38] + 2 * d[3] * d[34] * d[38] - 4 * d[8] * d[28] * d[39] + 2 * d[5] * d[29] * d[39] +
        2 * d[7] * d[29] * d[39] - 4 * d[8] * d[30] * d[39] + 2 * d[2] * d[32] * d[39] + 2 * d[6] * d[32] * d[39] +
        2 * d[5] * d[33] * d[39] + 2 * d[7] * d[33] * d[39] + 2 * d[2] * d[34] * d[39] + 2 * d[6] * d[34] * d[39] -
        4 * d[1] * d[35] * d[39] - 4 * d[3] * d[35] * d[39] + 8 * d[8] * d[27] * d[40] - 4 * d[2] * d[29] * d[40] -
        4 * d[6] * d[29] * d[40] - 4 * d[2] * d[33] * d[40] - 4 * d[6] * d[33] * d[40] + 8 * d[0] * d[35] * d[40] -
        4 * d[5] * d[27] * d[41] - 4 * d[7] * d[27] * d[41] + 2 * d[2] * d[28] * d[41] + 2 * d[6] * d[28] * d[41] +
        2 * d[1] * d[29] * d[41] + 2 * d[3] * d[29] * d[41] + 2 * d[2] * d[30] * d[41] + 2 * d[6] * d[30] * d[41] -
        4 * d[0] * d[32] * d[41] + 2 * d[1] * d[33] * d[41] + 2 * d[3] * d[33] * d[41] - 4 * d[0] * d[34] * d[41] +
        2 * d[5] * d[28] * d[42] + 2 * d[7] * d[28] * d[42] - 4 * d[4] * d[29] * d[42] + 2 * d[5] * d[30] * d[42] +
        2 * d[7] * d[30] * d[42] - 4 * d[2] * d[31] * d[42] - 4 * d[6] * d[31] * d[42] + 2 * d[1] * d[32] * d[42] +
        2 * d[3] * d[32] * d[42] - 4 * d[4] * d[33] * d[42] + 2 * d[1] * d[34] * d[42] + 2 * d[3] * d[34] * d[42] -
        4 * d[5] * d[27] * d[43] - 4 * d[7] * d[27] * d[43] + 2 * d[2] * d[28] * d[43] + 2 * d[6] * d[28] * d[43] +
        2 * d[1] * d[29] * d[43] + 2 * d[3] * d[29] * d[43] + 2 * d[2] * d[30] * d[43] + 2 * d[6] * d[30] * d[43] -
        4 * d[0] * d[32] * d[43] + 2 * d[1] * d[33] * d[43] + 2 * d[3] * d[33] * d[43] - 4 * d[0] * d[34] * d[43] +
        8 * d[4] * d[27] * d[44] - 4 * d[1] * d[28] * d[44] - 4 * d[3] * d[28] * d[44] - 4 * d[1] * d[30] * d[44] -
        4 * d[3] * d[30] * d[44] + 8 * d[0] * d[31] * d[44];
    coeffs[62] =
        8 * d[17] * d[31] * d[36] - 4 * d[14] * d[32] * d[36] - 4 * d[16] * d[32] * d[36] - 4 * d[14] * d[34] * d[36] -
        4 * d[16] * d[34] * d[36] + 8 * d[13] * d[35] * d[36] - 4 * d[17] * d[28] * d[37] + 2 * d[14] * d[29] * d[37] +
        2 * d[16] * d[29] * d[37] - 4 * d[17] * d[30] * d[37] + 2 * d[11] * d[32] * d[37] + 2 * d[15] * d[32] * d[37] +
        2 * d[14] * d[33] * d[37] + 2 * d[16] * d[33] * d[37] + 2 * d[11] * d[34] * d[37] + 2 * d[15] * d[34] * d[37] -
        4 * d[10] * d[35] * d[37] - 4 * d[12] * d[35] * d[37] + 2 * d[14] * d[28] * d[38] + 2 * d[16] * d[28] * d[38] -
        4 * d[13] * d[29] * d[38] + 2 * d[14] * d[30] * d[38] + 2 * d[16] * d[30] * d[38] - 4 * d[11] * d[31] * d[38] -
        4 * d[15] * d[31] * d[38] + 2 * d[10] * d[32] * d[38] + 2 * d[12] * d[32] * d[38] - 4 * d[13] * d[33] * d[38] +
        2 * d[10] * d[34] * d[38] + 2 * d[12] * d[34] * d[38] - 4 * d[17] * d[28] * d[39] + 2 * d[14] * d[29] * d[39] +
        2 * d[16] * d[29] * d[39] - 4 * d[17] * d[30] * d[39] + 2 * d[11] * d[32] * d[39] + 2 * d[15] * d[32] * d[39] +
        2 * d[14] * d[33] * d[39] + 2 * d[16] * d[33] * d[39] + 2 * d[11] * d[34] * d[39] + 2 * d[15] * d[34] * d[39] -
        4 * d[10] * d[35] * d[39] - 4 * d[12] * d[35] * d[39] + 8 * d[17] * d[27] * d[40] - 4 * d[11] * d[29] * d[40] -
        4 * d[15] * d[29] * d[40] - 4 * d[11] * d[33] * d[40] - 4 * d[15] * d[33] * d[40] + 8 * d[9] * d[35] * d[40] -
        4 * d[14] * d[27] * d[41] - 4 * d[16] * d[27] * d[41] + 2 * d[11] * d[28] * d[41] + 2 * d[15] * d[28] * d[41] +
        2 * d[10] * d[29] * d[41] + 2 * d[12] * d[29] * d[41] + 2 * d[11] * d[30] * d[41] + 2 * d[15] * d[30] * d[41] -
        4 * d[9] * d[32] * d[41] + 2 * d[10] * d[33] * d[41] + 2 * d[12] * d[33] * d[41] - 4 * d[9] * d[34] * d[41] +
        2 * d[14] * d[28] * d[42] + 2 * d[16] * d[28] * d[42] - 4 * d[13] * d[29] * d[42] + 2 * d[14] * d[30] * d[42] +
        2 * d[16] * d[30] * d[42] - 4 * d[11] * d[31] * d[42] - 4 * d[15] * d[31] * d[42] + 2 * d[10] * d[32] * d[42] +
        2 * d[12] * d[32] * d[42] - 4 * d[13] * d[33] * d[42] + 2 * d[10] * d[34] * d[42] + 2 * d[12] * d[34] * d[42] -
        4 * d[14] * d[27] * d[43] - 4 * d[16] * d[27] * d[43] + 2 * d[11] * d[28] * d[43] + 2 * d[15] * d[28] * d[43] +
        2 * d[10] * d[29] * d[43] + 2 * d[12] * d[29] * d[43] + 2 * d[11] * d[30] * d[43] + 2 * d[15] * d[30] * d[43] -
        4 * d[9] * d[32] * d[43] + 2 * d[10] * d[33] * d[43] + 2 * d[12] * d[33] * d[43] - 4 * d[9] * d[34] * d[43] +
        8 * d[13] * d[27] * d[44] - 4 * d[10] * d[28] * d[44] - 4 * d[12] * d[28] * d[44] - 4 * d[10] * d[30] * d[44] -
        4 * d[12] * d[30] * d[44] + 8 * d[9] * d[31] * d[44];
    coeffs[63] =
        8 * d[26] * d[31] * d[36] - 4 * d[23] * d[32] * d[36] - 4 * d[25] * d[32] * d[36] - 4 * d[23] * d[34] * d[36] -
        4 * d[25] * d[34] * d[36] + 8 * d[22] * d[35] * d[36] - 4 * d[26] * d[28] * d[37] + 2 * d[23] * d[29] * d[37] +
        2 * d[25] * d[29] * d[37] - 4 * d[26] * d[30] * d[37] + 2 * d[20] * d[32] * d[37] + 2 * d[24] * d[32] * d[37] +
        2 * d[23] * d[33] * d[37] + 2 * d[25] * d[33] * d[37] + 2 * d[20] * d[34] * d[37] + 2 * d[24] * d[34] * d[37] -
        4 * d[19] * d[35] * d[37] - 4 * d[21] * d[35] * d[37] + 2 * d[23] * d[28] * d[38] + 2 * d[25] * d[28] * d[38] -
        4 * d[22] * d[29] * d[38] + 2 * d[23] * d[30] * d[38] + 2 * d[25] * d[30] * d[38] - 4 * d[20] * d[31] * d[38] -
        4 * d[24] * d[31] * d[38] + 2 * d[19] * d[32] * d[38] + 2 * d[21] * d[32] * d[38] - 4 * d[22] * d[33] * d[38] +
        2 * d[19] * d[34] * d[38] + 2 * d[21] * d[34] * d[38] - 4 * d[26] * d[28] * d[39] + 2 * d[23] * d[29] * d[39] +
        2 * d[25] * d[29] * d[39] - 4 * d[26] * d[30] * d[39] + 2 * d[20] * d[32] * d[39] + 2 * d[24] * d[32] * d[39] +
        2 * d[23] * d[33] * d[39] + 2 * d[25] * d[33] * d[39] + 2 * d[20] * d[34] * d[39] + 2 * d[24] * d[34] * d[39] -
        4 * d[19] * d[35] * d[39] - 4 * d[21] * d[35] * d[39] + 8 * d[26] * d[27] * d[40] - 4 * d[20] * d[29] * d[40] -
        4 * d[24] * d[29] * d[40] - 4 * d[20] * d[33] * d[40] - 4 * d[24] * d[33] * d[40] + 8 * d[18] * d[35] * d[40] -
        4 * d[23] * d[27] * d[41] - 4 * d[25] * d[27] * d[41] + 2 * d[20] * d[28] * d[41] + 2 * d[24] * d[28] * d[41] +
        2 * d[19] * d[29] * d[41] + 2 * d[21] * d[29] * d[41] + 2 * d[20] * d[30] * d[41] + 2 * d[24] * d[30] * d[41] -
        4 * d[18] * d[32] * d[41] + 2 * d[19] * d[33] * d[41] + 2 * d[21] * d[33] * d[41] - 4 * d[18] * d[34] * d[41] +
        2 * d[23] * d[28] * d[42] + 2 * d[25] * d[28] * d[42] - 4 * d[22] * d[29] * d[42] + 2 * d[23] * d[30] * d[42] +
        2 * d[25] * d[30] * d[42] - 4 * d[20] * d[31] * d[42] - 4 * d[24] * d[31] * d[42] + 2 * d[19] * d[32] * d[42] +
        2 * d[21] * d[32] * d[42] - 4 * d[22] * d[33] * d[42] + 2 * d[19] * d[34] * d[42] + 2 * d[21] * d[34] * d[42] -
        4 * d[23] * d[27] * d[43] - 4 * d[25] * d[27] * d[43] + 2 * d[20] * d[28] * d[43] + 2 * d[24] * d[28] * d[43] +
        2 * d[19] * d[29] * d[43] + 2 * d[21] * d[29] * d[43] + 2 * d[20] * d[30] * d[43] + 2 * d[24] * d[30] * d[43] -
        4 * d[18] * d[32] * d[43] + 2 * d[19] * d[33] * d[43] + 2 * d[21] * d[33] * d[43] - 4 * d[18] * d[34] * d[43] +
        8 * d[22] * d[27] * d[44] - 4 * d[19] * d[28] * d[44] - 4 * d[21] * d[28] * d[44] - 4 * d[19] * d[30] * d[44] -
        4 * d[21] * d[30] * d[44] + 8 * d[18] * d[31] * d[44];
    coeffs[64] =
        -2 * std::pow(d[32], 2) * d[36] - 4 * d[32] * d[34] * d[36] - 2 * std::pow(d[34], 2) * d[36] +
        8 * d[31] * d[35] * d[36] + 2 * d[29] * d[32] * d[37] + 2 * d[32] * d[33] * d[37] + 2 * d[29] * d[34] * d[37] +
        2 * d[33] * d[34] * d[37] - 4 * d[28] * d[35] * d[37] - 4 * d[30] * d[35] * d[37] - 4 * d[29] * d[31] * d[38] +
        2 * d[28] * d[32] * d[38] + 2 * d[30] * d[32] * d[38] - 4 * d[31] * d[33] * d[38] + 2 * d[28] * d[34] * d[38] +
        2 * d[30] * d[34] * d[38] + 2 * d[29] * d[32] * d[39] + 2 * d[32] * d[33] * d[39] + 2 * d[29] * d[34] * d[39] +
        2 * d[33] * d[34] * d[39] - 4 * d[28] * d[35] * d[39] - 4 * d[30] * d[35] * d[39] -
        2 * std::pow(d[29], 2) * d[40] - 4 * d[29] * d[33] * d[40] - 2 * std::pow(d[33], 2) * d[40] +
        8 * d[27] * d[35] * d[40] + 2 * d[28] * d[29] * d[41] + 2 * d[29] * d[30] * d[41] - 4 * d[27] * d[32] * d[41] +
        2 * d[28] * d[33] * d[41] + 2 * d[30] * d[33] * d[41] - 4 * d[27] * d[34] * d[41] - 4 * d[29] * d[31] * d[42] +
        2 * d[28] * d[32] * d[42] + 2 * d[30] * d[32] * d[42] - 4 * d[31] * d[33] * d[42] + 2 * d[28] * d[34] * d[42] +
        2 * d[30] * d[34] * d[42] + 2 * d[28] * d[29] * d[43] + 2 * d[29] * d[30] * d[43] - 4 * d[27] * d[32] * d[43] +
        2 * d[28] * d[33] * d[43] + 2 * d[30] * d[33] * d[43] - 4 * d[27] * d[34] * d[43] -
        2 * std::pow(d[28], 2) * d[44] - 4 * d[28] * d[30] * d[44] - 2 * std::pow(d[30], 2) * d[44] +
        8 * d[27] * d[31] * d[44];
    coeffs[65] =
        -2 * d[8] * std::pow(d[37], 2) + 2 * d[5] * d[37] * d[38] + 2 * d[7] * d[37] * d[38] -
        2 * d[4] * std::pow(d[38], 2) - 4 * d[8] * d[37] * d[39] + 2 * d[5] * d[38] * d[39] + 2 * d[7] * d[38] * d[39] -
        2 * d[8] * std::pow(d[39], 2) + 8 * d[8] * d[36] * d[40] - 4 * d[2] * d[38] * d[40] - 4 * d[6] * d[38] * d[40] -
        4 * d[5] * d[36] * d[41] - 4 * d[7] * d[36] * d[41] + 2 * d[2] * d[37] * d[41] + 2 * d[6] * d[37] * d[41] +
        2 * d[1] * d[38] * d[41] + 2 * d[3] * d[38] * d[41] + 2 * d[2] * d[39] * d[41] + 2 * d[6] * d[39] * d[41] -
        2 * d[0] * std::pow(d[41], 2) + 2 * d[5] * d[37] * d[42] + 2 * d[7] * d[37] * d[42] - 4 * d[4] * d[38] * d[42] +
        2 * d[5] * d[39] * d[42] + 2 * d[7] * d[39] * d[42] - 4 * d[2] * d[40] * d[42] - 4 * d[6] * d[40] * d[42] +
        2 * d[1] * d[41] * d[42] + 2 * d[3] * d[41] * d[42] - 2 * d[4] * std::pow(d[42], 2) - 4 * d[5] * d[36] * d[43] -
        4 * d[7] * d[36] * d[43] + 2 * d[2] * d[37] * d[43] + 2 * d[6] * d[37] * d[43] + 2 * d[1] * d[38] * d[43] +
        2 * d[3] * d[38] * d[43] + 2 * d[2] * d[39] * d[43] + 2 * d[6] * d[39] * d[43] - 4 * d[0] * d[41] * d[43] +
        2 * d[1] * d[42] * d[43] + 2 * d[3] * d[42] * d[43] - 2 * d[0] * std::pow(d[43], 2) + 8 * d[4] * d[36] * d[44] -
        4 * d[1] * d[37] * d[44] - 4 * d[3] * d[37] * d[44] - 4 * d[1] * d[39] * d[44] - 4 * d[3] * d[39] * d[44] +
        8 * d[0] * d[40] * d[44];
    coeffs[66] =
        -2 * d[17] * std::pow(d[37], 2) + 2 * d[14] * d[37] * d[38] + 2 * d[16] * d[37] * d[38] -
        2 * d[13] * std::pow(d[38], 2) - 4 * d[17] * d[37] * d[39] + 2 * d[14] * d[38] * d[39] +
        2 * d[16] * d[38] * d[39] - 2 * d[17] * std::pow(d[39], 2) + 8 * d[17] * d[36] * d[40] -
        4 * d[11] * d[38] * d[40] - 4 * d[15] * d[38] * d[40] - 4 * d[14] * d[36] * d[41] - 4 * d[16] * d[36] * d[41] +
        2 * d[11] * d[37] * d[41] + 2 * d[15] * d[37] * d[41] + 2 * d[10] * d[38] * d[41] + 2 * d[12] * d[38] * d[41] +
        2 * d[11] * d[39] * d[41] + 2 * d[15] * d[39] * d[41] - 2 * d[9] * std::pow(d[41], 2) +
        2 * d[14] * d[37] * d[42] + 2 * d[16] * d[37] * d[42] - 4 * d[13] * d[38] * d[42] + 2 * d[14] * d[39] * d[42] +
        2 * d[16] * d[39] * d[42] - 4 * d[11] * d[40] * d[42] - 4 * d[15] * d[40] * d[42] + 2 * d[10] * d[41] * d[42] +
        2 * d[12] * d[41] * d[42] - 2 * d[13] * std::pow(d[42], 2) - 4 * d[14] * d[36] * d[43] -
        4 * d[16] * d[36] * d[43] + 2 * d[11] * d[37] * d[43] + 2 * d[15] * d[37] * d[43] + 2 * d[10] * d[38] * d[43] +
        2 * d[12] * d[38] * d[43] + 2 * d[11] * d[39] * d[43] + 2 * d[15] * d[39] * d[43] - 4 * d[9] * d[41] * d[43] +
        2 * d[10] * d[42] * d[43] + 2 * d[12] * d[42] * d[43] - 2 * d[9] * std::pow(d[43], 2) +
        8 * d[13] * d[36] * d[44] - 4 * d[10] * d[37] * d[44] - 4 * d[12] * d[37] * d[44] - 4 * d[10] * d[39] * d[44] -
        4 * d[12] * d[39] * d[44] + 8 * d[9] * d[40] * d[44];
    coeffs[67] =
        -2 * d[26] * std::pow(d[37], 2) + 2 * d[23] * d[37] * d[38] + 2 * d[25] * d[37] * d[38] -
        2 * d[22] * std::pow(d[38], 2) - 4 * d[26] * d[37] * d[39] + 2 * d[23] * d[38] * d[39] +
        2 * d[25] * d[38] * d[39] - 2 * d[26] * std::pow(d[39], 2) + 8 * d[26] * d[36] * d[40] -
        4 * d[20] * d[38] * d[40] - 4 * d[24] * d[38] * d[40] - 4 * d[23] * d[36] * d[41] - 4 * d[25] * d[36] * d[41] +
        2 * d[20] * d[37] * d[41] + 2 * d[24] * d[37] * d[41] + 2 * d[19] * d[38] * d[41] + 2 * d[21] * d[38] * d[41] +
        2 * d[20] * d[39] * d[41] + 2 * d[24] * d[39] * d[41] - 2 * d[18] * std::pow(d[41], 2) +
        2 * d[23] * d[37] * d[42] + 2 * d[25] * d[37] * d[42] - 4 * d[22] * d[38] * d[42] + 2 * d[23] * d[39] * d[42] +
        2 * d[25] * d[39] * d[42] - 4 * d[20] * d[40] * d[42] - 4 * d[24] * d[40] * d[42] + 2 * d[19] * d[41] * d[42] +
        2 * d[21] * d[41] * d[42] - 2 * d[22] * std::pow(d[42], 2) - 4 * d[23] * d[36] * d[43] -
        4 * d[25] * d[36] * d[43] + 2 * d[20] * d[37] * d[43] + 2 * d[24] * d[37] * d[43] + 2 * d[19] * d[38] * d[43] +
        2 * d[21] * d[38] * d[43] + 2 * d[20] * d[39] * d[43] + 2 * d[24] * d[39] * d[43] - 4 * d[18] * d[41] * d[43] +
        2 * d[19] * d[42] * d[43] + 2 * d[21] * d[42] * d[43] - 2 * d[18] * std::pow(d[43], 2) +
        8 * d[22] * d[36] * d[44] - 4 * d[19] * d[37] * d[44] - 4 * d[21] * d[37] * d[44] - 4 * d[19] * d[39] * d[44] -
        4 * d[21] * d[39] * d[44] + 8 * d[18] * d[40] * d[44];
    coeffs[68] =
        -2 * d[35] * std::pow(d[37], 2) + 2 * d[32] * d[37] * d[38] + 2 * d[34] * d[37] * d[38] -
        2 * d[31] * std::pow(d[38], 2) - 4 * d[35] * d[37] * d[39] + 2 * d[32] * d[38] * d[39] +
        2 * d[34] * d[38] * d[39] - 2 * d[35] * std::pow(d[39], 2) + 8 * d[35] * d[36] * d[40] -
        4 * d[29] * d[38] * d[40] - 4 * d[33] * d[38] * d[40] - 4 * d[32] * d[36] * d[41] - 4 * d[34] * d[36] * d[41] +
        2 * d[29] * d[37] * d[41] + 2 * d[33] * d[37] * d[41] + 2 * d[28] * d[38] * d[41] + 2 * d[30] * d[38] * d[41] +
        2 * d[29] * d[39] * d[41] + 2 * d[33] * d[39] * d[41] - 2 * d[27] * std::pow(d[41], 2) +
        2 * d[32] * d[37] * d[42] + 2 * d[34] * d[37] * d[42] - 4 * d[31] * d[38] * d[42] + 2 * d[32] * d[39] * d[42] +
        2 * d[34] * d[39] * d[42] - 4 * d[29] * d[40] * d[42] - 4 * d[33] * d[40] * d[42] + 2 * d[28] * d[41] * d[42] +
        2 * d[30] * d[41] * d[42] - 2 * d[31] * std::pow(d[42], 2) - 4 * d[32] * d[36] * d[43] -
        4 * d[34] * d[36] * d[43] + 2 * d[29] * d[37] * d[43] + 2 * d[33] * d[37] * d[43] + 2 * d[28] * d[38] * d[43] +
        2 * d[30] * d[38] * d[43] + 2 * d[29] * d[39] * d[43] + 2 * d[33] * d[39] * d[43] - 4 * d[27] * d[41] * d[43] +
        2 * d[28] * d[42] * d[43] + 2 * d[30] * d[42] * d[43] - 2 * d[27] * std::pow(d[43], 2) +
        8 * d[31] * d[36] * d[44] - 4 * d[28] * d[37] * d[44] - 4 * d[30] * d[37] * d[44] - 4 * d[28] * d[39] * d[44] -
        4 * d[30] * d[39] * d[44] + 8 * d[27] * d[40] * d[44];
    coeffs[69] = -2 * std::pow(d[38], 2) * d[40] + 2 * d[37] * d[38] * d[41] + 2 * d[38] * d[39] * d[41] -
                 2 * d[36] * std::pow(d[41], 2) - 4 * d[38] * d[40] * d[42] + 2 * d[37] * d[41] * d[42] +
                 2 * d[39] * d[41] * d[42] - 2 * d[40] * std::pow(d[42], 2) + 2 * d[37] * d[38] * d[43] +
                 2 * d[38] * d[39] * d[43] - 4 * d[36] * d[41] * d[43] + 2 * d[37] * d[42] * d[43] +
                 2 * d[39] * d[42] * d[43] - 2 * d[36] * std::pow(d[43], 2) - 2 * std::pow(d[37], 2) * d[44] -
                 4 * d[37] * d[39] * d[44] - 2 * std::pow(d[39], 2) * d[44] + 8 * d[36] * d[40] * d[44];
    coeffs[70] = std::pow(d[0], 3) + d[0] * std::pow(d[1], 2) + d[0] * std::pow(d[2], 2) + d[0] * std::pow(d[3], 2) +
                 2 * d[1] * d[3] * d[4] - d[0] * std::pow(d[4], 2) + 2 * d[2] * d[3] * d[5] - d[0] * std::pow(d[5], 2) +
                 d[0] * std::pow(d[6], 2) + 2 * d[1] * d[6] * d[7] - d[0] * std::pow(d[7], 2) + 2 * d[2] * d[6] * d[8] -
                 d[0] * std::pow(d[8], 2);
    coeffs[71] = 3 * std::pow(d[0], 2) * d[9] + std::pow(d[1], 2) * d[9] + std::pow(d[2], 2) * d[9] +
                 std::pow(d[3], 2) * d[9] - std::pow(d[4], 2) * d[9] - std::pow(d[5], 2) * d[9] +
                 std::pow(d[6], 2) * d[9] - std::pow(d[7], 2) * d[9] - std::pow(d[8], 2) * d[9] +
                 2 * d[0] * d[1] * d[10] + 2 * d[3] * d[4] * d[10] + 2 * d[6] * d[7] * d[10] + 2 * d[0] * d[2] * d[11] +
                 2 * d[3] * d[5] * d[11] + 2 * d[6] * d[8] * d[11] + 2 * d[0] * d[3] * d[12] + 2 * d[1] * d[4] * d[12] +
                 2 * d[2] * d[5] * d[12] + 2 * d[1] * d[3] * d[13] - 2 * d[0] * d[4] * d[13] + 2 * d[2] * d[3] * d[14] -
                 2 * d[0] * d[5] * d[14] + 2 * d[0] * d[6] * d[15] + 2 * d[1] * d[7] * d[15] + 2 * d[2] * d[8] * d[15] +
                 2 * d[1] * d[6] * d[16] - 2 * d[0] * d[7] * d[16] + 2 * d[2] * d[6] * d[17] - 2 * d[0] * d[8] * d[17];
    coeffs[72] =
        3 * d[0] * std::pow(d[9], 2) + 2 * d[1] * d[9] * d[10] + d[0] * std::pow(d[10], 2) + 2 * d[2] * d[9] * d[11] +
        d[0] * std::pow(d[11], 2) + 2 * d[3] * d[9] * d[12] + 2 * d[4] * d[10] * d[12] + 2 * d[5] * d[11] * d[12] +
        d[0] * std::pow(d[12], 2) - 2 * d[4] * d[9] * d[13] + 2 * d[3] * d[10] * d[13] + 2 * d[1] * d[12] * d[13] -
        d[0] * std::pow(d[13], 2) - 2 * d[5] * d[9] * d[14] + 2 * d[3] * d[11] * d[14] + 2 * d[2] * d[12] * d[14] -
        d[0] * std::pow(d[14], 2) + 2 * d[6] * d[9] * d[15] + 2 * d[7] * d[10] * d[15] + 2 * d[8] * d[11] * d[15] +
        d[0] * std::pow(d[15], 2) - 2 * d[7] * d[9] * d[16] + 2 * d[6] * d[10] * d[16] + 2 * d[1] * d[15] * d[16] -
        d[0] * std::pow(d[16], 2) - 2 * d[8] * d[9] * d[17] + 2 * d[6] * d[11] * d[17] + 2 * d[2] * d[15] * d[17] -
        d[0] * std::pow(d[17], 2);
    coeffs[73] = std::pow(d[9], 3) + d[9] * std::pow(d[10], 2) + d[9] * std::pow(d[11], 2) + d[9] * std::pow(d[12], 2) +
                 2 * d[10] * d[12] * d[13] - d[9] * std::pow(d[13], 2) + 2 * d[11] * d[12] * d[14] -
                 d[9] * std::pow(d[14], 2) + d[9] * std::pow(d[15], 2) + 2 * d[10] * d[15] * d[16] -
                 d[9] * std::pow(d[16], 2) + 2 * d[11] * d[15] * d[17] - d[9] * std::pow(d[17], 2);
    coeffs[74] = 3 * std::pow(d[0], 2) * d[18] + std::pow(d[1], 2) * d[18] + std::pow(d[2], 2) * d[18] +
                 std::pow(d[3], 2) * d[18] - std::pow(d[4], 2) * d[18] - std::pow(d[5], 2) * d[18] +
                 std::pow(d[6], 2) * d[18] - std::pow(d[7], 2) * d[18] - std::pow(d[8], 2) * d[18] +
                 2 * d[0] * d[1] * d[19] + 2 * d[3] * d[4] * d[19] + 2 * d[6] * d[7] * d[19] + 2 * d[0] * d[2] * d[20] +
                 2 * d[3] * d[5] * d[20] + 2 * d[6] * d[8] * d[20] + 2 * d[0] * d[3] * d[21] + 2 * d[1] * d[4] * d[21] +
                 2 * d[2] * d[5] * d[21] + 2 * d[1] * d[3] * d[22] - 2 * d[0] * d[4] * d[22] + 2 * d[2] * d[3] * d[23] -
                 2 * d[0] * d[5] * d[23] + 2 * d[0] * d[6] * d[24] + 2 * d[1] * d[7] * d[24] + 2 * d[2] * d[8] * d[24] +
                 2 * d[1] * d[6] * d[25] - 2 * d[0] * d[7] * d[25] + 2 * d[2] * d[6] * d[26] - 2 * d[0] * d[8] * d[26];
    coeffs[75] =
        6 * d[0] * d[9] * d[18] + 2 * d[1] * d[10] * d[18] + 2 * d[2] * d[11] * d[18] + 2 * d[3] * d[12] * d[18] -
        2 * d[4] * d[13] * d[18] - 2 * d[5] * d[14] * d[18] + 2 * d[6] * d[15] * d[18] - 2 * d[7] * d[16] * d[18] -
        2 * d[8] * d[17] * d[18] + 2 * d[1] * d[9] * d[19] + 2 * d[0] * d[10] * d[19] + 2 * d[4] * d[12] * d[19] +
        2 * d[3] * d[13] * d[19] + 2 * d[7] * d[15] * d[19] + 2 * d[6] * d[16] * d[19] + 2 * d[2] * d[9] * d[20] +
        2 * d[0] * d[11] * d[20] + 2 * d[5] * d[12] * d[20] + 2 * d[3] * d[14] * d[20] + 2 * d[8] * d[15] * d[20] +
        2 * d[6] * d[17] * d[20] + 2 * d[3] * d[9] * d[21] + 2 * d[4] * d[10] * d[21] + 2 * d[5] * d[11] * d[21] +
        2 * d[0] * d[12] * d[21] + 2 * d[1] * d[13] * d[21] + 2 * d[2] * d[14] * d[21] - 2 * d[4] * d[9] * d[22] +
        2 * d[3] * d[10] * d[22] + 2 * d[1] * d[12] * d[22] - 2 * d[0] * d[13] * d[22] - 2 * d[5] * d[9] * d[23] +
        2 * d[3] * d[11] * d[23] + 2 * d[2] * d[12] * d[23] - 2 * d[0] * d[14] * d[23] + 2 * d[6] * d[9] * d[24] +
        2 * d[7] * d[10] * d[24] + 2 * d[8] * d[11] * d[24] + 2 * d[0] * d[15] * d[24] + 2 * d[1] * d[16] * d[24] +
        2 * d[2] * d[17] * d[24] - 2 * d[7] * d[9] * d[25] + 2 * d[6] * d[10] * d[25] + 2 * d[1] * d[15] * d[25] -
        2 * d[0] * d[16] * d[25] - 2 * d[8] * d[9] * d[26] + 2 * d[6] * d[11] * d[26] + 2 * d[2] * d[15] * d[26] -
        2 * d[0] * d[17] * d[26];
    coeffs[76] =
        3 * std::pow(d[9], 2) * d[18] + std::pow(d[10], 2) * d[18] + std::pow(d[11], 2) * d[18] +
        std::pow(d[12], 2) * d[18] - std::pow(d[13], 2) * d[18] - std::pow(d[14], 2) * d[18] +
        std::pow(d[15], 2) * d[18] - std::pow(d[16], 2) * d[18] - std::pow(d[17], 2) * d[18] +
        2 * d[9] * d[10] * d[19] + 2 * d[12] * d[13] * d[19] + 2 * d[15] * d[16] * d[19] + 2 * d[9] * d[11] * d[20] +
        2 * d[12] * d[14] * d[20] + 2 * d[15] * d[17] * d[20] + 2 * d[9] * d[12] * d[21] + 2 * d[10] * d[13] * d[21] +
        2 * d[11] * d[14] * d[21] + 2 * d[10] * d[12] * d[22] - 2 * d[9] * d[13] * d[22] + 2 * d[11] * d[12] * d[23] -
        2 * d[9] * d[14] * d[23] + 2 * d[9] * d[15] * d[24] + 2 * d[10] * d[16] * d[24] + 2 * d[11] * d[17] * d[24] +
        2 * d[10] * d[15] * d[25] - 2 * d[9] * d[16] * d[25] + 2 * d[11] * d[15] * d[26] - 2 * d[9] * d[17] * d[26];
    coeffs[77] =
        3 * d[0] * std::pow(d[18], 2) + 2 * d[1] * d[18] * d[19] + d[0] * std::pow(d[19], 2) +
        2 * d[2] * d[18] * d[20] + d[0] * std::pow(d[20], 2) + 2 * d[3] * d[18] * d[21] + 2 * d[4] * d[19] * d[21] +
        2 * d[5] * d[20] * d[21] + d[0] * std::pow(d[21], 2) - 2 * d[4] * d[18] * d[22] + 2 * d[3] * d[19] * d[22] +
        2 * d[1] * d[21] * d[22] - d[0] * std::pow(d[22], 2) - 2 * d[5] * d[18] * d[23] + 2 * d[3] * d[20] * d[23] +
        2 * d[2] * d[21] * d[23] - d[0] * std::pow(d[23], 2) + 2 * d[6] * d[18] * d[24] + 2 * d[7] * d[19] * d[24] +
        2 * d[8] * d[20] * d[24] + d[0] * std::pow(d[24], 2) - 2 * d[7] * d[18] * d[25] + 2 * d[6] * d[19] * d[25] +
        2 * d[1] * d[24] * d[25] - d[0] * std::pow(d[25], 2) - 2 * d[8] * d[18] * d[26] + 2 * d[6] * d[20] * d[26] +
        2 * d[2] * d[24] * d[26] - d[0] * std::pow(d[26], 2);
    coeffs[78] =
        3 * d[9] * std::pow(d[18], 2) + 2 * d[10] * d[18] * d[19] + d[9] * std::pow(d[19], 2) +
        2 * d[11] * d[18] * d[20] + d[9] * std::pow(d[20], 2) + 2 * d[12] * d[18] * d[21] + 2 * d[13] * d[19] * d[21] +
        2 * d[14] * d[20] * d[21] + d[9] * std::pow(d[21], 2) - 2 * d[13] * d[18] * d[22] + 2 * d[12] * d[19] * d[22] +
        2 * d[10] * d[21] * d[22] - d[9] * std::pow(d[22], 2) - 2 * d[14] * d[18] * d[23] + 2 * d[12] * d[20] * d[23] +
        2 * d[11] * d[21] * d[23] - d[9] * std::pow(d[23], 2) + 2 * d[15] * d[18] * d[24] + 2 * d[16] * d[19] * d[24] +
        2 * d[17] * d[20] * d[24] + d[9] * std::pow(d[24], 2) - 2 * d[16] * d[18] * d[25] + 2 * d[15] * d[19] * d[25] +
        2 * d[10] * d[24] * d[25] - d[9] * std::pow(d[25], 2) - 2 * d[17] * d[18] * d[26] + 2 * d[15] * d[20] * d[26] +
        2 * d[11] * d[24] * d[26] - d[9] * std::pow(d[26], 2);
    coeffs[79] = std::pow(d[18], 3) + d[18] * std::pow(d[19], 2) + d[18] * std::pow(d[20], 2) +
                 d[18] * std::pow(d[21], 2) + 2 * d[19] * d[21] * d[22] - d[18] * std::pow(d[22], 2) +
                 2 * d[20] * d[21] * d[23] - d[18] * std::pow(d[23], 2) + d[18] * std::pow(d[24], 2) +
                 2 * d[19] * d[24] * d[25] - d[18] * std::pow(d[25], 2) + 2 * d[20] * d[24] * d[26] -
                 d[18] * std::pow(d[26], 2);
    coeffs[80] = 3 * std::pow(d[0], 2) * d[27] + std::pow(d[1], 2) * d[27] + std::pow(d[2], 2) * d[27] +
                 std::pow(d[3], 2) * d[27] - std::pow(d[4], 2) * d[27] - std::pow(d[5], 2) * d[27] +
                 std::pow(d[6], 2) * d[27] - std::pow(d[7], 2) * d[27] - std::pow(d[8], 2) * d[27] +
                 2 * d[0] * d[1] * d[28] + 2 * d[3] * d[4] * d[28] + 2 * d[6] * d[7] * d[28] + 2 * d[0] * d[2] * d[29] +
                 2 * d[3] * d[5] * d[29] + 2 * d[6] * d[8] * d[29] + 2 * d[0] * d[3] * d[30] + 2 * d[1] * d[4] * d[30] +
                 2 * d[2] * d[5] * d[30] + 2 * d[1] * d[3] * d[31] - 2 * d[0] * d[4] * d[31] + 2 * d[2] * d[3] * d[32] -
                 2 * d[0] * d[5] * d[32] + 2 * d[0] * d[6] * d[33] + 2 * d[1] * d[7] * d[33] + 2 * d[2] * d[8] * d[33] +
                 2 * d[1] * d[6] * d[34] - 2 * d[0] * d[7] * d[34] + 2 * d[2] * d[6] * d[35] - 2 * d[0] * d[8] * d[35];
    coeffs[81] =
        6 * d[0] * d[9] * d[27] + 2 * d[1] * d[10] * d[27] + 2 * d[2] * d[11] * d[27] + 2 * d[3] * d[12] * d[27] -
        2 * d[4] * d[13] * d[27] - 2 * d[5] * d[14] * d[27] + 2 * d[6] * d[15] * d[27] - 2 * d[7] * d[16] * d[27] -
        2 * d[8] * d[17] * d[27] + 2 * d[1] * d[9] * d[28] + 2 * d[0] * d[10] * d[28] + 2 * d[4] * d[12] * d[28] +
        2 * d[3] * d[13] * d[28] + 2 * d[7] * d[15] * d[28] + 2 * d[6] * d[16] * d[28] + 2 * d[2] * d[9] * d[29] +
        2 * d[0] * d[11] * d[29] + 2 * d[5] * d[12] * d[29] + 2 * d[3] * d[14] * d[29] + 2 * d[8] * d[15] * d[29] +
        2 * d[6] * d[17] * d[29] + 2 * d[3] * d[9] * d[30] + 2 * d[4] * d[10] * d[30] + 2 * d[5] * d[11] * d[30] +
        2 * d[0] * d[12] * d[30] + 2 * d[1] * d[13] * d[30] + 2 * d[2] * d[14] * d[30] - 2 * d[4] * d[9] * d[31] +
        2 * d[3] * d[10] * d[31] + 2 * d[1] * d[12] * d[31] - 2 * d[0] * d[13] * d[31] - 2 * d[5] * d[9] * d[32] +
        2 * d[3] * d[11] * d[32] + 2 * d[2] * d[12] * d[32] - 2 * d[0] * d[14] * d[32] + 2 * d[6] * d[9] * d[33] +
        2 * d[7] * d[10] * d[33] + 2 * d[8] * d[11] * d[33] + 2 * d[0] * d[15] * d[33] + 2 * d[1] * d[16] * d[33] +
        2 * d[2] * d[17] * d[33] - 2 * d[7] * d[9] * d[34] + 2 * d[6] * d[10] * d[34] + 2 * d[1] * d[15] * d[34] -
        2 * d[0] * d[16] * d[34] - 2 * d[8] * d[9] * d[35] + 2 * d[6] * d[11] * d[35] + 2 * d[2] * d[15] * d[35] -
        2 * d[0] * d[17] * d[35];
    coeffs[82] =
        3 * std::pow(d[9], 2) * d[27] + std::pow(d[10], 2) * d[27] + std::pow(d[11], 2) * d[27] +
        std::pow(d[12], 2) * d[27] - std::pow(d[13], 2) * d[27] - std::pow(d[14], 2) * d[27] +
        std::pow(d[15], 2) * d[27] - std::pow(d[16], 2) * d[27] - std::pow(d[17], 2) * d[27] +
        2 * d[9] * d[10] * d[28] + 2 * d[12] * d[13] * d[28] + 2 * d[15] * d[16] * d[28] + 2 * d[9] * d[11] * d[29] +
        2 * d[12] * d[14] * d[29] + 2 * d[15] * d[17] * d[29] + 2 * d[9] * d[12] * d[30] + 2 * d[10] * d[13] * d[30] +
        2 * d[11] * d[14] * d[30] + 2 * d[10] * d[12] * d[31] - 2 * d[9] * d[13] * d[31] + 2 * d[11] * d[12] * d[32] -
        2 * d[9] * d[14] * d[32] + 2 * d[9] * d[15] * d[33] + 2 * d[10] * d[16] * d[33] + 2 * d[11] * d[17] * d[33] +
        2 * d[10] * d[15] * d[34] - 2 * d[9] * d[16] * d[34] + 2 * d[11] * d[15] * d[35] - 2 * d[9] * d[17] * d[35];
    coeffs[83] =
        6 * d[0] * d[18] * d[27] + 2 * d[1] * d[19] * d[27] + 2 * d[2] * d[20] * d[27] + 2 * d[3] * d[21] * d[27] -
        2 * d[4] * d[22] * d[27] - 2 * d[5] * d[23] * d[27] + 2 * d[6] * d[24] * d[27] - 2 * d[7] * d[25] * d[27] -
        2 * d[8] * d[26] * d[27] + 2 * d[1] * d[18] * d[28] + 2 * d[0] * d[19] * d[28] + 2 * d[4] * d[21] * d[28] +
        2 * d[3] * d[22] * d[28] + 2 * d[7] * d[24] * d[28] + 2 * d[6] * d[25] * d[28] + 2 * d[2] * d[18] * d[29] +
        2 * d[0] * d[20] * d[29] + 2 * d[5] * d[21] * d[29] + 2 * d[3] * d[23] * d[29] + 2 * d[8] * d[24] * d[29] +
        2 * d[6] * d[26] * d[29] + 2 * d[3] * d[18] * d[30] + 2 * d[4] * d[19] * d[30] + 2 * d[5] * d[20] * d[30] +
        2 * d[0] * d[21] * d[30] + 2 * d[1] * d[22] * d[30] + 2 * d[2] * d[23] * d[30] - 2 * d[4] * d[18] * d[31] +
        2 * d[3] * d[19] * d[31] + 2 * d[1] * d[21] * d[31] - 2 * d[0] * d[22] * d[31] - 2 * d[5] * d[18] * d[32] +
        2 * d[3] * d[20] * d[32] + 2 * d[2] * d[21] * d[32] - 2 * d[0] * d[23] * d[32] + 2 * d[6] * d[18] * d[33] +
        2 * d[7] * d[19] * d[33] + 2 * d[8] * d[20] * d[33] + 2 * d[0] * d[24] * d[33] + 2 * d[1] * d[25] * d[33] +
        2 * d[2] * d[26] * d[33] - 2 * d[7] * d[18] * d[34] + 2 * d[6] * d[19] * d[34] + 2 * d[1] * d[24] * d[34] -
        2 * d[0] * d[25] * d[34] - 2 * d[8] * d[18] * d[35] + 2 * d[6] * d[20] * d[35] + 2 * d[2] * d[24] * d[35] -
        2 * d[0] * d[26] * d[35];
    coeffs[84] =
        6 * d[9] * d[18] * d[27] + 2 * d[10] * d[19] * d[27] + 2 * d[11] * d[20] * d[27] + 2 * d[12] * d[21] * d[27] -
        2 * d[13] * d[22] * d[27] - 2 * d[14] * d[23] * d[27] + 2 * d[15] * d[24] * d[27] - 2 * d[16] * d[25] * d[27] -
        2 * d[17] * d[26] * d[27] + 2 * d[10] * d[18] * d[28] + 2 * d[9] * d[19] * d[28] + 2 * d[13] * d[21] * d[28] +
        2 * d[12] * d[22] * d[28] + 2 * d[16] * d[24] * d[28] + 2 * d[15] * d[25] * d[28] + 2 * d[11] * d[18] * d[29] +
        2 * d[9] * d[20] * d[29] + 2 * d[14] * d[21] * d[29] + 2 * d[12] * d[23] * d[29] + 2 * d[17] * d[24] * d[29] +
        2 * d[15] * d[26] * d[29] + 2 * d[12] * d[18] * d[30] + 2 * d[13] * d[19] * d[30] + 2 * d[14] * d[20] * d[30] +
        2 * d[9] * d[21] * d[30] + 2 * d[10] * d[22] * d[30] + 2 * d[11] * d[23] * d[30] - 2 * d[13] * d[18] * d[31] +
        2 * d[12] * d[19] * d[31] + 2 * d[10] * d[21] * d[31] - 2 * d[9] * d[22] * d[31] - 2 * d[14] * d[18] * d[32] +
        2 * d[12] * d[20] * d[32] + 2 * d[11] * d[21] * d[32] - 2 * d[9] * d[23] * d[32] + 2 * d[15] * d[18] * d[33] +
        2 * d[16] * d[19] * d[33] + 2 * d[17] * d[20] * d[33] + 2 * d[9] * d[24] * d[33] + 2 * d[10] * d[25] * d[33] +
        2 * d[11] * d[26] * d[33] - 2 * d[16] * d[18] * d[34] + 2 * d[15] * d[19] * d[34] + 2 * d[10] * d[24] * d[34] -
        2 * d[9] * d[25] * d[34] - 2 * d[17] * d[18] * d[35] + 2 * d[15] * d[20] * d[35] + 2 * d[11] * d[24] * d[35] -
        2 * d[9] * d[26] * d[35];
    coeffs[85] =
        3 * std::pow(d[18], 2) * d[27] + std::pow(d[19], 2) * d[27] + std::pow(d[20], 2) * d[27] +
        std::pow(d[21], 2) * d[27] - std::pow(d[22], 2) * d[27] - std::pow(d[23], 2) * d[27] +
        std::pow(d[24], 2) * d[27] - std::pow(d[25], 2) * d[27] - std::pow(d[26], 2) * d[27] +
        2 * d[18] * d[19] * d[28] + 2 * d[21] * d[22] * d[28] + 2 * d[24] * d[25] * d[28] + 2 * d[18] * d[20] * d[29] +
        2 * d[21] * d[23] * d[29] + 2 * d[24] * d[26] * d[29] + 2 * d[18] * d[21] * d[30] + 2 * d[19] * d[22] * d[30] +
        2 * d[20] * d[23] * d[30] + 2 * d[19] * d[21] * d[31] - 2 * d[18] * d[22] * d[31] + 2 * d[20] * d[21] * d[32] -
        2 * d[18] * d[23] * d[32] + 2 * d[18] * d[24] * d[33] + 2 * d[19] * d[25] * d[33] + 2 * d[20] * d[26] * d[33] +
        2 * d[19] * d[24] * d[34] - 2 * d[18] * d[25] * d[34] + 2 * d[20] * d[24] * d[35] - 2 * d[18] * d[26] * d[35];
    coeffs[86] =
        3 * d[0] * std::pow(d[27], 2) + 2 * d[1] * d[27] * d[28] + d[0] * std::pow(d[28], 2) +
        2 * d[2] * d[27] * d[29] + d[0] * std::pow(d[29], 2) + 2 * d[3] * d[27] * d[30] + 2 * d[4] * d[28] * d[30] +
        2 * d[5] * d[29] * d[30] + d[0] * std::pow(d[30], 2) - 2 * d[4] * d[27] * d[31] + 2 * d[3] * d[28] * d[31] +
        2 * d[1] * d[30] * d[31] - d[0] * std::pow(d[31], 2) - 2 * d[5] * d[27] * d[32] + 2 * d[3] * d[29] * d[32] +
        2 * d[2] * d[30] * d[32] - d[0] * std::pow(d[32], 2) + 2 * d[6] * d[27] * d[33] + 2 * d[7] * d[28] * d[33] +
        2 * d[8] * d[29] * d[33] + d[0] * std::pow(d[33], 2) - 2 * d[7] * d[27] * d[34] + 2 * d[6] * d[28] * d[34] +
        2 * d[1] * d[33] * d[34] - d[0] * std::pow(d[34], 2) - 2 * d[8] * d[27] * d[35] + 2 * d[6] * d[29] * d[35] +
        2 * d[2] * d[33] * d[35] - d[0] * std::pow(d[35], 2);
    coeffs[87] =
        3 * d[9] * std::pow(d[27], 2) + 2 * d[10] * d[27] * d[28] + d[9] * std::pow(d[28], 2) +
        2 * d[11] * d[27] * d[29] + d[9] * std::pow(d[29], 2) + 2 * d[12] * d[27] * d[30] + 2 * d[13] * d[28] * d[30] +
        2 * d[14] * d[29] * d[30] + d[9] * std::pow(d[30], 2) - 2 * d[13] * d[27] * d[31] + 2 * d[12] * d[28] * d[31] +
        2 * d[10] * d[30] * d[31] - d[9] * std::pow(d[31], 2) - 2 * d[14] * d[27] * d[32] + 2 * d[12] * d[29] * d[32] +
        2 * d[11] * d[30] * d[32] - d[9] * std::pow(d[32], 2) + 2 * d[15] * d[27] * d[33] + 2 * d[16] * d[28] * d[33] +
        2 * d[17] * d[29] * d[33] + d[9] * std::pow(d[33], 2) - 2 * d[16] * d[27] * d[34] + 2 * d[15] * d[28] * d[34] +
        2 * d[10] * d[33] * d[34] - d[9] * std::pow(d[34], 2) - 2 * d[17] * d[27] * d[35] + 2 * d[15] * d[29] * d[35] +
        2 * d[11] * d[33] * d[35] - d[9] * std::pow(d[35], 2);
    coeffs[88] =
        3 * d[18] * std::pow(d[27], 2) + 2 * d[19] * d[27] * d[28] + d[18] * std::pow(d[28], 2) +
        2 * d[20] * d[27] * d[29] + d[18] * std::pow(d[29], 2) + 2 * d[21] * d[27] * d[30] + 2 * d[22] * d[28] * d[30] +
        2 * d[23] * d[29] * d[30] + d[18] * std::pow(d[30], 2) - 2 * d[22] * d[27] * d[31] + 2 * d[21] * d[28] * d[31] +
        2 * d[19] * d[30] * d[31] - d[18] * std::pow(d[31], 2) - 2 * d[23] * d[27] * d[32] + 2 * d[21] * d[29] * d[32] +
        2 * d[20] * d[30] * d[32] - d[18] * std::pow(d[32], 2) + 2 * d[24] * d[27] * d[33] + 2 * d[25] * d[28] * d[33] +
        2 * d[26] * d[29] * d[33] + d[18] * std::pow(d[33], 2) - 2 * d[25] * d[27] * d[34] + 2 * d[24] * d[28] * d[34] +
        2 * d[19] * d[33] * d[34] - d[18] * std::pow(d[34], 2) - 2 * d[26] * d[27] * d[35] + 2 * d[24] * d[29] * d[35] +
        2 * d[20] * d[33] * d[35] - d[18] * std::pow(d[35], 2);
    coeffs[89] = std::pow(d[27], 3) + d[27] * std::pow(d[28], 2) + d[27] * std::pow(d[29], 2) +
                 d[27] * std::pow(d[30], 2) + 2 * d[28] * d[30] * d[31] - d[27] * std::pow(d[31], 2) +
                 2 * d[29] * d[30] * d[32] - d[27] * std::pow(d[32], 2) + d[27] * std::pow(d[33], 2) +
                 2 * d[28] * d[33] * d[34] - d[27] * std::pow(d[34], 2) + 2 * d[29] * d[33] * d[35] -
                 d[27] * std::pow(d[35], 2);
    coeffs[90] = 3 * std::pow(d[0], 2) * d[36] + std::pow(d[1], 2) * d[36] + std::pow(d[2], 2) * d[36] +
                 std::pow(d[3], 2) * d[36] - std::pow(d[4], 2) * d[36] - std::pow(d[5], 2) * d[36] +
                 std::pow(d[6], 2) * d[36] - std::pow(d[7], 2) * d[36] - std::pow(d[8], 2) * d[36] +
                 2 * d[0] * d[1] * d[37] + 2 * d[3] * d[4] * d[37] + 2 * d[6] * d[7] * d[37] + 2 * d[0] * d[2] * d[38] +
                 2 * d[3] * d[5] * d[38] + 2 * d[6] * d[8] * d[38] + 2 * d[0] * d[3] * d[39] + 2 * d[1] * d[4] * d[39] +
                 2 * d[2] * d[5] * d[39] + 2 * d[1] * d[3] * d[40] - 2 * d[0] * d[4] * d[40] + 2 * d[2] * d[3] * d[41] -
                 2 * d[0] * d[5] * d[41] + 2 * d[0] * d[6] * d[42] + 2 * d[1] * d[7] * d[42] + 2 * d[2] * d[8] * d[42] +
                 2 * d[1] * d[6] * d[43] - 2 * d[0] * d[7] * d[43] + 2 * d[2] * d[6] * d[44] - 2 * d[0] * d[8] * d[44];
    coeffs[91] =
        6 * d[0] * d[9] * d[36] + 2 * d[1] * d[10] * d[36] + 2 * d[2] * d[11] * d[36] + 2 * d[3] * d[12] * d[36] -
        2 * d[4] * d[13] * d[36] - 2 * d[5] * d[14] * d[36] + 2 * d[6] * d[15] * d[36] - 2 * d[7] * d[16] * d[36] -
        2 * d[8] * d[17] * d[36] + 2 * d[1] * d[9] * d[37] + 2 * d[0] * d[10] * d[37] + 2 * d[4] * d[12] * d[37] +
        2 * d[3] * d[13] * d[37] + 2 * d[7] * d[15] * d[37] + 2 * d[6] * d[16] * d[37] + 2 * d[2] * d[9] * d[38] +
        2 * d[0] * d[11] * d[38] + 2 * d[5] * d[12] * d[38] + 2 * d[3] * d[14] * d[38] + 2 * d[8] * d[15] * d[38] +
        2 * d[6] * d[17] * d[38] + 2 * d[3] * d[9] * d[39] + 2 * d[4] * d[10] * d[39] + 2 * d[5] * d[11] * d[39] +
        2 * d[0] * d[12] * d[39] + 2 * d[1] * d[13] * d[39] + 2 * d[2] * d[14] * d[39] - 2 * d[4] * d[9] * d[40] +
        2 * d[3] * d[10] * d[40] + 2 * d[1] * d[12] * d[40] - 2 * d[0] * d[13] * d[40] - 2 * d[5] * d[9] * d[41] +
        2 * d[3] * d[11] * d[41] + 2 * d[2] * d[12] * d[41] - 2 * d[0] * d[14] * d[41] + 2 * d[6] * d[9] * d[42] +
        2 * d[7] * d[10] * d[42] + 2 * d[8] * d[11] * d[42] + 2 * d[0] * d[15] * d[42] + 2 * d[1] * d[16] * d[42] +
        2 * d[2] * d[17] * d[42] - 2 * d[7] * d[9] * d[43] + 2 * d[6] * d[10] * d[43] + 2 * d[1] * d[15] * d[43] -
        2 * d[0] * d[16] * d[43] - 2 * d[8] * d[9] * d[44] + 2 * d[6] * d[11] * d[44] + 2 * d[2] * d[15] * d[44] -
        2 * d[0] * d[17] * d[44];
    coeffs[92] =
        3 * std::pow(d[9], 2) * d[36] + std::pow(d[10], 2) * d[36] + std::pow(d[11], 2) * d[36] +
        std::pow(d[12], 2) * d[36] - std::pow(d[13], 2) * d[36] - std::pow(d[14], 2) * d[36] +
        std::pow(d[15], 2) * d[36] - std::pow(d[16], 2) * d[36] - std::pow(d[17], 2) * d[36] +
        2 * d[9] * d[10] * d[37] + 2 * d[12] * d[13] * d[37] + 2 * d[15] * d[16] * d[37] + 2 * d[9] * d[11] * d[38] +
        2 * d[12] * d[14] * d[38] + 2 * d[15] * d[17] * d[38] + 2 * d[9] * d[12] * d[39] + 2 * d[10] * d[13] * d[39] +
        2 * d[11] * d[14] * d[39] + 2 * d[10] * d[12] * d[40] - 2 * d[9] * d[13] * d[40] + 2 * d[11] * d[12] * d[41] -
        2 * d[9] * d[14] * d[41] + 2 * d[9] * d[15] * d[42] + 2 * d[10] * d[16] * d[42] + 2 * d[11] * d[17] * d[42] +
        2 * d[10] * d[15] * d[43] - 2 * d[9] * d[16] * d[43] + 2 * d[11] * d[15] * d[44] - 2 * d[9] * d[17] * d[44];
    coeffs[93] =
        6 * d[0] * d[18] * d[36] + 2 * d[1] * d[19] * d[36] + 2 * d[2] * d[20] * d[36] + 2 * d[3] * d[21] * d[36] -
        2 * d[4] * d[22] * d[36] - 2 * d[5] * d[23] * d[36] + 2 * d[6] * d[24] * d[36] - 2 * d[7] * d[25] * d[36] -
        2 * d[8] * d[26] * d[36] + 2 * d[1] * d[18] * d[37] + 2 * d[0] * d[19] * d[37] + 2 * d[4] * d[21] * d[37] +
        2 * d[3] * d[22] * d[37] + 2 * d[7] * d[24] * d[37] + 2 * d[6] * d[25] * d[37] + 2 * d[2] * d[18] * d[38] +
        2 * d[0] * d[20] * d[38] + 2 * d[5] * d[21] * d[38] + 2 * d[3] * d[23] * d[38] + 2 * d[8] * d[24] * d[38] +
        2 * d[6] * d[26] * d[38] + 2 * d[3] * d[18] * d[39] + 2 * d[4] * d[19] * d[39] + 2 * d[5] * d[20] * d[39] +
        2 * d[0] * d[21] * d[39] + 2 * d[1] * d[22] * d[39] + 2 * d[2] * d[23] * d[39] - 2 * d[4] * d[18] * d[40] +
        2 * d[3] * d[19] * d[40] + 2 * d[1] * d[21] * d[40] - 2 * d[0] * d[22] * d[40] - 2 * d[5] * d[18] * d[41] +
        2 * d[3] * d[20] * d[41] + 2 * d[2] * d[21] * d[41] - 2 * d[0] * d[23] * d[41] + 2 * d[6] * d[18] * d[42] +
        2 * d[7] * d[19] * d[42] + 2 * d[8] * d[20] * d[42] + 2 * d[0] * d[24] * d[42] + 2 * d[1] * d[25] * d[42] +
        2 * d[2] * d[26] * d[42] - 2 * d[7] * d[18] * d[43] + 2 * d[6] * d[19] * d[43] + 2 * d[1] * d[24] * d[43] -
        2 * d[0] * d[25] * d[43] - 2 * d[8] * d[18] * d[44] + 2 * d[6] * d[20] * d[44] + 2 * d[2] * d[24] * d[44] -
        2 * d[0] * d[26] * d[44];
    coeffs[94] =
        6 * d[9] * d[18] * d[36] + 2 * d[10] * d[19] * d[36] + 2 * d[11] * d[20] * d[36] + 2 * d[12] * d[21] * d[36] -
        2 * d[13] * d[22] * d[36] - 2 * d[14] * d[23] * d[36] + 2 * d[15] * d[24] * d[36] - 2 * d[16] * d[25] * d[36] -
        2 * d[17] * d[26] * d[36] + 2 * d[10] * d[18] * d[37] + 2 * d[9] * d[19] * d[37] + 2 * d[13] * d[21] * d[37] +
        2 * d[12] * d[22] * d[37] + 2 * d[16] * d[24] * d[37] + 2 * d[15] * d[25] * d[37] + 2 * d[11] * d[18] * d[38] +
        2 * d[9] * d[20] * d[38] + 2 * d[14] * d[21] * d[38] + 2 * d[12] * d[23] * d[38] + 2 * d[17] * d[24] * d[38] +
        2 * d[15] * d[26] * d[38] + 2 * d[12] * d[18] * d[39] + 2 * d[13] * d[19] * d[39] + 2 * d[14] * d[20] * d[39] +
        2 * d[9] * d[21] * d[39] + 2 * d[10] * d[22] * d[39] + 2 * d[11] * d[23] * d[39] - 2 * d[13] * d[18] * d[40] +
        2 * d[12] * d[19] * d[40] + 2 * d[10] * d[21] * d[40] - 2 * d[9] * d[22] * d[40] - 2 * d[14] * d[18] * d[41] +
        2 * d[12] * d[20] * d[41] + 2 * d[11] * d[21] * d[41] - 2 * d[9] * d[23] * d[41] + 2 * d[15] * d[18] * d[42] +
        2 * d[16] * d[19] * d[42] + 2 * d[17] * d[20] * d[42] + 2 * d[9] * d[24] * d[42] + 2 * d[10] * d[25] * d[42] +
        2 * d[11] * d[26] * d[42] - 2 * d[16] * d[18] * d[43] + 2 * d[15] * d[19] * d[43] + 2 * d[10] * d[24] * d[43] -
        2 * d[9] * d[25] * d[43] - 2 * d[17] * d[18] * d[44] + 2 * d[15] * d[20] * d[44] + 2 * d[11] * d[24] * d[44] -
        2 * d[9] * d[26] * d[44];
    coeffs[95] =
        3 * std::pow(d[18], 2) * d[36] + std::pow(d[19], 2) * d[36] + std::pow(d[20], 2) * d[36] +
        std::pow(d[21], 2) * d[36] - std::pow(d[22], 2) * d[36] - std::pow(d[23], 2) * d[36] +
        std::pow(d[24], 2) * d[36] - std::pow(d[25], 2) * d[36] - std::pow(d[26], 2) * d[36] +
        2 * d[18] * d[19] * d[37] + 2 * d[21] * d[22] * d[37] + 2 * d[24] * d[25] * d[37] + 2 * d[18] * d[20] * d[38] +
        2 * d[21] * d[23] * d[38] + 2 * d[24] * d[26] * d[38] + 2 * d[18] * d[21] * d[39] + 2 * d[19] * d[22] * d[39] +
        2 * d[20] * d[23] * d[39] + 2 * d[19] * d[21] * d[40] - 2 * d[18] * d[22] * d[40] + 2 * d[20] * d[21] * d[41] -
        2 * d[18] * d[23] * d[41] + 2 * d[18] * d[24] * d[42] + 2 * d[19] * d[25] * d[42] + 2 * d[20] * d[26] * d[42] +
        2 * d[19] * d[24] * d[43] - 2 * d[18] * d[25] * d[43] + 2 * d[20] * d[24] * d[44] - 2 * d[18] * d[26] * d[44];
    coeffs[96] =
        6 * d[0] * d[27] * d[36] + 2 * d[1] * d[28] * d[36] + 2 * d[2] * d[29] * d[36] + 2 * d[3] * d[30] * d[36] -
        2 * d[4] * d[31] * d[36] - 2 * d[5] * d[32] * d[36] + 2 * d[6] * d[33] * d[36] - 2 * d[7] * d[34] * d[36] -
        2 * d[8] * d[35] * d[36] + 2 * d[1] * d[27] * d[37] + 2 * d[0] * d[28] * d[37] + 2 * d[4] * d[30] * d[37] +
        2 * d[3] * d[31] * d[37] + 2 * d[7] * d[33] * d[37] + 2 * d[6] * d[34] * d[37] + 2 * d[2] * d[27] * d[38] +
        2 * d[0] * d[29] * d[38] + 2 * d[5] * d[30] * d[38] + 2 * d[3] * d[32] * d[38] + 2 * d[8] * d[33] * d[38] +
        2 * d[6] * d[35] * d[38] + 2 * d[3] * d[27] * d[39] + 2 * d[4] * d[28] * d[39] + 2 * d[5] * d[29] * d[39] +
        2 * d[0] * d[30] * d[39] + 2 * d[1] * d[31] * d[39] + 2 * d[2] * d[32] * d[39] - 2 * d[4] * d[27] * d[40] +
        2 * d[3] * d[28] * d[40] + 2 * d[1] * d[30] * d[40] - 2 * d[0] * d[31] * d[40] - 2 * d[5] * d[27] * d[41] +
        2 * d[3] * d[29] * d[41] + 2 * d[2] * d[30] * d[41] - 2 * d[0] * d[32] * d[41] + 2 * d[6] * d[27] * d[42] +
        2 * d[7] * d[28] * d[42] + 2 * d[8] * d[29] * d[42] + 2 * d[0] * d[33] * d[42] + 2 * d[1] * d[34] * d[42] +
        2 * d[2] * d[35] * d[42] - 2 * d[7] * d[27] * d[43] + 2 * d[6] * d[28] * d[43] + 2 * d[1] * d[33] * d[43] -
        2 * d[0] * d[34] * d[43] - 2 * d[8] * d[27] * d[44] + 2 * d[6] * d[29] * d[44] + 2 * d[2] * d[33] * d[44] -
        2 * d[0] * d[35] * d[44];
    coeffs[97] =
        6 * d[9] * d[27] * d[36] + 2 * d[10] * d[28] * d[36] + 2 * d[11] * d[29] * d[36] + 2 * d[12] * d[30] * d[36] -
        2 * d[13] * d[31] * d[36] - 2 * d[14] * d[32] * d[36] + 2 * d[15] * d[33] * d[36] - 2 * d[16] * d[34] * d[36] -
        2 * d[17] * d[35] * d[36] + 2 * d[10] * d[27] * d[37] + 2 * d[9] * d[28] * d[37] + 2 * d[13] * d[30] * d[37] +
        2 * d[12] * d[31] * d[37] + 2 * d[16] * d[33] * d[37] + 2 * d[15] * d[34] * d[37] + 2 * d[11] * d[27] * d[38] +
        2 * d[9] * d[29] * d[38] + 2 * d[14] * d[30] * d[38] + 2 * d[12] * d[32] * d[38] + 2 * d[17] * d[33] * d[38] +
        2 * d[15] * d[35] * d[38] + 2 * d[12] * d[27] * d[39] + 2 * d[13] * d[28] * d[39] + 2 * d[14] * d[29] * d[39] +
        2 * d[9] * d[30] * d[39] + 2 * d[10] * d[31] * d[39] + 2 * d[11] * d[32] * d[39] - 2 * d[13] * d[27] * d[40] +
        2 * d[12] * d[28] * d[40] + 2 * d[10] * d[30] * d[40] - 2 * d[9] * d[31] * d[40] - 2 * d[14] * d[27] * d[41] +
        2 * d[12] * d[29] * d[41] + 2 * d[11] * d[30] * d[41] - 2 * d[9] * d[32] * d[41] + 2 * d[15] * d[27] * d[42] +
        2 * d[16] * d[28] * d[42] + 2 * d[17] * d[29] * d[42] + 2 * d[9] * d[33] * d[42] + 2 * d[10] * d[34] * d[42] +
        2 * d[11] * d[35] * d[42] - 2 * d[16] * d[27] * d[43] + 2 * d[15] * d[28] * d[43] + 2 * d[10] * d[33] * d[43] -
        2 * d[9] * d[34] * d[43] - 2 * d[17] * d[27] * d[44] + 2 * d[15] * d[29] * d[44] + 2 * d[11] * d[33] * d[44] -
        2 * d[9] * d[35] * d[44];
    coeffs[98] =
        6 * d[18] * d[27] * d[36] + 2 * d[19] * d[28] * d[36] + 2 * d[20] * d[29] * d[36] + 2 * d[21] * d[30] * d[36] -
        2 * d[22] * d[31] * d[36] - 2 * d[23] * d[32] * d[36] + 2 * d[24] * d[33] * d[36] - 2 * d[25] * d[34] * d[36] -
        2 * d[26] * d[35] * d[36] + 2 * d[19] * d[27] * d[37] + 2 * d[18] * d[28] * d[37] + 2 * d[22] * d[30] * d[37] +
        2 * d[21] * d[31] * d[37] + 2 * d[25] * d[33] * d[37] + 2 * d[24] * d[34] * d[37] + 2 * d[20] * d[27] * d[38] +
        2 * d[18] * d[29] * d[38] + 2 * d[23] * d[30] * d[38] + 2 * d[21] * d[32] * d[38] + 2 * d[26] * d[33] * d[38] +
        2 * d[24] * d[35] * d[38] + 2 * d[21] * d[27] * d[39] + 2 * d[22] * d[28] * d[39] + 2 * d[23] * d[29] * d[39] +
        2 * d[18] * d[30] * d[39] + 2 * d[19] * d[31] * d[39] + 2 * d[20] * d[32] * d[39] - 2 * d[22] * d[27] * d[40] +
        2 * d[21] * d[28] * d[40] + 2 * d[19] * d[30] * d[40] - 2 * d[18] * d[31] * d[40] - 2 * d[23] * d[27] * d[41] +
        2 * d[21] * d[29] * d[41] + 2 * d[20] * d[30] * d[41] - 2 * d[18] * d[32] * d[41] + 2 * d[24] * d[27] * d[42] +
        2 * d[25] * d[28] * d[42] + 2 * d[26] * d[29] * d[42] + 2 * d[18] * d[33] * d[42] + 2 * d[19] * d[34] * d[42] +
        2 * d[20] * d[35] * d[42] - 2 * d[25] * d[27] * d[43] + 2 * d[24] * d[28] * d[43] + 2 * d[19] * d[33] * d[43] -
        2 * d[18] * d[34] * d[43] - 2 * d[26] * d[27] * d[44] + 2 * d[24] * d[29] * d[44] + 2 * d[20] * d[33] * d[44] -
        2 * d[18] * d[35] * d[44];
    coeffs[99] =
        3 * std::pow(d[27], 2) * d[36] + std::pow(d[28], 2) * d[36] + std::pow(d[29], 2) * d[36] +
        std::pow(d[30], 2) * d[36] - std::pow(d[31], 2) * d[36] - std::pow(d[32], 2) * d[36] +
        std::pow(d[33], 2) * d[36] - std::pow(d[34], 2) * d[36] - std::pow(d[35], 2) * d[36] +
        2 * d[27] * d[28] * d[37] + 2 * d[30] * d[31] * d[37] + 2 * d[33] * d[34] * d[37] + 2 * d[27] * d[29] * d[38] +
        2 * d[30] * d[32] * d[38] + 2 * d[33] * d[35] * d[38] + 2 * d[27] * d[30] * d[39] + 2 * d[28] * d[31] * d[39] +
        2 * d[29] * d[32] * d[39] + 2 * d[28] * d[30] * d[40] - 2 * d[27] * d[31] * d[40] + 2 * d[29] * d[30] * d[41] -
        2 * d[27] * d[32] * d[41] + 2 * d[27] * d[33] * d[42] + 2 * d[28] * d[34] * d[42] + 2 * d[29] * d[35] * d[42] +
        2 * d[28] * d[33] * d[43] - 2 * d[27] * d[34] * d[43] + 2 * d[29] * d[33] * d[44] - 2 * d[27] * d[35] * d[44];
    coeffs[100] =
        3 * d[0] * std::pow(d[36], 2) + 2 * d[1] * d[36] * d[37] + d[0] * std::pow(d[37], 2) +
        2 * d[2] * d[36] * d[38] + d[0] * std::pow(d[38], 2) + 2 * d[3] * d[36] * d[39] + 2 * d[4] * d[37] * d[39] +
        2 * d[5] * d[38] * d[39] + d[0] * std::pow(d[39], 2) - 2 * d[4] * d[36] * d[40] + 2 * d[3] * d[37] * d[40] +
        2 * d[1] * d[39] * d[40] - d[0] * std::pow(d[40], 2) - 2 * d[5] * d[36] * d[41] + 2 * d[3] * d[38] * d[41] +
        2 * d[2] * d[39] * d[41] - d[0] * std::pow(d[41], 2) + 2 * d[6] * d[36] * d[42] + 2 * d[7] * d[37] * d[42] +
        2 * d[8] * d[38] * d[42] + d[0] * std::pow(d[42], 2) - 2 * d[7] * d[36] * d[43] + 2 * d[6] * d[37] * d[43] +
        2 * d[1] * d[42] * d[43] - d[0] * std::pow(d[43], 2) - 2 * d[8] * d[36] * d[44] + 2 * d[6] * d[38] * d[44] +
        2 * d[2] * d[42] * d[44] - d[0] * std::pow(d[44], 2);
    coeffs[101] =
        3 * d[9] * std::pow(d[36], 2) + 2 * d[10] * d[36] * d[37] + d[9] * std::pow(d[37], 2) +
        2 * d[11] * d[36] * d[38] + d[9] * std::pow(d[38], 2) + 2 * d[12] * d[36] * d[39] + 2 * d[13] * d[37] * d[39] +
        2 * d[14] * d[38] * d[39] + d[9] * std::pow(d[39], 2) - 2 * d[13] * d[36] * d[40] + 2 * d[12] * d[37] * d[40] +
        2 * d[10] * d[39] * d[40] - d[9] * std::pow(d[40], 2) - 2 * d[14] * d[36] * d[41] + 2 * d[12] * d[38] * d[41] +
        2 * d[11] * d[39] * d[41] - d[9] * std::pow(d[41], 2) + 2 * d[15] * d[36] * d[42] + 2 * d[16] * d[37] * d[42] +
        2 * d[17] * d[38] * d[42] + d[9] * std::pow(d[42], 2) - 2 * d[16] * d[36] * d[43] + 2 * d[15] * d[37] * d[43] +
        2 * d[10] * d[42] * d[43] - d[9] * std::pow(d[43], 2) - 2 * d[17] * d[36] * d[44] + 2 * d[15] * d[38] * d[44] +
        2 * d[11] * d[42] * d[44] - d[9] * std::pow(d[44], 2);
    coeffs[102] =
        3 * d[18] * std::pow(d[36], 2) + 2 * d[19] * d[36] * d[37] + d[18] * std::pow(d[37], 2) +
        2 * d[20] * d[36] * d[38] + d[18] * std::pow(d[38], 2) + 2 * d[21] * d[36] * d[39] + 2 * d[22] * d[37] * d[39] +
        2 * d[23] * d[38] * d[39] + d[18] * std::pow(d[39], 2) - 2 * d[22] * d[36] * d[40] + 2 * d[21] * d[37] * d[40] +
        2 * d[19] * d[39] * d[40] - d[18] * std::pow(d[40], 2) - 2 * d[23] * d[36] * d[41] + 2 * d[21] * d[38] * d[41] +
        2 * d[20] * d[39] * d[41] - d[18] * std::pow(d[41], 2) + 2 * d[24] * d[36] * d[42] + 2 * d[25] * d[37] * d[42] +
        2 * d[26] * d[38] * d[42] + d[18] * std::pow(d[42], 2) - 2 * d[25] * d[36] * d[43] + 2 * d[24] * d[37] * d[43] +
        2 * d[19] * d[42] * d[43] - d[18] * std::pow(d[43], 2) - 2 * d[26] * d[36] * d[44] + 2 * d[24] * d[38] * d[44] +
        2 * d[20] * d[42] * d[44] - d[18] * std::pow(d[44], 2);
    coeffs[103] =
        3 * d[27] * std::pow(d[36], 2) + 2 * d[28] * d[36] * d[37] + d[27] * std::pow(d[37], 2) +
        2 * d[29] * d[36] * d[38] + d[27] * std::pow(d[38], 2) + 2 * d[30] * d[36] * d[39] + 2 * d[31] * d[37] * d[39] +
        2 * d[32] * d[38] * d[39] + d[27] * std::pow(d[39], 2) - 2 * d[31] * d[36] * d[40] + 2 * d[30] * d[37] * d[40] +
        2 * d[28] * d[39] * d[40] - d[27] * std::pow(d[40], 2) - 2 * d[32] * d[36] * d[41] + 2 * d[30] * d[38] * d[41] +
        2 * d[29] * d[39] * d[41] - d[27] * std::pow(d[41], 2) + 2 * d[33] * d[36] * d[42] + 2 * d[34] * d[37] * d[42] +
        2 * d[35] * d[38] * d[42] + d[27] * std::pow(d[42], 2) - 2 * d[34] * d[36] * d[43] + 2 * d[33] * d[37] * d[43] +
        2 * d[28] * d[42] * d[43] - d[27] * std::pow(d[43], 2) - 2 * d[35] * d[36] * d[44] + 2 * d[33] * d[38] * d[44] +
        2 * d[29] * d[42] * d[44] - d[27] * std::pow(d[44], 2);
    coeffs[104] = std::pow(d[36], 3) + d[36] * std::pow(d[37], 2) + d[36] * std::pow(d[38], 2) +
                  d[36] * std::pow(d[39], 2) + 2 * d[37] * d[39] * d[40] - d[36] * std::pow(d[40], 2) +
                  2 * d[38] * d[39] * d[41] - d[36] * std::pow(d[41], 2) + d[36] * std::pow(d[42], 2) +
                  2 * d[37] * d[42] * d[43] - d[36] * std::pow(d[43], 2) + 2 * d[38] * d[42] * d[44] -
                  d[36] * std::pow(d[44], 2);
    coeffs[105] = std::pow(d[0], 2) * d[1] + std::pow(d[1], 3) + d[1] * std::pow(d[2], 2) - d[1] * std::pow(d[3], 2) +
                  2 * d[0] * d[3] * d[4] + d[1] * std::pow(d[4], 2) + 2 * d[2] * d[4] * d[5] -
                  d[1] * std::pow(d[5], 2) - d[1] * std::pow(d[6], 2) + 2 * d[0] * d[6] * d[7] +
                  d[1] * std::pow(d[7], 2) + 2 * d[2] * d[7] * d[8] - d[1] * std::pow(d[8], 2);
    coeffs[106] = 2 * d[0] * d[1] * d[9] + 2 * d[3] * d[4] * d[9] + 2 * d[6] * d[7] * d[9] + std::pow(d[0], 2) * d[10] +
                  3 * std::pow(d[1], 2) * d[10] + std::pow(d[2], 2) * d[10] - std::pow(d[3], 2) * d[10] +
                  std::pow(d[4], 2) * d[10] - std::pow(d[5], 2) * d[10] - std::pow(d[6], 2) * d[10] +
                  std::pow(d[7], 2) * d[10] - std::pow(d[8], 2) * d[10] + 2 * d[1] * d[2] * d[11] +
                  2 * d[4] * d[5] * d[11] + 2 * d[7] * d[8] * d[11] - 2 * d[1] * d[3] * d[12] +
                  2 * d[0] * d[4] * d[12] + 2 * d[0] * d[3] * d[13] + 2 * d[1] * d[4] * d[13] +
                  2 * d[2] * d[5] * d[13] + 2 * d[2] * d[4] * d[14] - 2 * d[1] * d[5] * d[14] -
                  2 * d[1] * d[6] * d[15] + 2 * d[0] * d[7] * d[15] + 2 * d[0] * d[6] * d[16] +
                  2 * d[1] * d[7] * d[16] + 2 * d[2] * d[8] * d[16] + 2 * d[2] * d[7] * d[17] - 2 * d[1] * d[8] * d[17];
    coeffs[107] =
        d[1] * std::pow(d[9], 2) + 2 * d[0] * d[9] * d[10] + 3 * d[1] * std::pow(d[10], 2) + 2 * d[2] * d[10] * d[11] +
        d[1] * std::pow(d[11], 2) + 2 * d[4] * d[9] * d[12] - 2 * d[3] * d[10] * d[12] - d[1] * std::pow(d[12], 2) +
        2 * d[3] * d[9] * d[13] + 2 * d[4] * d[10] * d[13] + 2 * d[5] * d[11] * d[13] + 2 * d[0] * d[12] * d[13] +
        d[1] * std::pow(d[13], 2) - 2 * d[5] * d[10] * d[14] + 2 * d[4] * d[11] * d[14] + 2 * d[2] * d[13] * d[14] -
        d[1] * std::pow(d[14], 2) + 2 * d[7] * d[9] * d[15] - 2 * d[6] * d[10] * d[15] - d[1] * std::pow(d[15], 2) +
        2 * d[6] * d[9] * d[16] + 2 * d[7] * d[10] * d[16] + 2 * d[8] * d[11] * d[16] + 2 * d[0] * d[15] * d[16] +
        d[1] * std::pow(d[16], 2) - 2 * d[8] * d[10] * d[17] + 2 * d[7] * d[11] * d[17] + 2 * d[2] * d[16] * d[17] -
        d[1] * std::pow(d[17], 2);
    coeffs[108] = std::pow(d[9], 2) * d[10] + std::pow(d[10], 3) + d[10] * std::pow(d[11], 2) -
                  d[10] * std::pow(d[12], 2) + 2 * d[9] * d[12] * d[13] + d[10] * std::pow(d[13], 2) +
                  2 * d[11] * d[13] * d[14] - d[10] * std::pow(d[14], 2) - d[10] * std::pow(d[15], 2) +
                  2 * d[9] * d[15] * d[16] + d[10] * std::pow(d[16], 2) + 2 * d[11] * d[16] * d[17] -
                  d[10] * std::pow(d[17], 2);
    coeffs[109] =
        2 * d[0] * d[1] * d[18] + 2 * d[3] * d[4] * d[18] + 2 * d[6] * d[7] * d[18] + std::pow(d[0], 2) * d[19] +
        3 * std::pow(d[1], 2) * d[19] + std::pow(d[2], 2) * d[19] - std::pow(d[3], 2) * d[19] +
        std::pow(d[4], 2) * d[19] - std::pow(d[5], 2) * d[19] - std::pow(d[6], 2) * d[19] + std::pow(d[7], 2) * d[19] -
        std::pow(d[8], 2) * d[19] + 2 * d[1] * d[2] * d[20] + 2 * d[4] * d[5] * d[20] + 2 * d[7] * d[8] * d[20] -
        2 * d[1] * d[3] * d[21] + 2 * d[0] * d[4] * d[21] + 2 * d[0] * d[3] * d[22] + 2 * d[1] * d[4] * d[22] +
        2 * d[2] * d[5] * d[22] + 2 * d[2] * d[4] * d[23] - 2 * d[1] * d[5] * d[23] - 2 * d[1] * d[6] * d[24] +
        2 * d[0] * d[7] * d[24] + 2 * d[0] * d[6] * d[25] + 2 * d[1] * d[7] * d[25] + 2 * d[2] * d[8] * d[25] +
        2 * d[2] * d[7] * d[26] - 2 * d[1] * d[8] * d[26];
    coeffs[110] =
        2 * d[1] * d[9] * d[18] + 2 * d[0] * d[10] * d[18] + 2 * d[4] * d[12] * d[18] + 2 * d[3] * d[13] * d[18] +
        2 * d[7] * d[15] * d[18] + 2 * d[6] * d[16] * d[18] + 2 * d[0] * d[9] * d[19] + 6 * d[1] * d[10] * d[19] +
        2 * d[2] * d[11] * d[19] - 2 * d[3] * d[12] * d[19] + 2 * d[4] * d[13] * d[19] - 2 * d[5] * d[14] * d[19] -
        2 * d[6] * d[15] * d[19] + 2 * d[7] * d[16] * d[19] - 2 * d[8] * d[17] * d[19] + 2 * d[2] * d[10] * d[20] +
        2 * d[1] * d[11] * d[20] + 2 * d[5] * d[13] * d[20] + 2 * d[4] * d[14] * d[20] + 2 * d[8] * d[16] * d[20] +
        2 * d[7] * d[17] * d[20] + 2 * d[4] * d[9] * d[21] - 2 * d[3] * d[10] * d[21] - 2 * d[1] * d[12] * d[21] +
        2 * d[0] * d[13] * d[21] + 2 * d[3] * d[9] * d[22] + 2 * d[4] * d[10] * d[22] + 2 * d[5] * d[11] * d[22] +
        2 * d[0] * d[12] * d[22] + 2 * d[1] * d[13] * d[22] + 2 * d[2] * d[14] * d[22] - 2 * d[5] * d[10] * d[23] +
        2 * d[4] * d[11] * d[23] + 2 * d[2] * d[13] * d[23] - 2 * d[1] * d[14] * d[23] + 2 * d[7] * d[9] * d[24] -
        2 * d[6] * d[10] * d[24] - 2 * d[1] * d[15] * d[24] + 2 * d[0] * d[16] * d[24] + 2 * d[6] * d[9] * d[25] +
        2 * d[7] * d[10] * d[25] + 2 * d[8] * d[11] * d[25] + 2 * d[0] * d[15] * d[25] + 2 * d[1] * d[16] * d[25] +
        2 * d[2] * d[17] * d[25] - 2 * d[8] * d[10] * d[26] + 2 * d[7] * d[11] * d[26] + 2 * d[2] * d[16] * d[26] -
        2 * d[1] * d[17] * d[26];
    coeffs[111] =
        2 * d[9] * d[10] * d[18] + 2 * d[12] * d[13] * d[18] + 2 * d[15] * d[16] * d[18] + std::pow(d[9], 2) * d[19] +
        3 * std::pow(d[10], 2) * d[19] + std::pow(d[11], 2) * d[19] - std::pow(d[12], 2) * d[19] +
        std::pow(d[13], 2) * d[19] - std::pow(d[14], 2) * d[19] - std::pow(d[15], 2) * d[19] +
        std::pow(d[16], 2) * d[19] - std::pow(d[17], 2) * d[19] + 2 * d[10] * d[11] * d[20] +
        2 * d[13] * d[14] * d[20] + 2 * d[16] * d[17] * d[20] - 2 * d[10] * d[12] * d[21] + 2 * d[9] * d[13] * d[21] +
        2 * d[9] * d[12] * d[22] + 2 * d[10] * d[13] * d[22] + 2 * d[11] * d[14] * d[22] + 2 * d[11] * d[13] * d[23] -
        2 * d[10] * d[14] * d[23] - 2 * d[10] * d[15] * d[24] + 2 * d[9] * d[16] * d[24] + 2 * d[9] * d[15] * d[25] +
        2 * d[10] * d[16] * d[25] + 2 * d[11] * d[17] * d[25] + 2 * d[11] * d[16] * d[26] - 2 * d[10] * d[17] * d[26];
    coeffs[112] =
        d[1] * std::pow(d[18], 2) + 2 * d[0] * d[18] * d[19] + 3 * d[1] * std::pow(d[19], 2) +
        2 * d[2] * d[19] * d[20] + d[1] * std::pow(d[20], 2) + 2 * d[4] * d[18] * d[21] - 2 * d[3] * d[19] * d[21] -
        d[1] * std::pow(d[21], 2) + 2 * d[3] * d[18] * d[22] + 2 * d[4] * d[19] * d[22] + 2 * d[5] * d[20] * d[22] +
        2 * d[0] * d[21] * d[22] + d[1] * std::pow(d[22], 2) - 2 * d[5] * d[19] * d[23] + 2 * d[4] * d[20] * d[23] +
        2 * d[2] * d[22] * d[23] - d[1] * std::pow(d[23], 2) + 2 * d[7] * d[18] * d[24] - 2 * d[6] * d[19] * d[24] -
        d[1] * std::pow(d[24], 2) + 2 * d[6] * d[18] * d[25] + 2 * d[7] * d[19] * d[25] + 2 * d[8] * d[20] * d[25] +
        2 * d[0] * d[24] * d[25] + d[1] * std::pow(d[25], 2) - 2 * d[8] * d[19] * d[26] + 2 * d[7] * d[20] * d[26] +
        2 * d[2] * d[25] * d[26] - d[1] * std::pow(d[26], 2);
    coeffs[113] =
        d[10] * std::pow(d[18], 2) + 2 * d[9] * d[18] * d[19] + 3 * d[10] * std::pow(d[19], 2) +
        2 * d[11] * d[19] * d[20] + d[10] * std::pow(d[20], 2) + 2 * d[13] * d[18] * d[21] - 2 * d[12] * d[19] * d[21] -
        d[10] * std::pow(d[21], 2) + 2 * d[12] * d[18] * d[22] + 2 * d[13] * d[19] * d[22] + 2 * d[14] * d[20] * d[22] +
        2 * d[9] * d[21] * d[22] + d[10] * std::pow(d[22], 2) - 2 * d[14] * d[19] * d[23] + 2 * d[13] * d[20] * d[23] +
        2 * d[11] * d[22] * d[23] - d[10] * std::pow(d[23], 2) + 2 * d[16] * d[18] * d[24] - 2 * d[15] * d[19] * d[24] -
        d[10] * std::pow(d[24], 2) + 2 * d[15] * d[18] * d[25] + 2 * d[16] * d[19] * d[25] + 2 * d[17] * d[20] * d[25] +
        2 * d[9] * d[24] * d[25] + d[10] * std::pow(d[25], 2) - 2 * d[17] * d[19] * d[26] + 2 * d[16] * d[20] * d[26] +
        2 * d[11] * d[25] * d[26] - d[10] * std::pow(d[26], 2);
    coeffs[114] = std::pow(d[18], 2) * d[19] + std::pow(d[19], 3) + d[19] * std::pow(d[20], 2) -
                  d[19] * std::pow(d[21], 2) + 2 * d[18] * d[21] * d[22] + d[19] * std::pow(d[22], 2) +
                  2 * d[20] * d[22] * d[23] - d[19] * std::pow(d[23], 2) - d[19] * std::pow(d[24], 2) +
                  2 * d[18] * d[24] * d[25] + d[19] * std::pow(d[25], 2) + 2 * d[20] * d[25] * d[26] -
                  d[19] * std::pow(d[26], 2);
    coeffs[115] =
        2 * d[0] * d[1] * d[27] + 2 * d[3] * d[4] * d[27] + 2 * d[6] * d[7] * d[27] + std::pow(d[0], 2) * d[28] +
        3 * std::pow(d[1], 2) * d[28] + std::pow(d[2], 2) * d[28] - std::pow(d[3], 2) * d[28] +
        std::pow(d[4], 2) * d[28] - std::pow(d[5], 2) * d[28] - std::pow(d[6], 2) * d[28] + std::pow(d[7], 2) * d[28] -
        std::pow(d[8], 2) * d[28] + 2 * d[1] * d[2] * d[29] + 2 * d[4] * d[5] * d[29] + 2 * d[7] * d[8] * d[29] -
        2 * d[1] * d[3] * d[30] + 2 * d[0] * d[4] * d[30] + 2 * d[0] * d[3] * d[31] + 2 * d[1] * d[4] * d[31] +
        2 * d[2] * d[5] * d[31] + 2 * d[2] * d[4] * d[32] - 2 * d[1] * d[5] * d[32] - 2 * d[1] * d[6] * d[33] +
        2 * d[0] * d[7] * d[33] + 2 * d[0] * d[6] * d[34] + 2 * d[1] * d[7] * d[34] + 2 * d[2] * d[8] * d[34] +
        2 * d[2] * d[7] * d[35] - 2 * d[1] * d[8] * d[35];
    coeffs[116] =
        2 * d[1] * d[9] * d[27] + 2 * d[0] * d[10] * d[27] + 2 * d[4] * d[12] * d[27] + 2 * d[3] * d[13] * d[27] +
        2 * d[7] * d[15] * d[27] + 2 * d[6] * d[16] * d[27] + 2 * d[0] * d[9] * d[28] + 6 * d[1] * d[10] * d[28] +
        2 * d[2] * d[11] * d[28] - 2 * d[3] * d[12] * d[28] + 2 * d[4] * d[13] * d[28] - 2 * d[5] * d[14] * d[28] -
        2 * d[6] * d[15] * d[28] + 2 * d[7] * d[16] * d[28] - 2 * d[8] * d[17] * d[28] + 2 * d[2] * d[10] * d[29] +
        2 * d[1] * d[11] * d[29] + 2 * d[5] * d[13] * d[29] + 2 * d[4] * d[14] * d[29] + 2 * d[8] * d[16] * d[29] +
        2 * d[7] * d[17] * d[29] + 2 * d[4] * d[9] * d[30] - 2 * d[3] * d[10] * d[30] - 2 * d[1] * d[12] * d[30] +
        2 * d[0] * d[13] * d[30] + 2 * d[3] * d[9] * d[31] + 2 * d[4] * d[10] * d[31] + 2 * d[5] * d[11] * d[31] +
        2 * d[0] * d[12] * d[31] + 2 * d[1] * d[13] * d[31] + 2 * d[2] * d[14] * d[31] - 2 * d[5] * d[10] * d[32] +
        2 * d[4] * d[11] * d[32] + 2 * d[2] * d[13] * d[32] - 2 * d[1] * d[14] * d[32] + 2 * d[7] * d[9] * d[33] -
        2 * d[6] * d[10] * d[33] - 2 * d[1] * d[15] * d[33] + 2 * d[0] * d[16] * d[33] + 2 * d[6] * d[9] * d[34] +
        2 * d[7] * d[10] * d[34] + 2 * d[8] * d[11] * d[34] + 2 * d[0] * d[15] * d[34] + 2 * d[1] * d[16] * d[34] +
        2 * d[2] * d[17] * d[34] - 2 * d[8] * d[10] * d[35] + 2 * d[7] * d[11] * d[35] + 2 * d[2] * d[16] * d[35] -
        2 * d[1] * d[17] * d[35];
    coeffs[117] =
        2 * d[9] * d[10] * d[27] + 2 * d[12] * d[13] * d[27] + 2 * d[15] * d[16] * d[27] + std::pow(d[9], 2) * d[28] +
        3 * std::pow(d[10], 2) * d[28] + std::pow(d[11], 2) * d[28] - std::pow(d[12], 2) * d[28] +
        std::pow(d[13], 2) * d[28] - std::pow(d[14], 2) * d[28] - std::pow(d[15], 2) * d[28] +
        std::pow(d[16], 2) * d[28] - std::pow(d[17], 2) * d[28] + 2 * d[10] * d[11] * d[29] +
        2 * d[13] * d[14] * d[29] + 2 * d[16] * d[17] * d[29] - 2 * d[10] * d[12] * d[30] + 2 * d[9] * d[13] * d[30] +
        2 * d[9] * d[12] * d[31] + 2 * d[10] * d[13] * d[31] + 2 * d[11] * d[14] * d[31] + 2 * d[11] * d[13] * d[32] -
        2 * d[10] * d[14] * d[32] - 2 * d[10] * d[15] * d[33] + 2 * d[9] * d[16] * d[33] + 2 * d[9] * d[15] * d[34] +
        2 * d[10] * d[16] * d[34] + 2 * d[11] * d[17] * d[34] + 2 * d[11] * d[16] * d[35] - 2 * d[10] * d[17] * d[35];
    coeffs[118] =
        2 * d[1] * d[18] * d[27] + 2 * d[0] * d[19] * d[27] + 2 * d[4] * d[21] * d[27] + 2 * d[3] * d[22] * d[27] +
        2 * d[7] * d[24] * d[27] + 2 * d[6] * d[25] * d[27] + 2 * d[0] * d[18] * d[28] + 6 * d[1] * d[19] * d[28] +
        2 * d[2] * d[20] * d[28] - 2 * d[3] * d[21] * d[28] + 2 * d[4] * d[22] * d[28] - 2 * d[5] * d[23] * d[28] -
        2 * d[6] * d[24] * d[28] + 2 * d[7] * d[25] * d[28] - 2 * d[8] * d[26] * d[28] + 2 * d[2] * d[19] * d[29] +
        2 * d[1] * d[20] * d[29] + 2 * d[5] * d[22] * d[29] + 2 * d[4] * d[23] * d[29] + 2 * d[8] * d[25] * d[29] +
        2 * d[7] * d[26] * d[29] + 2 * d[4] * d[18] * d[30] - 2 * d[3] * d[19] * d[30] - 2 * d[1] * d[21] * d[30] +
        2 * d[0] * d[22] * d[30] + 2 * d[3] * d[18] * d[31] + 2 * d[4] * d[19] * d[31] + 2 * d[5] * d[20] * d[31] +
        2 * d[0] * d[21] * d[31] + 2 * d[1] * d[22] * d[31] + 2 * d[2] * d[23] * d[31] - 2 * d[5] * d[19] * d[32] +
        2 * d[4] * d[20] * d[32] + 2 * d[2] * d[22] * d[32] - 2 * d[1] * d[23] * d[32] + 2 * d[7] * d[18] * d[33] -
        2 * d[6] * d[19] * d[33] - 2 * d[1] * d[24] * d[33] + 2 * d[0] * d[25] * d[33] + 2 * d[6] * d[18] * d[34] +
        2 * d[7] * d[19] * d[34] + 2 * d[8] * d[20] * d[34] + 2 * d[0] * d[24] * d[34] + 2 * d[1] * d[25] * d[34] +
        2 * d[2] * d[26] * d[34] - 2 * d[8] * d[19] * d[35] + 2 * d[7] * d[20] * d[35] + 2 * d[2] * d[25] * d[35] -
        2 * d[1] * d[26] * d[35];
    coeffs[119] =
        2 * d[10] * d[18] * d[27] + 2 * d[9] * d[19] * d[27] + 2 * d[13] * d[21] * d[27] + 2 * d[12] * d[22] * d[27] +
        2 * d[16] * d[24] * d[27] + 2 * d[15] * d[25] * d[27] + 2 * d[9] * d[18] * d[28] + 6 * d[10] * d[19] * d[28] +
        2 * d[11] * d[20] * d[28] - 2 * d[12] * d[21] * d[28] + 2 * d[13] * d[22] * d[28] - 2 * d[14] * d[23] * d[28] -
        2 * d[15] * d[24] * d[28] + 2 * d[16] * d[25] * d[28] - 2 * d[17] * d[26] * d[28] + 2 * d[11] * d[19] * d[29] +
        2 * d[10] * d[20] * d[29] + 2 * d[14] * d[22] * d[29] + 2 * d[13] * d[23] * d[29] + 2 * d[17] * d[25] * d[29] +
        2 * d[16] * d[26] * d[29] + 2 * d[13] * d[18] * d[30] - 2 * d[12] * d[19] * d[30] - 2 * d[10] * d[21] * d[30] +
        2 * d[9] * d[22] * d[30] + 2 * d[12] * d[18] * d[31] + 2 * d[13] * d[19] * d[31] + 2 * d[14] * d[20] * d[31] +
        2 * d[9] * d[21] * d[31] + 2 * d[10] * d[22] * d[31] + 2 * d[11] * d[23] * d[31] - 2 * d[14] * d[19] * d[32] +
        2 * d[13] * d[20] * d[32] + 2 * d[11] * d[22] * d[32] - 2 * d[10] * d[23] * d[32] + 2 * d[16] * d[18] * d[33] -
        2 * d[15] * d[19] * d[33] - 2 * d[10] * d[24] * d[33] + 2 * d[9] * d[25] * d[33] + 2 * d[15] * d[18] * d[34] +
        2 * d[16] * d[19] * d[34] + 2 * d[17] * d[20] * d[34] + 2 * d[9] * d[24] * d[34] + 2 * d[10] * d[25] * d[34] +
        2 * d[11] * d[26] * d[34] - 2 * d[17] * d[19] * d[35] + 2 * d[16] * d[20] * d[35] + 2 * d[11] * d[25] * d[35] -
        2 * d[10] * d[26] * d[35];
    coeffs[120] =
        2 * d[18] * d[19] * d[27] + 2 * d[21] * d[22] * d[27] + 2 * d[24] * d[25] * d[27] + std::pow(d[18], 2) * d[28] +
        3 * std::pow(d[19], 2) * d[28] + std::pow(d[20], 2) * d[28] - std::pow(d[21], 2) * d[28] +
        std::pow(d[22], 2) * d[28] - std::pow(d[23], 2) * d[28] - std::pow(d[24], 2) * d[28] +
        std::pow(d[25], 2) * d[28] - std::pow(d[26], 2) * d[28] + 2 * d[19] * d[20] * d[29] +
        2 * d[22] * d[23] * d[29] + 2 * d[25] * d[26] * d[29] - 2 * d[19] * d[21] * d[30] + 2 * d[18] * d[22] * d[30] +
        2 * d[18] * d[21] * d[31] + 2 * d[19] * d[22] * d[31] + 2 * d[20] * d[23] * d[31] + 2 * d[20] * d[22] * d[32] -
        2 * d[19] * d[23] * d[32] - 2 * d[19] * d[24] * d[33] + 2 * d[18] * d[25] * d[33] + 2 * d[18] * d[24] * d[34] +
        2 * d[19] * d[25] * d[34] + 2 * d[20] * d[26] * d[34] + 2 * d[20] * d[25] * d[35] - 2 * d[19] * d[26] * d[35];
    coeffs[121] =
        d[1] * std::pow(d[27], 2) + 2 * d[0] * d[27] * d[28] + 3 * d[1] * std::pow(d[28], 2) +
        2 * d[2] * d[28] * d[29] + d[1] * std::pow(d[29], 2) + 2 * d[4] * d[27] * d[30] - 2 * d[3] * d[28] * d[30] -
        d[1] * std::pow(d[30], 2) + 2 * d[3] * d[27] * d[31] + 2 * d[4] * d[28] * d[31] + 2 * d[5] * d[29] * d[31] +
        2 * d[0] * d[30] * d[31] + d[1] * std::pow(d[31], 2) - 2 * d[5] * d[28] * d[32] + 2 * d[4] * d[29] * d[32] +
        2 * d[2] * d[31] * d[32] - d[1] * std::pow(d[32], 2) + 2 * d[7] * d[27] * d[33] - 2 * d[6] * d[28] * d[33] -
        d[1] * std::pow(d[33], 2) + 2 * d[6] * d[27] * d[34] + 2 * d[7] * d[28] * d[34] + 2 * d[8] * d[29] * d[34] +
        2 * d[0] * d[33] * d[34] + d[1] * std::pow(d[34], 2) - 2 * d[8] * d[28] * d[35] + 2 * d[7] * d[29] * d[35] +
        2 * d[2] * d[34] * d[35] - d[1] * std::pow(d[35], 2);
    coeffs[122] =
        d[10] * std::pow(d[27], 2) + 2 * d[9] * d[27] * d[28] + 3 * d[10] * std::pow(d[28], 2) +
        2 * d[11] * d[28] * d[29] + d[10] * std::pow(d[29], 2) + 2 * d[13] * d[27] * d[30] - 2 * d[12] * d[28] * d[30] -
        d[10] * std::pow(d[30], 2) + 2 * d[12] * d[27] * d[31] + 2 * d[13] * d[28] * d[31] + 2 * d[14] * d[29] * d[31] +
        2 * d[9] * d[30] * d[31] + d[10] * std::pow(d[31], 2) - 2 * d[14] * d[28] * d[32] + 2 * d[13] * d[29] * d[32] +
        2 * d[11] * d[31] * d[32] - d[10] * std::pow(d[32], 2) + 2 * d[16] * d[27] * d[33] - 2 * d[15] * d[28] * d[33] -
        d[10] * std::pow(d[33], 2) + 2 * d[15] * d[27] * d[34] + 2 * d[16] * d[28] * d[34] + 2 * d[17] * d[29] * d[34] +
        2 * d[9] * d[33] * d[34] + d[10] * std::pow(d[34], 2) - 2 * d[17] * d[28] * d[35] + 2 * d[16] * d[29] * d[35] +
        2 * d[11] * d[34] * d[35] - d[10] * std::pow(d[35], 2);
    coeffs[123] =
        d[19] * std::pow(d[27], 2) + 2 * d[18] * d[27] * d[28] + 3 * d[19] * std::pow(d[28], 2) +
        2 * d[20] * d[28] * d[29] + d[19] * std::pow(d[29], 2) + 2 * d[22] * d[27] * d[30] - 2 * d[21] * d[28] * d[30] -
        d[19] * std::pow(d[30], 2) + 2 * d[21] * d[27] * d[31] + 2 * d[22] * d[28] * d[31] + 2 * d[23] * d[29] * d[31] +
        2 * d[18] * d[30] * d[31] + d[19] * std::pow(d[31], 2) - 2 * d[23] * d[28] * d[32] + 2 * d[22] * d[29] * d[32] +
        2 * d[20] * d[31] * d[32] - d[19] * std::pow(d[32], 2) + 2 * d[25] * d[27] * d[33] - 2 * d[24] * d[28] * d[33] -
        d[19] * std::pow(d[33], 2) + 2 * d[24] * d[27] * d[34] + 2 * d[25] * d[28] * d[34] + 2 * d[26] * d[29] * d[34] +
        2 * d[18] * d[33] * d[34] + d[19] * std::pow(d[34], 2) - 2 * d[26] * d[28] * d[35] + 2 * d[25] * d[29] * d[35] +
        2 * d[20] * d[34] * d[35] - d[19] * std::pow(d[35], 2);
    coeffs[124] = std::pow(d[27], 2) * d[28] + std::pow(d[28], 3) + d[28] * std::pow(d[29], 2) -
                  d[28] * std::pow(d[30], 2) + 2 * d[27] * d[30] * d[31] + d[28] * std::pow(d[31], 2) +
                  2 * d[29] * d[31] * d[32] - d[28] * std::pow(d[32], 2) - d[28] * std::pow(d[33], 2) +
                  2 * d[27] * d[33] * d[34] + d[28] * std::pow(d[34], 2) + 2 * d[29] * d[34] * d[35] -
                  d[28] * std::pow(d[35], 2);
    coeffs[125] =
        2 * d[0] * d[1] * d[36] + 2 * d[3] * d[4] * d[36] + 2 * d[6] * d[7] * d[36] + std::pow(d[0], 2) * d[37] +
        3 * std::pow(d[1], 2) * d[37] + std::pow(d[2], 2) * d[37] - std::pow(d[3], 2) * d[37] +
        std::pow(d[4], 2) * d[37] - std::pow(d[5], 2) * d[37] - std::pow(d[6], 2) * d[37] + std::pow(d[7], 2) * d[37] -
        std::pow(d[8], 2) * d[37] + 2 * d[1] * d[2] * d[38] + 2 * d[4] * d[5] * d[38] + 2 * d[7] * d[8] * d[38] -
        2 * d[1] * d[3] * d[39] + 2 * d[0] * d[4] * d[39] + 2 * d[0] * d[3] * d[40] + 2 * d[1] * d[4] * d[40] +
        2 * d[2] * d[5] * d[40] + 2 * d[2] * d[4] * d[41] - 2 * d[1] * d[5] * d[41] - 2 * d[1] * d[6] * d[42] +
        2 * d[0] * d[7] * d[42] + 2 * d[0] * d[6] * d[43] + 2 * d[1] * d[7] * d[43] + 2 * d[2] * d[8] * d[43] +
        2 * d[2] * d[7] * d[44] - 2 * d[1] * d[8] * d[44];
    coeffs[126] =
        2 * d[1] * d[9] * d[36] + 2 * d[0] * d[10] * d[36] + 2 * d[4] * d[12] * d[36] + 2 * d[3] * d[13] * d[36] +
        2 * d[7] * d[15] * d[36] + 2 * d[6] * d[16] * d[36] + 2 * d[0] * d[9] * d[37] + 6 * d[1] * d[10] * d[37] +
        2 * d[2] * d[11] * d[37] - 2 * d[3] * d[12] * d[37] + 2 * d[4] * d[13] * d[37] - 2 * d[5] * d[14] * d[37] -
        2 * d[6] * d[15] * d[37] + 2 * d[7] * d[16] * d[37] - 2 * d[8] * d[17] * d[37] + 2 * d[2] * d[10] * d[38] +
        2 * d[1] * d[11] * d[38] + 2 * d[5] * d[13] * d[38] + 2 * d[4] * d[14] * d[38] + 2 * d[8] * d[16] * d[38] +
        2 * d[7] * d[17] * d[38] + 2 * d[4] * d[9] * d[39] - 2 * d[3] * d[10] * d[39] - 2 * d[1] * d[12] * d[39] +
        2 * d[0] * d[13] * d[39] + 2 * d[3] * d[9] * d[40] + 2 * d[4] * d[10] * d[40] + 2 * d[5] * d[11] * d[40] +
        2 * d[0] * d[12] * d[40] + 2 * d[1] * d[13] * d[40] + 2 * d[2] * d[14] * d[40] - 2 * d[5] * d[10] * d[41] +
        2 * d[4] * d[11] * d[41] + 2 * d[2] * d[13] * d[41] - 2 * d[1] * d[14] * d[41] + 2 * d[7] * d[9] * d[42] -
        2 * d[6] * d[10] * d[42] - 2 * d[1] * d[15] * d[42] + 2 * d[0] * d[16] * d[42] + 2 * d[6] * d[9] * d[43] +
        2 * d[7] * d[10] * d[43] + 2 * d[8] * d[11] * d[43] + 2 * d[0] * d[15] * d[43] + 2 * d[1] * d[16] * d[43] +
        2 * d[2] * d[17] * d[43] - 2 * d[8] * d[10] * d[44] + 2 * d[7] * d[11] * d[44] + 2 * d[2] * d[16] * d[44] -
        2 * d[1] * d[17] * d[44];
    coeffs[127] =
        2 * d[9] * d[10] * d[36] + 2 * d[12] * d[13] * d[36] + 2 * d[15] * d[16] * d[36] + std::pow(d[9], 2) * d[37] +
        3 * std::pow(d[10], 2) * d[37] + std::pow(d[11], 2) * d[37] - std::pow(d[12], 2) * d[37] +
        std::pow(d[13], 2) * d[37] - std::pow(d[14], 2) * d[37] - std::pow(d[15], 2) * d[37] +
        std::pow(d[16], 2) * d[37] - std::pow(d[17], 2) * d[37] + 2 * d[10] * d[11] * d[38] +
        2 * d[13] * d[14] * d[38] + 2 * d[16] * d[17] * d[38] - 2 * d[10] * d[12] * d[39] + 2 * d[9] * d[13] * d[39] +
        2 * d[9] * d[12] * d[40] + 2 * d[10] * d[13] * d[40] + 2 * d[11] * d[14] * d[40] + 2 * d[11] * d[13] * d[41] -
        2 * d[10] * d[14] * d[41] - 2 * d[10] * d[15] * d[42] + 2 * d[9] * d[16] * d[42] + 2 * d[9] * d[15] * d[43] +
        2 * d[10] * d[16] * d[43] + 2 * d[11] * d[17] * d[43] + 2 * d[11] * d[16] * d[44] - 2 * d[10] * d[17] * d[44];
    coeffs[128] =
        2 * d[1] * d[18] * d[36] + 2 * d[0] * d[19] * d[36] + 2 * d[4] * d[21] * d[36] + 2 * d[3] * d[22] * d[36] +
        2 * d[7] * d[24] * d[36] + 2 * d[6] * d[25] * d[36] + 2 * d[0] * d[18] * d[37] + 6 * d[1] * d[19] * d[37] +
        2 * d[2] * d[20] * d[37] - 2 * d[3] * d[21] * d[37] + 2 * d[4] * d[22] * d[37] - 2 * d[5] * d[23] * d[37] -
        2 * d[6] * d[24] * d[37] + 2 * d[7] * d[25] * d[37] - 2 * d[8] * d[26] * d[37] + 2 * d[2] * d[19] * d[38] +
        2 * d[1] * d[20] * d[38] + 2 * d[5] * d[22] * d[38] + 2 * d[4] * d[23] * d[38] + 2 * d[8] * d[25] * d[38] +
        2 * d[7] * d[26] * d[38] + 2 * d[4] * d[18] * d[39] - 2 * d[3] * d[19] * d[39] - 2 * d[1] * d[21] * d[39] +
        2 * d[0] * d[22] * d[39] + 2 * d[3] * d[18] * d[40] + 2 * d[4] * d[19] * d[40] + 2 * d[5] * d[20] * d[40] +
        2 * d[0] * d[21] * d[40] + 2 * d[1] * d[22] * d[40] + 2 * d[2] * d[23] * d[40] - 2 * d[5] * d[19] * d[41] +
        2 * d[4] * d[20] * d[41] + 2 * d[2] * d[22] * d[41] - 2 * d[1] * d[23] * d[41] + 2 * d[7] * d[18] * d[42] -
        2 * d[6] * d[19] * d[42] - 2 * d[1] * d[24] * d[42] + 2 * d[0] * d[25] * d[42] + 2 * d[6] * d[18] * d[43] +
        2 * d[7] * d[19] * d[43] + 2 * d[8] * d[20] * d[43] + 2 * d[0] * d[24] * d[43] + 2 * d[1] * d[25] * d[43] +
        2 * d[2] * d[26] * d[43] - 2 * d[8] * d[19] * d[44] + 2 * d[7] * d[20] * d[44] + 2 * d[2] * d[25] * d[44] -
        2 * d[1] * d[26] * d[44];
    coeffs[129] =
        2 * d[10] * d[18] * d[36] + 2 * d[9] * d[19] * d[36] + 2 * d[13] * d[21] * d[36] + 2 * d[12] * d[22] * d[36] +
        2 * d[16] * d[24] * d[36] + 2 * d[15] * d[25] * d[36] + 2 * d[9] * d[18] * d[37] + 6 * d[10] * d[19] * d[37] +
        2 * d[11] * d[20] * d[37] - 2 * d[12] * d[21] * d[37] + 2 * d[13] * d[22] * d[37] - 2 * d[14] * d[23] * d[37] -
        2 * d[15] * d[24] * d[37] + 2 * d[16] * d[25] * d[37] - 2 * d[17] * d[26] * d[37] + 2 * d[11] * d[19] * d[38] +
        2 * d[10] * d[20] * d[38] + 2 * d[14] * d[22] * d[38] + 2 * d[13] * d[23] * d[38] + 2 * d[17] * d[25] * d[38] +
        2 * d[16] * d[26] * d[38] + 2 * d[13] * d[18] * d[39] - 2 * d[12] * d[19] * d[39] - 2 * d[10] * d[21] * d[39] +
        2 * d[9] * d[22] * d[39] + 2 * d[12] * d[18] * d[40] + 2 * d[13] * d[19] * d[40] + 2 * d[14] * d[20] * d[40] +
        2 * d[9] * d[21] * d[40] + 2 * d[10] * d[22] * d[40] + 2 * d[11] * d[23] * d[40] - 2 * d[14] * d[19] * d[41] +
        2 * d[13] * d[20] * d[41] + 2 * d[11] * d[22] * d[41] - 2 * d[10] * d[23] * d[41] + 2 * d[16] * d[18] * d[42] -
        2 * d[15] * d[19] * d[42] - 2 * d[10] * d[24] * d[42] + 2 * d[9] * d[25] * d[42] + 2 * d[15] * d[18] * d[43] +
        2 * d[16] * d[19] * d[43] + 2 * d[17] * d[20] * d[43] + 2 * d[9] * d[24] * d[43] + 2 * d[10] * d[25] * d[43] +
        2 * d[11] * d[26] * d[43] - 2 * d[17] * d[19] * d[44] + 2 * d[16] * d[20] * d[44] + 2 * d[11] * d[25] * d[44] -
        2 * d[10] * d[26] * d[44];
    coeffs[130] =
        2 * d[18] * d[19] * d[36] + 2 * d[21] * d[22] * d[36] + 2 * d[24] * d[25] * d[36] + std::pow(d[18], 2) * d[37] +
        3 * std::pow(d[19], 2) * d[37] + std::pow(d[20], 2) * d[37] - std::pow(d[21], 2) * d[37] +
        std::pow(d[22], 2) * d[37] - std::pow(d[23], 2) * d[37] - std::pow(d[24], 2) * d[37] +
        std::pow(d[25], 2) * d[37] - std::pow(d[26], 2) * d[37] + 2 * d[19] * d[20] * d[38] +
        2 * d[22] * d[23] * d[38] + 2 * d[25] * d[26] * d[38] - 2 * d[19] * d[21] * d[39] + 2 * d[18] * d[22] * d[39] +
        2 * d[18] * d[21] * d[40] + 2 * d[19] * d[22] * d[40] + 2 * d[20] * d[23] * d[40] + 2 * d[20] * d[22] * d[41] -
        2 * d[19] * d[23] * d[41] - 2 * d[19] * d[24] * d[42] + 2 * d[18] * d[25] * d[42] + 2 * d[18] * d[24] * d[43] +
        2 * d[19] * d[25] * d[43] + 2 * d[20] * d[26] * d[43] + 2 * d[20] * d[25] * d[44] - 2 * d[19] * d[26] * d[44];
    coeffs[131] =
        2 * d[1] * d[27] * d[36] + 2 * d[0] * d[28] * d[36] + 2 * d[4] * d[30] * d[36] + 2 * d[3] * d[31] * d[36] +
        2 * d[7] * d[33] * d[36] + 2 * d[6] * d[34] * d[36] + 2 * d[0] * d[27] * d[37] + 6 * d[1] * d[28] * d[37] +
        2 * d[2] * d[29] * d[37] - 2 * d[3] * d[30] * d[37] + 2 * d[4] * d[31] * d[37] - 2 * d[5] * d[32] * d[37] -
        2 * d[6] * d[33] * d[37] + 2 * d[7] * d[34] * d[37] - 2 * d[8] * d[35] * d[37] + 2 * d[2] * d[28] * d[38] +
        2 * d[1] * d[29] * d[38] + 2 * d[5] * d[31] * d[38] + 2 * d[4] * d[32] * d[38] + 2 * d[8] * d[34] * d[38] +
        2 * d[7] * d[35] * d[38] + 2 * d[4] * d[27] * d[39] - 2 * d[3] * d[28] * d[39] - 2 * d[1] * d[30] * d[39] +
        2 * d[0] * d[31] * d[39] + 2 * d[3] * d[27] * d[40] + 2 * d[4] * d[28] * d[40] + 2 * d[5] * d[29] * d[40] +
        2 * d[0] * d[30] * d[40] + 2 * d[1] * d[31] * d[40] + 2 * d[2] * d[32] * d[40] - 2 * d[5] * d[28] * d[41] +
        2 * d[4] * d[29] * d[41] + 2 * d[2] * d[31] * d[41] - 2 * d[1] * d[32] * d[41] + 2 * d[7] * d[27] * d[42] -
        2 * d[6] * d[28] * d[42] - 2 * d[1] * d[33] * d[42] + 2 * d[0] * d[34] * d[42] + 2 * d[6] * d[27] * d[43] +
        2 * d[7] * d[28] * d[43] + 2 * d[8] * d[29] * d[43] + 2 * d[0] * d[33] * d[43] + 2 * d[1] * d[34] * d[43] +
        2 * d[2] * d[35] * d[43] - 2 * d[8] * d[28] * d[44] + 2 * d[7] * d[29] * d[44] + 2 * d[2] * d[34] * d[44] -
        2 * d[1] * d[35] * d[44];
    coeffs[132] =
        2 * d[10] * d[27] * d[36] + 2 * d[9] * d[28] * d[36] + 2 * d[13] * d[30] * d[36] + 2 * d[12] * d[31] * d[36] +
        2 * d[16] * d[33] * d[36] + 2 * d[15] * d[34] * d[36] + 2 * d[9] * d[27] * d[37] + 6 * d[10] * d[28] * d[37] +
        2 * d[11] * d[29] * d[37] - 2 * d[12] * d[30] * d[37] + 2 * d[13] * d[31] * d[37] - 2 * d[14] * d[32] * d[37] -
        2 * d[15] * d[33] * d[37] + 2 * d[16] * d[34] * d[37] - 2 * d[17] * d[35] * d[37] + 2 * d[11] * d[28] * d[38] +
        2 * d[10] * d[29] * d[38] + 2 * d[14] * d[31] * d[38] + 2 * d[13] * d[32] * d[38] + 2 * d[17] * d[34] * d[38] +
        2 * d[16] * d[35] * d[38] + 2 * d[13] * d[27] * d[39] - 2 * d[12] * d[28] * d[39] - 2 * d[10] * d[30] * d[39] +
        2 * d[9] * d[31] * d[39] + 2 * d[12] * d[27] * d[40] + 2 * d[13] * d[28] * d[40] + 2 * d[14] * d[29] * d[40] +
        2 * d[9] * d[30] * d[40] + 2 * d[10] * d[31] * d[40] + 2 * d[11] * d[32] * d[40] - 2 * d[14] * d[28] * d[41] +
        2 * d[13] * d[29] * d[41] + 2 * d[11] * d[31] * d[41] - 2 * d[10] * d[32] * d[41] + 2 * d[16] * d[27] * d[42] -
        2 * d[15] * d[28] * d[42] - 2 * d[10] * d[33] * d[42] + 2 * d[9] * d[34] * d[42] + 2 * d[15] * d[27] * d[43] +
        2 * d[16] * d[28] * d[43] + 2 * d[17] * d[29] * d[43] + 2 * d[9] * d[33] * d[43] + 2 * d[10] * d[34] * d[43] +
        2 * d[11] * d[35] * d[43] - 2 * d[17] * d[28] * d[44] + 2 * d[16] * d[29] * d[44] + 2 * d[11] * d[34] * d[44] -
        2 * d[10] * d[35] * d[44];
    coeffs[133] =
        2 * d[19] * d[27] * d[36] + 2 * d[18] * d[28] * d[36] + 2 * d[22] * d[30] * d[36] + 2 * d[21] * d[31] * d[36] +
        2 * d[25] * d[33] * d[36] + 2 * d[24] * d[34] * d[36] + 2 * d[18] * d[27] * d[37] + 6 * d[19] * d[28] * d[37] +
        2 * d[20] * d[29] * d[37] - 2 * d[21] * d[30] * d[37] + 2 * d[22] * d[31] * d[37] - 2 * d[23] * d[32] * d[37] -
        2 * d[24] * d[33] * d[37] + 2 * d[25] * d[34] * d[37] - 2 * d[26] * d[35] * d[37] + 2 * d[20] * d[28] * d[38] +
        2 * d[19] * d[29] * d[38] + 2 * d[23] * d[31] * d[38] + 2 * d[22] * d[32] * d[38] + 2 * d[26] * d[34] * d[38] +
        2 * d[25] * d[35] * d[38] + 2 * d[22] * d[27] * d[39] - 2 * d[21] * d[28] * d[39] - 2 * d[19] * d[30] * d[39] +
        2 * d[18] * d[31] * d[39] + 2 * d[21] * d[27] * d[40] + 2 * d[22] * d[28] * d[40] + 2 * d[23] * d[29] * d[40] +
        2 * d[18] * d[30] * d[40] + 2 * d[19] * d[31] * d[40] + 2 * d[20] * d[32] * d[40] - 2 * d[23] * d[28] * d[41] +
        2 * d[22] * d[29] * d[41] + 2 * d[20] * d[31] * d[41] - 2 * d[19] * d[32] * d[41] + 2 * d[25] * d[27] * d[42] -
        2 * d[24] * d[28] * d[42] - 2 * d[19] * d[33] * d[42] + 2 * d[18] * d[34] * d[42] + 2 * d[24] * d[27] * d[43] +
        2 * d[25] * d[28] * d[43] + 2 * d[26] * d[29] * d[43] + 2 * d[18] * d[33] * d[43] + 2 * d[19] * d[34] * d[43] +
        2 * d[20] * d[35] * d[43] - 2 * d[26] * d[28] * d[44] + 2 * d[25] * d[29] * d[44] + 2 * d[20] * d[34] * d[44] -
        2 * d[19] * d[35] * d[44];
    coeffs[134] =
        2 * d[27] * d[28] * d[36] + 2 * d[30] * d[31] * d[36] + 2 * d[33] * d[34] * d[36] + std::pow(d[27], 2) * d[37] +
        3 * std::pow(d[28], 2) * d[37] + std::pow(d[29], 2) * d[37] - std::pow(d[30], 2) * d[37] +
        std::pow(d[31], 2) * d[37] - std::pow(d[32], 2) * d[37] - std::pow(d[33], 2) * d[37] +
        std::pow(d[34], 2) * d[37] - std::pow(d[35], 2) * d[37] + 2 * d[28] * d[29] * d[38] +
        2 * d[31] * d[32] * d[38] + 2 * d[34] * d[35] * d[38] - 2 * d[28] * d[30] * d[39] + 2 * d[27] * d[31] * d[39] +
        2 * d[27] * d[30] * d[40] + 2 * d[28] * d[31] * d[40] + 2 * d[29] * d[32] * d[40] + 2 * d[29] * d[31] * d[41] -
        2 * d[28] * d[32] * d[41] - 2 * d[28] * d[33] * d[42] + 2 * d[27] * d[34] * d[42] + 2 * d[27] * d[33] * d[43] +
        2 * d[28] * d[34] * d[43] + 2 * d[29] * d[35] * d[43] + 2 * d[29] * d[34] * d[44] - 2 * d[28] * d[35] * d[44];
    coeffs[135] =
        d[1] * std::pow(d[36], 2) + 2 * d[0] * d[36] * d[37] + 3 * d[1] * std::pow(d[37], 2) +
        2 * d[2] * d[37] * d[38] + d[1] * std::pow(d[38], 2) + 2 * d[4] * d[36] * d[39] - 2 * d[3] * d[37] * d[39] -
        d[1] * std::pow(d[39], 2) + 2 * d[3] * d[36] * d[40] + 2 * d[4] * d[37] * d[40] + 2 * d[5] * d[38] * d[40] +
        2 * d[0] * d[39] * d[40] + d[1] * std::pow(d[40], 2) - 2 * d[5] * d[37] * d[41] + 2 * d[4] * d[38] * d[41] +
        2 * d[2] * d[40] * d[41] - d[1] * std::pow(d[41], 2) + 2 * d[7] * d[36] * d[42] - 2 * d[6] * d[37] * d[42] -
        d[1] * std::pow(d[42], 2) + 2 * d[6] * d[36] * d[43] + 2 * d[7] * d[37] * d[43] + 2 * d[8] * d[38] * d[43] +
        2 * d[0] * d[42] * d[43] + d[1] * std::pow(d[43], 2) - 2 * d[8] * d[37] * d[44] + 2 * d[7] * d[38] * d[44] +
        2 * d[2] * d[43] * d[44] - d[1] * std::pow(d[44], 2);
    coeffs[136] =
        d[10] * std::pow(d[36], 2) + 2 * d[9] * d[36] * d[37] + 3 * d[10] * std::pow(d[37], 2) +
        2 * d[11] * d[37] * d[38] + d[10] * std::pow(d[38], 2) + 2 * d[13] * d[36] * d[39] - 2 * d[12] * d[37] * d[39] -
        d[10] * std::pow(d[39], 2) + 2 * d[12] * d[36] * d[40] + 2 * d[13] * d[37] * d[40] + 2 * d[14] * d[38] * d[40] +
        2 * d[9] * d[39] * d[40] + d[10] * std::pow(d[40], 2) - 2 * d[14] * d[37] * d[41] + 2 * d[13] * d[38] * d[41] +
        2 * d[11] * d[40] * d[41] - d[10] * std::pow(d[41], 2) + 2 * d[16] * d[36] * d[42] - 2 * d[15] * d[37] * d[42] -
        d[10] * std::pow(d[42], 2) + 2 * d[15] * d[36] * d[43] + 2 * d[16] * d[37] * d[43] + 2 * d[17] * d[38] * d[43] +
        2 * d[9] * d[42] * d[43] + d[10] * std::pow(d[43], 2) - 2 * d[17] * d[37] * d[44] + 2 * d[16] * d[38] * d[44] +
        2 * d[11] * d[43] * d[44] - d[10] * std::pow(d[44], 2);
    coeffs[137] =
        d[19] * std::pow(d[36], 2) + 2 * d[18] * d[36] * d[37] + 3 * d[19] * std::pow(d[37], 2) +
        2 * d[20] * d[37] * d[38] + d[19] * std::pow(d[38], 2) + 2 * d[22] * d[36] * d[39] - 2 * d[21] * d[37] * d[39] -
        d[19] * std::pow(d[39], 2) + 2 * d[21] * d[36] * d[40] + 2 * d[22] * d[37] * d[40] + 2 * d[23] * d[38] * d[40] +
        2 * d[18] * d[39] * d[40] + d[19] * std::pow(d[40], 2) - 2 * d[23] * d[37] * d[41] + 2 * d[22] * d[38] * d[41] +
        2 * d[20] * d[40] * d[41] - d[19] * std::pow(d[41], 2) + 2 * d[25] * d[36] * d[42] - 2 * d[24] * d[37] * d[42] -
        d[19] * std::pow(d[42], 2) + 2 * d[24] * d[36] * d[43] + 2 * d[25] * d[37] * d[43] + 2 * d[26] * d[38] * d[43] +
        2 * d[18] * d[42] * d[43] + d[19] * std::pow(d[43], 2) - 2 * d[26] * d[37] * d[44] + 2 * d[25] * d[38] * d[44] +
        2 * d[20] * d[43] * d[44] - d[19] * std::pow(d[44], 2);
    coeffs[138] =
        d[28] * std::pow(d[36], 2) + 2 * d[27] * d[36] * d[37] + 3 * d[28] * std::pow(d[37], 2) +
        2 * d[29] * d[37] * d[38] + d[28] * std::pow(d[38], 2) + 2 * d[31] * d[36] * d[39] - 2 * d[30] * d[37] * d[39] -
        d[28] * std::pow(d[39], 2) + 2 * d[30] * d[36] * d[40] + 2 * d[31] * d[37] * d[40] + 2 * d[32] * d[38] * d[40] +
        2 * d[27] * d[39] * d[40] + d[28] * std::pow(d[40], 2) - 2 * d[32] * d[37] * d[41] + 2 * d[31] * d[38] * d[41] +
        2 * d[29] * d[40] * d[41] - d[28] * std::pow(d[41], 2) + 2 * d[34] * d[36] * d[42] - 2 * d[33] * d[37] * d[42] -
        d[28] * std::pow(d[42], 2) + 2 * d[33] * d[36] * d[43] + 2 * d[34] * d[37] * d[43] + 2 * d[35] * d[38] * d[43] +
        2 * d[27] * d[42] * d[43] + d[28] * std::pow(d[43], 2) - 2 * d[35] * d[37] * d[44] + 2 * d[34] * d[38] * d[44] +
        2 * d[29] * d[43] * d[44] - d[28] * std::pow(d[44], 2);
    coeffs[139] = std::pow(d[36], 2) * d[37] + std::pow(d[37], 3) + d[37] * std::pow(d[38], 2) -
                  d[37] * std::pow(d[39], 2) + 2 * d[36] * d[39] * d[40] + d[37] * std::pow(d[40], 2) +
                  2 * d[38] * d[40] * d[41] - d[37] * std::pow(d[41], 2) - d[37] * std::pow(d[42], 2) +
                  2 * d[36] * d[42] * d[43] + d[37] * std::pow(d[43], 2) + 2 * d[38] * d[43] * d[44] -
                  d[37] * std::pow(d[44], 2);
    coeffs[140] = std::pow(d[0], 2) * d[2] + std::pow(d[1], 2) * d[2] + std::pow(d[2], 3) - d[2] * std::pow(d[3], 2) -
                  d[2] * std::pow(d[4], 2) + 2 * d[0] * d[3] * d[5] + 2 * d[1] * d[4] * d[5] +
                  d[2] * std::pow(d[5], 2) - d[2] * std::pow(d[6], 2) - d[2] * std::pow(d[7], 2) +
                  2 * d[0] * d[6] * d[8] + 2 * d[1] * d[7] * d[8] + d[2] * std::pow(d[8], 2);
    coeffs[141] = 2 * d[0] * d[2] * d[9] + 2 * d[3] * d[5] * d[9] + 2 * d[6] * d[8] * d[9] + 2 * d[1] * d[2] * d[10] +
                  2 * d[4] * d[5] * d[10] + 2 * d[7] * d[8] * d[10] + std::pow(d[0], 2) * d[11] +
                  std::pow(d[1], 2) * d[11] + 3 * std::pow(d[2], 2) * d[11] - std::pow(d[3], 2) * d[11] -
                  std::pow(d[4], 2) * d[11] + std::pow(d[5], 2) * d[11] - std::pow(d[6], 2) * d[11] -
                  std::pow(d[7], 2) * d[11] + std::pow(d[8], 2) * d[11] - 2 * d[2] * d[3] * d[12] +
                  2 * d[0] * d[5] * d[12] - 2 * d[2] * d[4] * d[13] + 2 * d[1] * d[5] * d[13] +
                  2 * d[0] * d[3] * d[14] + 2 * d[1] * d[4] * d[14] + 2 * d[2] * d[5] * d[14] -
                  2 * d[2] * d[6] * d[15] + 2 * d[0] * d[8] * d[15] - 2 * d[2] * d[7] * d[16] +
                  2 * d[1] * d[8] * d[16] + 2 * d[0] * d[6] * d[17] + 2 * d[1] * d[7] * d[17] + 2 * d[2] * d[8] * d[17];
    coeffs[142] =
        d[2] * std::pow(d[9], 2) + d[2] * std::pow(d[10], 2) + 2 * d[0] * d[9] * d[11] + 2 * d[1] * d[10] * d[11] +
        3 * d[2] * std::pow(d[11], 2) + 2 * d[5] * d[9] * d[12] - 2 * d[3] * d[11] * d[12] - d[2] * std::pow(d[12], 2) +
        2 * d[5] * d[10] * d[13] - 2 * d[4] * d[11] * d[13] - d[2] * std::pow(d[13], 2) + 2 * d[3] * d[9] * d[14] +
        2 * d[4] * d[10] * d[14] + 2 * d[5] * d[11] * d[14] + 2 * d[0] * d[12] * d[14] + 2 * d[1] * d[13] * d[14] +
        d[2] * std::pow(d[14], 2) + 2 * d[8] * d[9] * d[15] - 2 * d[6] * d[11] * d[15] - d[2] * std::pow(d[15], 2) +
        2 * d[8] * d[10] * d[16] - 2 * d[7] * d[11] * d[16] - d[2] * std::pow(d[16], 2) + 2 * d[6] * d[9] * d[17] +
        2 * d[7] * d[10] * d[17] + 2 * d[8] * d[11] * d[17] + 2 * d[0] * d[15] * d[17] + 2 * d[1] * d[16] * d[17] +
        d[2] * std::pow(d[17], 2);
    coeffs[143] = std::pow(d[9], 2) * d[11] + std::pow(d[10], 2) * d[11] + std::pow(d[11], 3) -
                  d[11] * std::pow(d[12], 2) - d[11] * std::pow(d[13], 2) + 2 * d[9] * d[12] * d[14] +
                  2 * d[10] * d[13] * d[14] + d[11] * std::pow(d[14], 2) - d[11] * std::pow(d[15], 2) -
                  d[11] * std::pow(d[16], 2) + 2 * d[9] * d[15] * d[17] + 2 * d[10] * d[16] * d[17] +
                  d[11] * std::pow(d[17], 2);
    coeffs[144] =
        2 * d[0] * d[2] * d[18] + 2 * d[3] * d[5] * d[18] + 2 * d[6] * d[8] * d[18] + 2 * d[1] * d[2] * d[19] +
        2 * d[4] * d[5] * d[19] + 2 * d[7] * d[8] * d[19] + std::pow(d[0], 2) * d[20] + std::pow(d[1], 2) * d[20] +
        3 * std::pow(d[2], 2) * d[20] - std::pow(d[3], 2) * d[20] - std::pow(d[4], 2) * d[20] +
        std::pow(d[5], 2) * d[20] - std::pow(d[6], 2) * d[20] - std::pow(d[7], 2) * d[20] + std::pow(d[8], 2) * d[20] -
        2 * d[2] * d[3] * d[21] + 2 * d[0] * d[5] * d[21] - 2 * d[2] * d[4] * d[22] + 2 * d[1] * d[5] * d[22] +
        2 * d[0] * d[3] * d[23] + 2 * d[1] * d[4] * d[23] + 2 * d[2] * d[5] * d[23] - 2 * d[2] * d[6] * d[24] +
        2 * d[0] * d[8] * d[24] - 2 * d[2] * d[7] * d[25] + 2 * d[1] * d[8] * d[25] + 2 * d[0] * d[6] * d[26] +
        2 * d[1] * d[7] * d[26] + 2 * d[2] * d[8] * d[26];
    coeffs[145] =
        2 * d[2] * d[9] * d[18] + 2 * d[0] * d[11] * d[18] + 2 * d[5] * d[12] * d[18] + 2 * d[3] * d[14] * d[18] +
        2 * d[8] * d[15] * d[18] + 2 * d[6] * d[17] * d[18] + 2 * d[2] * d[10] * d[19] + 2 * d[1] * d[11] * d[19] +
        2 * d[5] * d[13] * d[19] + 2 * d[4] * d[14] * d[19] + 2 * d[8] * d[16] * d[19] + 2 * d[7] * d[17] * d[19] +
        2 * d[0] * d[9] * d[20] + 2 * d[1] * d[10] * d[20] + 6 * d[2] * d[11] * d[20] - 2 * d[3] * d[12] * d[20] -
        2 * d[4] * d[13] * d[20] + 2 * d[5] * d[14] * d[20] - 2 * d[6] * d[15] * d[20] - 2 * d[7] * d[16] * d[20] +
        2 * d[8] * d[17] * d[20] + 2 * d[5] * d[9] * d[21] - 2 * d[3] * d[11] * d[21] - 2 * d[2] * d[12] * d[21] +
        2 * d[0] * d[14] * d[21] + 2 * d[5] * d[10] * d[22] - 2 * d[4] * d[11] * d[22] - 2 * d[2] * d[13] * d[22] +
        2 * d[1] * d[14] * d[22] + 2 * d[3] * d[9] * d[23] + 2 * d[4] * d[10] * d[23] + 2 * d[5] * d[11] * d[23] +
        2 * d[0] * d[12] * d[23] + 2 * d[1] * d[13] * d[23] + 2 * d[2] * d[14] * d[23] + 2 * d[8] * d[9] * d[24] -
        2 * d[6] * d[11] * d[24] - 2 * d[2] * d[15] * d[24] + 2 * d[0] * d[17] * d[24] + 2 * d[8] * d[10] * d[25] -
        2 * d[7] * d[11] * d[25] - 2 * d[2] * d[16] * d[25] + 2 * d[1] * d[17] * d[25] + 2 * d[6] * d[9] * d[26] +
        2 * d[7] * d[10] * d[26] + 2 * d[8] * d[11] * d[26] + 2 * d[0] * d[15] * d[26] + 2 * d[1] * d[16] * d[26] +
        2 * d[2] * d[17] * d[26];
    coeffs[146] =
        2 * d[9] * d[11] * d[18] + 2 * d[12] * d[14] * d[18] + 2 * d[15] * d[17] * d[18] + 2 * d[10] * d[11] * d[19] +
        2 * d[13] * d[14] * d[19] + 2 * d[16] * d[17] * d[19] + std::pow(d[9], 2) * d[20] + std::pow(d[10], 2) * d[20] +
        3 * std::pow(d[11], 2) * d[20] - std::pow(d[12], 2) * d[20] - std::pow(d[13], 2) * d[20] +
        std::pow(d[14], 2) * d[20] - std::pow(d[15], 2) * d[20] - std::pow(d[16], 2) * d[20] +
        std::pow(d[17], 2) * d[20] - 2 * d[11] * d[12] * d[21] + 2 * d[9] * d[14] * d[21] - 2 * d[11] * d[13] * d[22] +
        2 * d[10] * d[14] * d[22] + 2 * d[9] * d[12] * d[23] + 2 * d[10] * d[13] * d[23] + 2 * d[11] * d[14] * d[23] -
        2 * d[11] * d[15] * d[24] + 2 * d[9] * d[17] * d[24] - 2 * d[11] * d[16] * d[25] + 2 * d[10] * d[17] * d[25] +
        2 * d[9] * d[15] * d[26] + 2 * d[10] * d[16] * d[26] + 2 * d[11] * d[17] * d[26];
    coeffs[147] =
        d[2] * std::pow(d[18], 2) + d[2] * std::pow(d[19], 2) + 2 * d[0] * d[18] * d[20] + 2 * d[1] * d[19] * d[20] +
        3 * d[2] * std::pow(d[20], 2) + 2 * d[5] * d[18] * d[21] - 2 * d[3] * d[20] * d[21] -
        d[2] * std::pow(d[21], 2) + 2 * d[5] * d[19] * d[22] - 2 * d[4] * d[20] * d[22] - d[2] * std::pow(d[22], 2) +
        2 * d[3] * d[18] * d[23] + 2 * d[4] * d[19] * d[23] + 2 * d[5] * d[20] * d[23] + 2 * d[0] * d[21] * d[23] +
        2 * d[1] * d[22] * d[23] + d[2] * std::pow(d[23], 2) + 2 * d[8] * d[18] * d[24] - 2 * d[6] * d[20] * d[24] -
        d[2] * std::pow(d[24], 2) + 2 * d[8] * d[19] * d[25] - 2 * d[7] * d[20] * d[25] - d[2] * std::pow(d[25], 2) +
        2 * d[6] * d[18] * d[26] + 2 * d[7] * d[19] * d[26] + 2 * d[8] * d[20] * d[26] + 2 * d[0] * d[24] * d[26] +
        2 * d[1] * d[25] * d[26] + d[2] * std::pow(d[26], 2);
    coeffs[148] =
        d[11] * std::pow(d[18], 2) + d[11] * std::pow(d[19], 2) + 2 * d[9] * d[18] * d[20] + 2 * d[10] * d[19] * d[20] +
        3 * d[11] * std::pow(d[20], 2) + 2 * d[14] * d[18] * d[21] - 2 * d[12] * d[20] * d[21] -
        d[11] * std::pow(d[21], 2) + 2 * d[14] * d[19] * d[22] - 2 * d[13] * d[20] * d[22] -
        d[11] * std::pow(d[22], 2) + 2 * d[12] * d[18] * d[23] + 2 * d[13] * d[19] * d[23] + 2 * d[14] * d[20] * d[23] +
        2 * d[9] * d[21] * d[23] + 2 * d[10] * d[22] * d[23] + d[11] * std::pow(d[23], 2) + 2 * d[17] * d[18] * d[24] -
        2 * d[15] * d[20] * d[24] - d[11] * std::pow(d[24], 2) + 2 * d[17] * d[19] * d[25] - 2 * d[16] * d[20] * d[25] -
        d[11] * std::pow(d[25], 2) + 2 * d[15] * d[18] * d[26] + 2 * d[16] * d[19] * d[26] + 2 * d[17] * d[20] * d[26] +
        2 * d[9] * d[24] * d[26] + 2 * d[10] * d[25] * d[26] + d[11] * std::pow(d[26], 2);
    coeffs[149] = std::pow(d[18], 2) * d[20] + std::pow(d[19], 2) * d[20] + std::pow(d[20], 3) -
                  d[20] * std::pow(d[21], 2) - d[20] * std::pow(d[22], 2) + 2 * d[18] * d[21] * d[23] +
                  2 * d[19] * d[22] * d[23] + d[20] * std::pow(d[23], 2) - d[20] * std::pow(d[24], 2) -
                  d[20] * std::pow(d[25], 2) + 2 * d[18] * d[24] * d[26] + 2 * d[19] * d[25] * d[26] +
                  d[20] * std::pow(d[26], 2);
    coeffs[150] =
        2 * d[0] * d[2] * d[27] + 2 * d[3] * d[5] * d[27] + 2 * d[6] * d[8] * d[27] + 2 * d[1] * d[2] * d[28] +
        2 * d[4] * d[5] * d[28] + 2 * d[7] * d[8] * d[28] + std::pow(d[0], 2) * d[29] + std::pow(d[1], 2) * d[29] +
        3 * std::pow(d[2], 2) * d[29] - std::pow(d[3], 2) * d[29] - std::pow(d[4], 2) * d[29] +
        std::pow(d[5], 2) * d[29] - std::pow(d[6], 2) * d[29] - std::pow(d[7], 2) * d[29] + std::pow(d[8], 2) * d[29] -
        2 * d[2] * d[3] * d[30] + 2 * d[0] * d[5] * d[30] - 2 * d[2] * d[4] * d[31] + 2 * d[1] * d[5] * d[31] +
        2 * d[0] * d[3] * d[32] + 2 * d[1] * d[4] * d[32] + 2 * d[2] * d[5] * d[32] - 2 * d[2] * d[6] * d[33] +
        2 * d[0] * d[8] * d[33] - 2 * d[2] * d[7] * d[34] + 2 * d[1] * d[8] * d[34] + 2 * d[0] * d[6] * d[35] +
        2 * d[1] * d[7] * d[35] + 2 * d[2] * d[8] * d[35];
    coeffs[151] =
        2 * d[2] * d[9] * d[27] + 2 * d[0] * d[11] * d[27] + 2 * d[5] * d[12] * d[27] + 2 * d[3] * d[14] * d[27] +
        2 * d[8] * d[15] * d[27] + 2 * d[6] * d[17] * d[27] + 2 * d[2] * d[10] * d[28] + 2 * d[1] * d[11] * d[28] +
        2 * d[5] * d[13] * d[28] + 2 * d[4] * d[14] * d[28] + 2 * d[8] * d[16] * d[28] + 2 * d[7] * d[17] * d[28] +
        2 * d[0] * d[9] * d[29] + 2 * d[1] * d[10] * d[29] + 6 * d[2] * d[11] * d[29] - 2 * d[3] * d[12] * d[29] -
        2 * d[4] * d[13] * d[29] + 2 * d[5] * d[14] * d[29] - 2 * d[6] * d[15] * d[29] - 2 * d[7] * d[16] * d[29] +
        2 * d[8] * d[17] * d[29] + 2 * d[5] * d[9] * d[30] - 2 * d[3] * d[11] * d[30] - 2 * d[2] * d[12] * d[30] +
        2 * d[0] * d[14] * d[30] + 2 * d[5] * d[10] * d[31] - 2 * d[4] * d[11] * d[31] - 2 * d[2] * d[13] * d[31] +
        2 * d[1] * d[14] * d[31] + 2 * d[3] * d[9] * d[32] + 2 * d[4] * d[10] * d[32] + 2 * d[5] * d[11] * d[32] +
        2 * d[0] * d[12] * d[32] + 2 * d[1] * d[13] * d[32] + 2 * d[2] * d[14] * d[32] + 2 * d[8] * d[9] * d[33] -
        2 * d[6] * d[11] * d[33] - 2 * d[2] * d[15] * d[33] + 2 * d[0] * d[17] * d[33] + 2 * d[8] * d[10] * d[34] -
        2 * d[7] * d[11] * d[34] - 2 * d[2] * d[16] * d[34] + 2 * d[1] * d[17] * d[34] + 2 * d[6] * d[9] * d[35] +
        2 * d[7] * d[10] * d[35] + 2 * d[8] * d[11] * d[35] + 2 * d[0] * d[15] * d[35] + 2 * d[1] * d[16] * d[35] +
        2 * d[2] * d[17] * d[35];
    coeffs[152] =
        2 * d[9] * d[11] * d[27] + 2 * d[12] * d[14] * d[27] + 2 * d[15] * d[17] * d[27] + 2 * d[10] * d[11] * d[28] +
        2 * d[13] * d[14] * d[28] + 2 * d[16] * d[17] * d[28] + std::pow(d[9], 2) * d[29] + std::pow(d[10], 2) * d[29] +
        3 * std::pow(d[11], 2) * d[29] - std::pow(d[12], 2) * d[29] - std::pow(d[13], 2) * d[29] +
        std::pow(d[14], 2) * d[29] - std::pow(d[15], 2) * d[29] - std::pow(d[16], 2) * d[29] +
        std::pow(d[17], 2) * d[29] - 2 * d[11] * d[12] * d[30] + 2 * d[9] * d[14] * d[30] - 2 * d[11] * d[13] * d[31] +
        2 * d[10] * d[14] * d[31] + 2 * d[9] * d[12] * d[32] + 2 * d[10] * d[13] * d[32] + 2 * d[11] * d[14] * d[32] -
        2 * d[11] * d[15] * d[33] + 2 * d[9] * d[17] * d[33] - 2 * d[11] * d[16] * d[34] + 2 * d[10] * d[17] * d[34] +
        2 * d[9] * d[15] * d[35] + 2 * d[10] * d[16] * d[35] + 2 * d[11] * d[17] * d[35];
    coeffs[153] =
        2 * d[2] * d[18] * d[27] + 2 * d[0] * d[20] * d[27] + 2 * d[5] * d[21] * d[27] + 2 * d[3] * d[23] * d[27] +
        2 * d[8] * d[24] * d[27] + 2 * d[6] * d[26] * d[27] + 2 * d[2] * d[19] * d[28] + 2 * d[1] * d[20] * d[28] +
        2 * d[5] * d[22] * d[28] + 2 * d[4] * d[23] * d[28] + 2 * d[8] * d[25] * d[28] + 2 * d[7] * d[26] * d[28] +
        2 * d[0] * d[18] * d[29] + 2 * d[1] * d[19] * d[29] + 6 * d[2] * d[20] * d[29] - 2 * d[3] * d[21] * d[29] -
        2 * d[4] * d[22] * d[29] + 2 * d[5] * d[23] * d[29] - 2 * d[6] * d[24] * d[29] - 2 * d[7] * d[25] * d[29] +
        2 * d[8] * d[26] * d[29] + 2 * d[5] * d[18] * d[30] - 2 * d[3] * d[20] * d[30] - 2 * d[2] * d[21] * d[30] +
        2 * d[0] * d[23] * d[30] + 2 * d[5] * d[19] * d[31] - 2 * d[4] * d[20] * d[31] - 2 * d[2] * d[22] * d[31] +
        2 * d[1] * d[23] * d[31] + 2 * d[3] * d[18] * d[32] + 2 * d[4] * d[19] * d[32] + 2 * d[5] * d[20] * d[32] +
        2 * d[0] * d[21] * d[32] + 2 * d[1] * d[22] * d[32] + 2 * d[2] * d[23] * d[32] + 2 * d[8] * d[18] * d[33] -
        2 * d[6] * d[20] * d[33] - 2 * d[2] * d[24] * d[33] + 2 * d[0] * d[26] * d[33] + 2 * d[8] * d[19] * d[34] -
        2 * d[7] * d[20] * d[34] - 2 * d[2] * d[25] * d[34] + 2 * d[1] * d[26] * d[34] + 2 * d[6] * d[18] * d[35] +
        2 * d[7] * d[19] * d[35] + 2 * d[8] * d[20] * d[35] + 2 * d[0] * d[24] * d[35] + 2 * d[1] * d[25] * d[35] +
        2 * d[2] * d[26] * d[35];
    coeffs[154] =
        2 * d[11] * d[18] * d[27] + 2 * d[9] * d[20] * d[27] + 2 * d[14] * d[21] * d[27] + 2 * d[12] * d[23] * d[27] +
        2 * d[17] * d[24] * d[27] + 2 * d[15] * d[26] * d[27] + 2 * d[11] * d[19] * d[28] + 2 * d[10] * d[20] * d[28] +
        2 * d[14] * d[22] * d[28] + 2 * d[13] * d[23] * d[28] + 2 * d[17] * d[25] * d[28] + 2 * d[16] * d[26] * d[28] +
        2 * d[9] * d[18] * d[29] + 2 * d[10] * d[19] * d[29] + 6 * d[11] * d[20] * d[29] - 2 * d[12] * d[21] * d[29] -
        2 * d[13] * d[22] * d[29] + 2 * d[14] * d[23] * d[29] - 2 * d[15] * d[24] * d[29] - 2 * d[16] * d[25] * d[29] +
        2 * d[17] * d[26] * d[29] + 2 * d[14] * d[18] * d[30] - 2 * d[12] * d[20] * d[30] - 2 * d[11] * d[21] * d[30] +
        2 * d[9] * d[23] * d[30] + 2 * d[14] * d[19] * d[31] - 2 * d[13] * d[20] * d[31] - 2 * d[11] * d[22] * d[31] +
        2 * d[10] * d[23] * d[31] + 2 * d[12] * d[18] * d[32] + 2 * d[13] * d[19] * d[32] + 2 * d[14] * d[20] * d[32] +
        2 * d[9] * d[21] * d[32] + 2 * d[10] * d[22] * d[32] + 2 * d[11] * d[23] * d[32] + 2 * d[17] * d[18] * d[33] -
        2 * d[15] * d[20] * d[33] - 2 * d[11] * d[24] * d[33] + 2 * d[9] * d[26] * d[33] + 2 * d[17] * d[19] * d[34] -
        2 * d[16] * d[20] * d[34] - 2 * d[11] * d[25] * d[34] + 2 * d[10] * d[26] * d[34] + 2 * d[15] * d[18] * d[35] +
        2 * d[16] * d[19] * d[35] + 2 * d[17] * d[20] * d[35] + 2 * d[9] * d[24] * d[35] + 2 * d[10] * d[25] * d[35] +
        2 * d[11] * d[26] * d[35];
    coeffs[155] = 2 * d[18] * d[20] * d[27] + 2 * d[21] * d[23] * d[27] + 2 * d[24] * d[26] * d[27] +
                  2 * d[19] * d[20] * d[28] + 2 * d[22] * d[23] * d[28] + 2 * d[25] * d[26] * d[28] +
                  std::pow(d[18], 2) * d[29] + std::pow(d[19], 2) * d[29] + 3 * std::pow(d[20], 2) * d[29] -
                  std::pow(d[21], 2) * d[29] - std::pow(d[22], 2) * d[29] + std::pow(d[23], 2) * d[29] -
                  std::pow(d[24], 2) * d[29] - std::pow(d[25], 2) * d[29] + std::pow(d[26], 2) * d[29] -
                  2 * d[20] * d[21] * d[30] + 2 * d[18] * d[23] * d[30] - 2 * d[20] * d[22] * d[31] +
                  2 * d[19] * d[23] * d[31] + 2 * d[18] * d[21] * d[32] + 2 * d[19] * d[22] * d[32] +
                  2 * d[20] * d[23] * d[32] - 2 * d[20] * d[24] * d[33] + 2 * d[18] * d[26] * d[33] -
                  2 * d[20] * d[25] * d[34] + 2 * d[19] * d[26] * d[34] + 2 * d[18] * d[24] * d[35] +
                  2 * d[19] * d[25] * d[35] + 2 * d[20] * d[26] * d[35];
    coeffs[156] =
        d[2] * std::pow(d[27], 2) + d[2] * std::pow(d[28], 2) + 2 * d[0] * d[27] * d[29] + 2 * d[1] * d[28] * d[29] +
        3 * d[2] * std::pow(d[29], 2) + 2 * d[5] * d[27] * d[30] - 2 * d[3] * d[29] * d[30] -
        d[2] * std::pow(d[30], 2) + 2 * d[5] * d[28] * d[31] - 2 * d[4] * d[29] * d[31] - d[2] * std::pow(d[31], 2) +
        2 * d[3] * d[27] * d[32] + 2 * d[4] * d[28] * d[32] + 2 * d[5] * d[29] * d[32] + 2 * d[0] * d[30] * d[32] +
        2 * d[1] * d[31] * d[32] + d[2] * std::pow(d[32], 2) + 2 * d[8] * d[27] * d[33] - 2 * d[6] * d[29] * d[33] -
        d[2] * std::pow(d[33], 2) + 2 * d[8] * d[28] * d[34] - 2 * d[7] * d[29] * d[34] - d[2] * std::pow(d[34], 2) +
        2 * d[6] * d[27] * d[35] + 2 * d[7] * d[28] * d[35] + 2 * d[8] * d[29] * d[35] + 2 * d[0] * d[33] * d[35] +
        2 * d[1] * d[34] * d[35] + d[2] * std::pow(d[35], 2);
    coeffs[157] =
        d[11] * std::pow(d[27], 2) + d[11] * std::pow(d[28], 2) + 2 * d[9] * d[27] * d[29] + 2 * d[10] * d[28] * d[29] +
        3 * d[11] * std::pow(d[29], 2) + 2 * d[14] * d[27] * d[30] - 2 * d[12] * d[29] * d[30] -
        d[11] * std::pow(d[30], 2) + 2 * d[14] * d[28] * d[31] - 2 * d[13] * d[29] * d[31] -
        d[11] * std::pow(d[31], 2) + 2 * d[12] * d[27] * d[32] + 2 * d[13] * d[28] * d[32] + 2 * d[14] * d[29] * d[32] +
        2 * d[9] * d[30] * d[32] + 2 * d[10] * d[31] * d[32] + d[11] * std::pow(d[32], 2) + 2 * d[17] * d[27] * d[33] -
        2 * d[15] * d[29] * d[33] - d[11] * std::pow(d[33], 2) + 2 * d[17] * d[28] * d[34] - 2 * d[16] * d[29] * d[34] -
        d[11] * std::pow(d[34], 2) + 2 * d[15] * d[27] * d[35] + 2 * d[16] * d[28] * d[35] + 2 * d[17] * d[29] * d[35] +
        2 * d[9] * d[33] * d[35] + 2 * d[10] * d[34] * d[35] + d[11] * std::pow(d[35], 2);
    coeffs[158] =
        d[20] * std::pow(d[27], 2) + d[20] * std::pow(d[28], 2) + 2 * d[18] * d[27] * d[29] +
        2 * d[19] * d[28] * d[29] + 3 * d[20] * std::pow(d[29], 2) + 2 * d[23] * d[27] * d[30] -
        2 * d[21] * d[29] * d[30] - d[20] * std::pow(d[30], 2) + 2 * d[23] * d[28] * d[31] - 2 * d[22] * d[29] * d[31] -
        d[20] * std::pow(d[31], 2) + 2 * d[21] * d[27] * d[32] + 2 * d[22] * d[28] * d[32] + 2 * d[23] * d[29] * d[32] +
        2 * d[18] * d[30] * d[32] + 2 * d[19] * d[31] * d[32] + d[20] * std::pow(d[32], 2) + 2 * d[26] * d[27] * d[33] -
        2 * d[24] * d[29] * d[33] - d[20] * std::pow(d[33], 2) + 2 * d[26] * d[28] * d[34] - 2 * d[25] * d[29] * d[34] -
        d[20] * std::pow(d[34], 2) + 2 * d[24] * d[27] * d[35] + 2 * d[25] * d[28] * d[35] + 2 * d[26] * d[29] * d[35] +
        2 * d[18] * d[33] * d[35] + 2 * d[19] * d[34] * d[35] + d[20] * std::pow(d[35], 2);
    coeffs[159] = std::pow(d[27], 2) * d[29] + std::pow(d[28], 2) * d[29] + std::pow(d[29], 3) -
                  d[29] * std::pow(d[30], 2) - d[29] * std::pow(d[31], 2) + 2 * d[27] * d[30] * d[32] +
                  2 * d[28] * d[31] * d[32] + d[29] * std::pow(d[32], 2) - d[29] * std::pow(d[33], 2) -
                  d[29] * std::pow(d[34], 2) + 2 * d[27] * d[33] * d[35] + 2 * d[28] * d[34] * d[35] +
                  d[29] * std::pow(d[35], 2);
    coeffs[160] =
        2 * d[0] * d[2] * d[36] + 2 * d[3] * d[5] * d[36] + 2 * d[6] * d[8] * d[36] + 2 * d[1] * d[2] * d[37] +
        2 * d[4] * d[5] * d[37] + 2 * d[7] * d[8] * d[37] + std::pow(d[0], 2) * d[38] + std::pow(d[1], 2) * d[38] +
        3 * std::pow(d[2], 2) * d[38] - std::pow(d[3], 2) * d[38] - std::pow(d[4], 2) * d[38] +
        std::pow(d[5], 2) * d[38] - std::pow(d[6], 2) * d[38] - std::pow(d[7], 2) * d[38] + std::pow(d[8], 2) * d[38] -
        2 * d[2] * d[3] * d[39] + 2 * d[0] * d[5] * d[39] - 2 * d[2] * d[4] * d[40] + 2 * d[1] * d[5] * d[40] +
        2 * d[0] * d[3] * d[41] + 2 * d[1] * d[4] * d[41] + 2 * d[2] * d[5] * d[41] - 2 * d[2] * d[6] * d[42] +
        2 * d[0] * d[8] * d[42] - 2 * d[2] * d[7] * d[43] + 2 * d[1] * d[8] * d[43] + 2 * d[0] * d[6] * d[44] +
        2 * d[1] * d[7] * d[44] + 2 * d[2] * d[8] * d[44];
    coeffs[161] =
        2 * d[2] * d[9] * d[36] + 2 * d[0] * d[11] * d[36] + 2 * d[5] * d[12] * d[36] + 2 * d[3] * d[14] * d[36] +
        2 * d[8] * d[15] * d[36] + 2 * d[6] * d[17] * d[36] + 2 * d[2] * d[10] * d[37] + 2 * d[1] * d[11] * d[37] +
        2 * d[5] * d[13] * d[37] + 2 * d[4] * d[14] * d[37] + 2 * d[8] * d[16] * d[37] + 2 * d[7] * d[17] * d[37] +
        2 * d[0] * d[9] * d[38] + 2 * d[1] * d[10] * d[38] + 6 * d[2] * d[11] * d[38] - 2 * d[3] * d[12] * d[38] -
        2 * d[4] * d[13] * d[38] + 2 * d[5] * d[14] * d[38] - 2 * d[6] * d[15] * d[38] - 2 * d[7] * d[16] * d[38] +
        2 * d[8] * d[17] * d[38] + 2 * d[5] * d[9] * d[39] - 2 * d[3] * d[11] * d[39] - 2 * d[2] * d[12] * d[39] +
        2 * d[0] * d[14] * d[39] + 2 * d[5] * d[10] * d[40] - 2 * d[4] * d[11] * d[40] - 2 * d[2] * d[13] * d[40] +
        2 * d[1] * d[14] * d[40] + 2 * d[3] * d[9] * d[41] + 2 * d[4] * d[10] * d[41] + 2 * d[5] * d[11] * d[41] +
        2 * d[0] * d[12] * d[41] + 2 * d[1] * d[13] * d[41] + 2 * d[2] * d[14] * d[41] + 2 * d[8] * d[9] * d[42] -
        2 * d[6] * d[11] * d[42] - 2 * d[2] * d[15] * d[42] + 2 * d[0] * d[17] * d[42] + 2 * d[8] * d[10] * d[43] -
        2 * d[7] * d[11] * d[43] - 2 * d[2] * d[16] * d[43] + 2 * d[1] * d[17] * d[43] + 2 * d[6] * d[9] * d[44] +
        2 * d[7] * d[10] * d[44] + 2 * d[8] * d[11] * d[44] + 2 * d[0] * d[15] * d[44] + 2 * d[1] * d[16] * d[44] +
        2 * d[2] * d[17] * d[44];
    coeffs[162] =
        2 * d[9] * d[11] * d[36] + 2 * d[12] * d[14] * d[36] + 2 * d[15] * d[17] * d[36] + 2 * d[10] * d[11] * d[37] +
        2 * d[13] * d[14] * d[37] + 2 * d[16] * d[17] * d[37] + std::pow(d[9], 2) * d[38] + std::pow(d[10], 2) * d[38] +
        3 * std::pow(d[11], 2) * d[38] - std::pow(d[12], 2) * d[38] - std::pow(d[13], 2) * d[38] +
        std::pow(d[14], 2) * d[38] - std::pow(d[15], 2) * d[38] - std::pow(d[16], 2) * d[38] +
        std::pow(d[17], 2) * d[38] - 2 * d[11] * d[12] * d[39] + 2 * d[9] * d[14] * d[39] - 2 * d[11] * d[13] * d[40] +
        2 * d[10] * d[14] * d[40] + 2 * d[9] * d[12] * d[41] + 2 * d[10] * d[13] * d[41] + 2 * d[11] * d[14] * d[41] -
        2 * d[11] * d[15] * d[42] + 2 * d[9] * d[17] * d[42] - 2 * d[11] * d[16] * d[43] + 2 * d[10] * d[17] * d[43] +
        2 * d[9] * d[15] * d[44] + 2 * d[10] * d[16] * d[44] + 2 * d[11] * d[17] * d[44];
    coeffs[163] =
        2 * d[2] * d[18] * d[36] + 2 * d[0] * d[20] * d[36] + 2 * d[5] * d[21] * d[36] + 2 * d[3] * d[23] * d[36] +
        2 * d[8] * d[24] * d[36] + 2 * d[6] * d[26] * d[36] + 2 * d[2] * d[19] * d[37] + 2 * d[1] * d[20] * d[37] +
        2 * d[5] * d[22] * d[37] + 2 * d[4] * d[23] * d[37] + 2 * d[8] * d[25] * d[37] + 2 * d[7] * d[26] * d[37] +
        2 * d[0] * d[18] * d[38] + 2 * d[1] * d[19] * d[38] + 6 * d[2] * d[20] * d[38] - 2 * d[3] * d[21] * d[38] -
        2 * d[4] * d[22] * d[38] + 2 * d[5] * d[23] * d[38] - 2 * d[6] * d[24] * d[38] - 2 * d[7] * d[25] * d[38] +
        2 * d[8] * d[26] * d[38] + 2 * d[5] * d[18] * d[39] - 2 * d[3] * d[20] * d[39] - 2 * d[2] * d[21] * d[39] +
        2 * d[0] * d[23] * d[39] + 2 * d[5] * d[19] * d[40] - 2 * d[4] * d[20] * d[40] - 2 * d[2] * d[22] * d[40] +
        2 * d[1] * d[23] * d[40] + 2 * d[3] * d[18] * d[41] + 2 * d[4] * d[19] * d[41] + 2 * d[5] * d[20] * d[41] +
        2 * d[0] * d[21] * d[41] + 2 * d[1] * d[22] * d[41] + 2 * d[2] * d[23] * d[41] + 2 * d[8] * d[18] * d[42] -
        2 * d[6] * d[20] * d[42] - 2 * d[2] * d[24] * d[42] + 2 * d[0] * d[26] * d[42] + 2 * d[8] * d[19] * d[43] -
        2 * d[7] * d[20] * d[43] - 2 * d[2] * d[25] * d[43] + 2 * d[1] * d[26] * d[43] + 2 * d[6] * d[18] * d[44] +
        2 * d[7] * d[19] * d[44] + 2 * d[8] * d[20] * d[44] + 2 * d[0] * d[24] * d[44] + 2 * d[1] * d[25] * d[44] +
        2 * d[2] * d[26] * d[44];
    coeffs[164] =
        2 * d[11] * d[18] * d[36] + 2 * d[9] * d[20] * d[36] + 2 * d[14] * d[21] * d[36] + 2 * d[12] * d[23] * d[36] +
        2 * d[17] * d[24] * d[36] + 2 * d[15] * d[26] * d[36] + 2 * d[11] * d[19] * d[37] + 2 * d[10] * d[20] * d[37] +
        2 * d[14] * d[22] * d[37] + 2 * d[13] * d[23] * d[37] + 2 * d[17] * d[25] * d[37] + 2 * d[16] * d[26] * d[37] +
        2 * d[9] * d[18] * d[38] + 2 * d[10] * d[19] * d[38] + 6 * d[11] * d[20] * d[38] - 2 * d[12] * d[21] * d[38] -
        2 * d[13] * d[22] * d[38] + 2 * d[14] * d[23] * d[38] - 2 * d[15] * d[24] * d[38] - 2 * d[16] * d[25] * d[38] +
        2 * d[17] * d[26] * d[38] + 2 * d[14] * d[18] * d[39] - 2 * d[12] * d[20] * d[39] - 2 * d[11] * d[21] * d[39] +
        2 * d[9] * d[23] * d[39] + 2 * d[14] * d[19] * d[40] - 2 * d[13] * d[20] * d[40] - 2 * d[11] * d[22] * d[40] +
        2 * d[10] * d[23] * d[40] + 2 * d[12] * d[18] * d[41] + 2 * d[13] * d[19] * d[41] + 2 * d[14] * d[20] * d[41] +
        2 * d[9] * d[21] * d[41] + 2 * d[10] * d[22] * d[41] + 2 * d[11] * d[23] * d[41] + 2 * d[17] * d[18] * d[42] -
        2 * d[15] * d[20] * d[42] - 2 * d[11] * d[24] * d[42] + 2 * d[9] * d[26] * d[42] + 2 * d[17] * d[19] * d[43] -
        2 * d[16] * d[20] * d[43] - 2 * d[11] * d[25] * d[43] + 2 * d[10] * d[26] * d[43] + 2 * d[15] * d[18] * d[44] +
        2 * d[16] * d[19] * d[44] + 2 * d[17] * d[20] * d[44] + 2 * d[9] * d[24] * d[44] + 2 * d[10] * d[25] * d[44] +
        2 * d[11] * d[26] * d[44];
    coeffs[165] = 2 * d[18] * d[20] * d[36] + 2 * d[21] * d[23] * d[36] + 2 * d[24] * d[26] * d[36] +
                  2 * d[19] * d[20] * d[37] + 2 * d[22] * d[23] * d[37] + 2 * d[25] * d[26] * d[37] +
                  std::pow(d[18], 2) * d[38] + std::pow(d[19], 2) * d[38] + 3 * std::pow(d[20], 2) * d[38] -
                  std::pow(d[21], 2) * d[38] - std::pow(d[22], 2) * d[38] + std::pow(d[23], 2) * d[38] -
                  std::pow(d[24], 2) * d[38] - std::pow(d[25], 2) * d[38] + std::pow(d[26], 2) * d[38] -
                  2 * d[20] * d[21] * d[39] + 2 * d[18] * d[23] * d[39] - 2 * d[20] * d[22] * d[40] +
                  2 * d[19] * d[23] * d[40] + 2 * d[18] * d[21] * d[41] + 2 * d[19] * d[22] * d[41] +
                  2 * d[20] * d[23] * d[41] - 2 * d[20] * d[24] * d[42] + 2 * d[18] * d[26] * d[42] -
                  2 * d[20] * d[25] * d[43] + 2 * d[19] * d[26] * d[43] + 2 * d[18] * d[24] * d[44] +
                  2 * d[19] * d[25] * d[44] + 2 * d[20] * d[26] * d[44];
    coeffs[166] =
        2 * d[2] * d[27] * d[36] + 2 * d[0] * d[29] * d[36] + 2 * d[5] * d[30] * d[36] + 2 * d[3] * d[32] * d[36] +
        2 * d[8] * d[33] * d[36] + 2 * d[6] * d[35] * d[36] + 2 * d[2] * d[28] * d[37] + 2 * d[1] * d[29] * d[37] +
        2 * d[5] * d[31] * d[37] + 2 * d[4] * d[32] * d[37] + 2 * d[8] * d[34] * d[37] + 2 * d[7] * d[35] * d[37] +
        2 * d[0] * d[27] * d[38] + 2 * d[1] * d[28] * d[38] + 6 * d[2] * d[29] * d[38] - 2 * d[3] * d[30] * d[38] -
        2 * d[4] * d[31] * d[38] + 2 * d[5] * d[32] * d[38] - 2 * d[6] * d[33] * d[38] - 2 * d[7] * d[34] * d[38] +
        2 * d[8] * d[35] * d[38] + 2 * d[5] * d[27] * d[39] - 2 * d[3] * d[29] * d[39] - 2 * d[2] * d[30] * d[39] +
        2 * d[0] * d[32] * d[39] + 2 * d[5] * d[28] * d[40] - 2 * d[4] * d[29] * d[40] - 2 * d[2] * d[31] * d[40] +
        2 * d[1] * d[32] * d[40] + 2 * d[3] * d[27] * d[41] + 2 * d[4] * d[28] * d[41] + 2 * d[5] * d[29] * d[41] +
        2 * d[0] * d[30] * d[41] + 2 * d[1] * d[31] * d[41] + 2 * d[2] * d[32] * d[41] + 2 * d[8] * d[27] * d[42] -
        2 * d[6] * d[29] * d[42] - 2 * d[2] * d[33] * d[42] + 2 * d[0] * d[35] * d[42] + 2 * d[8] * d[28] * d[43] -
        2 * d[7] * d[29] * d[43] - 2 * d[2] * d[34] * d[43] + 2 * d[1] * d[35] * d[43] + 2 * d[6] * d[27] * d[44] +
        2 * d[7] * d[28] * d[44] + 2 * d[8] * d[29] * d[44] + 2 * d[0] * d[33] * d[44] + 2 * d[1] * d[34] * d[44] +
        2 * d[2] * d[35] * d[44];
    coeffs[167] =
        2 * d[11] * d[27] * d[36] + 2 * d[9] * d[29] * d[36] + 2 * d[14] * d[30] * d[36] + 2 * d[12] * d[32] * d[36] +
        2 * d[17] * d[33] * d[36] + 2 * d[15] * d[35] * d[36] + 2 * d[11] * d[28] * d[37] + 2 * d[10] * d[29] * d[37] +
        2 * d[14] * d[31] * d[37] + 2 * d[13] * d[32] * d[37] + 2 * d[17] * d[34] * d[37] + 2 * d[16] * d[35] * d[37] +
        2 * d[9] * d[27] * d[38] + 2 * d[10] * d[28] * d[38] + 6 * d[11] * d[29] * d[38] - 2 * d[12] * d[30] * d[38] -
        2 * d[13] * d[31] * d[38] + 2 * d[14] * d[32] * d[38] - 2 * d[15] * d[33] * d[38] - 2 * d[16] * d[34] * d[38] +
        2 * d[17] * d[35] * d[38] + 2 * d[14] * d[27] * d[39] - 2 * d[12] * d[29] * d[39] - 2 * d[11] * d[30] * d[39] +
        2 * d[9] * d[32] * d[39] + 2 * d[14] * d[28] * d[40] - 2 * d[13] * d[29] * d[40] - 2 * d[11] * d[31] * d[40] +
        2 * d[10] * d[32] * d[40] + 2 * d[12] * d[27] * d[41] + 2 * d[13] * d[28] * d[41] + 2 * d[14] * d[29] * d[41] +
        2 * d[9] * d[30] * d[41] + 2 * d[10] * d[31] * d[41] + 2 * d[11] * d[32] * d[41] + 2 * d[17] * d[27] * d[42] -
        2 * d[15] * d[29] * d[42] - 2 * d[11] * d[33] * d[42] + 2 * d[9] * d[35] * d[42] + 2 * d[17] * d[28] * d[43] -
        2 * d[16] * d[29] * d[43] - 2 * d[11] * d[34] * d[43] + 2 * d[10] * d[35] * d[43] + 2 * d[15] * d[27] * d[44] +
        2 * d[16] * d[28] * d[44] + 2 * d[17] * d[29] * d[44] + 2 * d[9] * d[33] * d[44] + 2 * d[10] * d[34] * d[44] +
        2 * d[11] * d[35] * d[44];
    coeffs[168] =
        2 * d[20] * d[27] * d[36] + 2 * d[18] * d[29] * d[36] + 2 * d[23] * d[30] * d[36] + 2 * d[21] * d[32] * d[36] +
        2 * d[26] * d[33] * d[36] + 2 * d[24] * d[35] * d[36] + 2 * d[20] * d[28] * d[37] + 2 * d[19] * d[29] * d[37] +
        2 * d[23] * d[31] * d[37] + 2 * d[22] * d[32] * d[37] + 2 * d[26] * d[34] * d[37] + 2 * d[25] * d[35] * d[37] +
        2 * d[18] * d[27] * d[38] + 2 * d[19] * d[28] * d[38] + 6 * d[20] * d[29] * d[38] - 2 * d[21] * d[30] * d[38] -
        2 * d[22] * d[31] * d[38] + 2 * d[23] * d[32] * d[38] - 2 * d[24] * d[33] * d[38] - 2 * d[25] * d[34] * d[38] +
        2 * d[26] * d[35] * d[38] + 2 * d[23] * d[27] * d[39] - 2 * d[21] * d[29] * d[39] - 2 * d[20] * d[30] * d[39] +
        2 * d[18] * d[32] * d[39] + 2 * d[23] * d[28] * d[40] - 2 * d[22] * d[29] * d[40] - 2 * d[20] * d[31] * d[40] +
        2 * d[19] * d[32] * d[40] + 2 * d[21] * d[27] * d[41] + 2 * d[22] * d[28] * d[41] + 2 * d[23] * d[29] * d[41] +
        2 * d[18] * d[30] * d[41] + 2 * d[19] * d[31] * d[41] + 2 * d[20] * d[32] * d[41] + 2 * d[26] * d[27] * d[42] -
        2 * d[24] * d[29] * d[42] - 2 * d[20] * d[33] * d[42] + 2 * d[18] * d[35] * d[42] + 2 * d[26] * d[28] * d[43] -
        2 * d[25] * d[29] * d[43] - 2 * d[20] * d[34] * d[43] + 2 * d[19] * d[35] * d[43] + 2 * d[24] * d[27] * d[44] +
        2 * d[25] * d[28] * d[44] + 2 * d[26] * d[29] * d[44] + 2 * d[18] * d[33] * d[44] + 2 * d[19] * d[34] * d[44] +
        2 * d[20] * d[35] * d[44];
    coeffs[169] = 2 * d[27] * d[29] * d[36] + 2 * d[30] * d[32] * d[36] + 2 * d[33] * d[35] * d[36] +
                  2 * d[28] * d[29] * d[37] + 2 * d[31] * d[32] * d[37] + 2 * d[34] * d[35] * d[37] +
                  std::pow(d[27], 2) * d[38] + std::pow(d[28], 2) * d[38] + 3 * std::pow(d[29], 2) * d[38] -
                  std::pow(d[30], 2) * d[38] - std::pow(d[31], 2) * d[38] + std::pow(d[32], 2) * d[38] -
                  std::pow(d[33], 2) * d[38] - std::pow(d[34], 2) * d[38] + std::pow(d[35], 2) * d[38] -
                  2 * d[29] * d[30] * d[39] + 2 * d[27] * d[32] * d[39] - 2 * d[29] * d[31] * d[40] +
                  2 * d[28] * d[32] * d[40] + 2 * d[27] * d[30] * d[41] + 2 * d[28] * d[31] * d[41] +
                  2 * d[29] * d[32] * d[41] - 2 * d[29] * d[33] * d[42] + 2 * d[27] * d[35] * d[42] -
                  2 * d[29] * d[34] * d[43] + 2 * d[28] * d[35] * d[43] + 2 * d[27] * d[33] * d[44] +
                  2 * d[28] * d[34] * d[44] + 2 * d[29] * d[35] * d[44];
    coeffs[170] =
        d[2] * std::pow(d[36], 2) + d[2] * std::pow(d[37], 2) + 2 * d[0] * d[36] * d[38] + 2 * d[1] * d[37] * d[38] +
        3 * d[2] * std::pow(d[38], 2) + 2 * d[5] * d[36] * d[39] - 2 * d[3] * d[38] * d[39] -
        d[2] * std::pow(d[39], 2) + 2 * d[5] * d[37] * d[40] - 2 * d[4] * d[38] * d[40] - d[2] * std::pow(d[40], 2) +
        2 * d[3] * d[36] * d[41] + 2 * d[4] * d[37] * d[41] + 2 * d[5] * d[38] * d[41] + 2 * d[0] * d[39] * d[41] +
        2 * d[1] * d[40] * d[41] + d[2] * std::pow(d[41], 2) + 2 * d[8] * d[36] * d[42] - 2 * d[6] * d[38] * d[42] -
        d[2] * std::pow(d[42], 2) + 2 * d[8] * d[37] * d[43] - 2 * d[7] * d[38] * d[43] - d[2] * std::pow(d[43], 2) +
        2 * d[6] * d[36] * d[44] + 2 * d[7] * d[37] * d[44] + 2 * d[8] * d[38] * d[44] + 2 * d[0] * d[42] * d[44] +
        2 * d[1] * d[43] * d[44] + d[2] * std::pow(d[44], 2);
    coeffs[171] =
        d[11] * std::pow(d[36], 2) + d[11] * std::pow(d[37], 2) + 2 * d[9] * d[36] * d[38] + 2 * d[10] * d[37] * d[38] +
        3 * d[11] * std::pow(d[38], 2) + 2 * d[14] * d[36] * d[39] - 2 * d[12] * d[38] * d[39] -
        d[11] * std::pow(d[39], 2) + 2 * d[14] * d[37] * d[40] - 2 * d[13] * d[38] * d[40] -
        d[11] * std::pow(d[40], 2) + 2 * d[12] * d[36] * d[41] + 2 * d[13] * d[37] * d[41] + 2 * d[14] * d[38] * d[41] +
        2 * d[9] * d[39] * d[41] + 2 * d[10] * d[40] * d[41] + d[11] * std::pow(d[41], 2) + 2 * d[17] * d[36] * d[42] -
        2 * d[15] * d[38] * d[42] - d[11] * std::pow(d[42], 2) + 2 * d[17] * d[37] * d[43] - 2 * d[16] * d[38] * d[43] -
        d[11] * std::pow(d[43], 2) + 2 * d[15] * d[36] * d[44] + 2 * d[16] * d[37] * d[44] + 2 * d[17] * d[38] * d[44] +
        2 * d[9] * d[42] * d[44] + 2 * d[10] * d[43] * d[44] + d[11] * std::pow(d[44], 2);
    coeffs[172] =
        d[20] * std::pow(d[36], 2) + d[20] * std::pow(d[37], 2) + 2 * d[18] * d[36] * d[38] +
        2 * d[19] * d[37] * d[38] + 3 * d[20] * std::pow(d[38], 2) + 2 * d[23] * d[36] * d[39] -
        2 * d[21] * d[38] * d[39] - d[20] * std::pow(d[39], 2) + 2 * d[23] * d[37] * d[40] - 2 * d[22] * d[38] * d[40] -
        d[20] * std::pow(d[40], 2) + 2 * d[21] * d[36] * d[41] + 2 * d[22] * d[37] * d[41] + 2 * d[23] * d[38] * d[41] +
        2 * d[18] * d[39] * d[41] + 2 * d[19] * d[40] * d[41] + d[20] * std::pow(d[41], 2) + 2 * d[26] * d[36] * d[42] -
        2 * d[24] * d[38] * d[42] - d[20] * std::pow(d[42], 2) + 2 * d[26] * d[37] * d[43] - 2 * d[25] * d[38] * d[43] -
        d[20] * std::pow(d[43], 2) + 2 * d[24] * d[36] * d[44] + 2 * d[25] * d[37] * d[44] + 2 * d[26] * d[38] * d[44] +
        2 * d[18] * d[42] * d[44] + 2 * d[19] * d[43] * d[44] + d[20] * std::pow(d[44], 2);
    coeffs[173] =
        d[29] * std::pow(d[36], 2) + d[29] * std::pow(d[37], 2) + 2 * d[27] * d[36] * d[38] +
        2 * d[28] * d[37] * d[38] + 3 * d[29] * std::pow(d[38], 2) + 2 * d[32] * d[36] * d[39] -
        2 * d[30] * d[38] * d[39] - d[29] * std::pow(d[39], 2) + 2 * d[32] * d[37] * d[40] - 2 * d[31] * d[38] * d[40] -
        d[29] * std::pow(d[40], 2) + 2 * d[30] * d[36] * d[41] + 2 * d[31] * d[37] * d[41] + 2 * d[32] * d[38] * d[41] +
        2 * d[27] * d[39] * d[41] + 2 * d[28] * d[40] * d[41] + d[29] * std::pow(d[41], 2) + 2 * d[35] * d[36] * d[42] -
        2 * d[33] * d[38] * d[42] - d[29] * std::pow(d[42], 2) + 2 * d[35] * d[37] * d[43] - 2 * d[34] * d[38] * d[43] -
        d[29] * std::pow(d[43], 2) + 2 * d[33] * d[36] * d[44] + 2 * d[34] * d[37] * d[44] + 2 * d[35] * d[38] * d[44] +
        2 * d[27] * d[42] * d[44] + 2 * d[28] * d[43] * d[44] + d[29] * std::pow(d[44], 2);
    coeffs[174] = std::pow(d[36], 2) * d[38] + std::pow(d[37], 2) * d[38] + std::pow(d[38], 3) -
                  d[38] * std::pow(d[39], 2) - d[38] * std::pow(d[40], 2) + 2 * d[36] * d[39] * d[41] +
                  2 * d[37] * d[40] * d[41] + d[38] * std::pow(d[41], 2) - d[38] * std::pow(d[42], 2) -
                  d[38] * std::pow(d[43], 2) + 2 * d[36] * d[42] * d[44] + 2 * d[37] * d[43] * d[44] +
                  d[38] * std::pow(d[44], 2);
    coeffs[175] = std::pow(d[0], 2) * d[3] - std::pow(d[1], 2) * d[3] - std::pow(d[2], 2) * d[3] + std::pow(d[3], 3) +
                  2 * d[0] * d[1] * d[4] + d[3] * std::pow(d[4], 2) + 2 * d[0] * d[2] * d[5] +
                  d[3] * std::pow(d[5], 2) + d[3] * std::pow(d[6], 2) + 2 * d[4] * d[6] * d[7] -
                  d[3] * std::pow(d[7], 2) + 2 * d[5] * d[6] * d[8] - d[3] * std::pow(d[8], 2);
    coeffs[176] = 2 * d[0] * d[3] * d[9] + 2 * d[1] * d[4] * d[9] + 2 * d[2] * d[5] * d[9] - 2 * d[1] * d[3] * d[10] +
                  2 * d[0] * d[4] * d[10] - 2 * d[2] * d[3] * d[11] + 2 * d[0] * d[5] * d[11] +
                  std::pow(d[0], 2) * d[12] - std::pow(d[1], 2) * d[12] - std::pow(d[2], 2) * d[12] +
                  3 * std::pow(d[3], 2) * d[12] + std::pow(d[4], 2) * d[12] + std::pow(d[5], 2) * d[12] +
                  std::pow(d[6], 2) * d[12] - std::pow(d[7], 2) * d[12] - std::pow(d[8], 2) * d[12] +
                  2 * d[0] * d[1] * d[13] + 2 * d[3] * d[4] * d[13] + 2 * d[6] * d[7] * d[13] +
                  2 * d[0] * d[2] * d[14] + 2 * d[3] * d[5] * d[14] + 2 * d[6] * d[8] * d[14] +
                  2 * d[3] * d[6] * d[15] + 2 * d[4] * d[7] * d[15] + 2 * d[5] * d[8] * d[15] +
                  2 * d[4] * d[6] * d[16] - 2 * d[3] * d[7] * d[16] + 2 * d[5] * d[6] * d[17] - 2 * d[3] * d[8] * d[17];
    coeffs[177] =
        d[3] * std::pow(d[9], 2) + 2 * d[4] * d[9] * d[10] - d[3] * std::pow(d[10], 2) + 2 * d[5] * d[9] * d[11] -
        d[3] * std::pow(d[11], 2) + 2 * d[0] * d[9] * d[12] - 2 * d[1] * d[10] * d[12] - 2 * d[2] * d[11] * d[12] +
        3 * d[3] * std::pow(d[12], 2) + 2 * d[1] * d[9] * d[13] + 2 * d[0] * d[10] * d[13] + 2 * d[4] * d[12] * d[13] +
        d[3] * std::pow(d[13], 2) + 2 * d[2] * d[9] * d[14] + 2 * d[0] * d[11] * d[14] + 2 * d[5] * d[12] * d[14] +
        d[3] * std::pow(d[14], 2) + 2 * d[6] * d[12] * d[15] + 2 * d[7] * d[13] * d[15] + 2 * d[8] * d[14] * d[15] +
        d[3] * std::pow(d[15], 2) - 2 * d[7] * d[12] * d[16] + 2 * d[6] * d[13] * d[16] + 2 * d[4] * d[15] * d[16] -
        d[3] * std::pow(d[16], 2) - 2 * d[8] * d[12] * d[17] + 2 * d[6] * d[14] * d[17] + 2 * d[5] * d[15] * d[17] -
        d[3] * std::pow(d[17], 2);
    coeffs[178] = std::pow(d[9], 2) * d[12] - std::pow(d[10], 2) * d[12] - std::pow(d[11], 2) * d[12] +
                  std::pow(d[12], 3) + 2 * d[9] * d[10] * d[13] + d[12] * std::pow(d[13], 2) +
                  2 * d[9] * d[11] * d[14] + d[12] * std::pow(d[14], 2) + d[12] * std::pow(d[15], 2) +
                  2 * d[13] * d[15] * d[16] - d[12] * std::pow(d[16], 2) + 2 * d[14] * d[15] * d[17] -
                  d[12] * std::pow(d[17], 2);
    coeffs[179] =
        2 * d[0] * d[3] * d[18] + 2 * d[1] * d[4] * d[18] + 2 * d[2] * d[5] * d[18] - 2 * d[1] * d[3] * d[19] +
        2 * d[0] * d[4] * d[19] - 2 * d[2] * d[3] * d[20] + 2 * d[0] * d[5] * d[20] + std::pow(d[0], 2) * d[21] -
        std::pow(d[1], 2) * d[21] - std::pow(d[2], 2) * d[21] + 3 * std::pow(d[3], 2) * d[21] +
        std::pow(d[4], 2) * d[21] + std::pow(d[5], 2) * d[21] + std::pow(d[6], 2) * d[21] - std::pow(d[7], 2) * d[21] -
        std::pow(d[8], 2) * d[21] + 2 * d[0] * d[1] * d[22] + 2 * d[3] * d[4] * d[22] + 2 * d[6] * d[7] * d[22] +
        2 * d[0] * d[2] * d[23] + 2 * d[3] * d[5] * d[23] + 2 * d[6] * d[8] * d[23] + 2 * d[3] * d[6] * d[24] +
        2 * d[4] * d[7] * d[24] + 2 * d[5] * d[8] * d[24] + 2 * d[4] * d[6] * d[25] - 2 * d[3] * d[7] * d[25] +
        2 * d[5] * d[6] * d[26] - 2 * d[3] * d[8] * d[26];
    coeffs[180] =
        2 * d[3] * d[9] * d[18] + 2 * d[4] * d[10] * d[18] + 2 * d[5] * d[11] * d[18] + 2 * d[0] * d[12] * d[18] +
        2 * d[1] * d[13] * d[18] + 2 * d[2] * d[14] * d[18] + 2 * d[4] * d[9] * d[19] - 2 * d[3] * d[10] * d[19] -
        2 * d[1] * d[12] * d[19] + 2 * d[0] * d[13] * d[19] + 2 * d[5] * d[9] * d[20] - 2 * d[3] * d[11] * d[20] -
        2 * d[2] * d[12] * d[20] + 2 * d[0] * d[14] * d[20] + 2 * d[0] * d[9] * d[21] - 2 * d[1] * d[10] * d[21] -
        2 * d[2] * d[11] * d[21] + 6 * d[3] * d[12] * d[21] + 2 * d[4] * d[13] * d[21] + 2 * d[5] * d[14] * d[21] +
        2 * d[6] * d[15] * d[21] - 2 * d[7] * d[16] * d[21] - 2 * d[8] * d[17] * d[21] + 2 * d[1] * d[9] * d[22] +
        2 * d[0] * d[10] * d[22] + 2 * d[4] * d[12] * d[22] + 2 * d[3] * d[13] * d[22] + 2 * d[7] * d[15] * d[22] +
        2 * d[6] * d[16] * d[22] + 2 * d[2] * d[9] * d[23] + 2 * d[0] * d[11] * d[23] + 2 * d[5] * d[12] * d[23] +
        2 * d[3] * d[14] * d[23] + 2 * d[8] * d[15] * d[23] + 2 * d[6] * d[17] * d[23] + 2 * d[6] * d[12] * d[24] +
        2 * d[7] * d[13] * d[24] + 2 * d[8] * d[14] * d[24] + 2 * d[3] * d[15] * d[24] + 2 * d[4] * d[16] * d[24] +
        2 * d[5] * d[17] * d[24] - 2 * d[7] * d[12] * d[25] + 2 * d[6] * d[13] * d[25] + 2 * d[4] * d[15] * d[25] -
        2 * d[3] * d[16] * d[25] - 2 * d[8] * d[12] * d[26] + 2 * d[6] * d[14] * d[26] + 2 * d[5] * d[15] * d[26] -
        2 * d[3] * d[17] * d[26];
    coeffs[181] =
        2 * d[9] * d[12] * d[18] + 2 * d[10] * d[13] * d[18] + 2 * d[11] * d[14] * d[18] - 2 * d[10] * d[12] * d[19] +
        2 * d[9] * d[13] * d[19] - 2 * d[11] * d[12] * d[20] + 2 * d[9] * d[14] * d[20] + std::pow(d[9], 2) * d[21] -
        std::pow(d[10], 2) * d[21] - std::pow(d[11], 2) * d[21] + 3 * std::pow(d[12], 2) * d[21] +
        std::pow(d[13], 2) * d[21] + std::pow(d[14], 2) * d[21] + std::pow(d[15], 2) * d[21] -
        std::pow(d[16], 2) * d[21] - std::pow(d[17], 2) * d[21] + 2 * d[9] * d[10] * d[22] + 2 * d[12] * d[13] * d[22] +
        2 * d[15] * d[16] * d[22] + 2 * d[9] * d[11] * d[23] + 2 * d[12] * d[14] * d[23] + 2 * d[15] * d[17] * d[23] +
        2 * d[12] * d[15] * d[24] + 2 * d[13] * d[16] * d[24] + 2 * d[14] * d[17] * d[24] + 2 * d[13] * d[15] * d[25] -
        2 * d[12] * d[16] * d[25] + 2 * d[14] * d[15] * d[26] - 2 * d[12] * d[17] * d[26];
    coeffs[182] =
        d[3] * std::pow(d[18], 2) + 2 * d[4] * d[18] * d[19] - d[3] * std::pow(d[19], 2) + 2 * d[5] * d[18] * d[20] -
        d[3] * std::pow(d[20], 2) + 2 * d[0] * d[18] * d[21] - 2 * d[1] * d[19] * d[21] - 2 * d[2] * d[20] * d[21] +
        3 * d[3] * std::pow(d[21], 2) + 2 * d[1] * d[18] * d[22] + 2 * d[0] * d[19] * d[22] + 2 * d[4] * d[21] * d[22] +
        d[3] * std::pow(d[22], 2) + 2 * d[2] * d[18] * d[23] + 2 * d[0] * d[20] * d[23] + 2 * d[5] * d[21] * d[23] +
        d[3] * std::pow(d[23], 2) + 2 * d[6] * d[21] * d[24] + 2 * d[7] * d[22] * d[24] + 2 * d[8] * d[23] * d[24] +
        d[3] * std::pow(d[24], 2) - 2 * d[7] * d[21] * d[25] + 2 * d[6] * d[22] * d[25] + 2 * d[4] * d[24] * d[25] -
        d[3] * std::pow(d[25], 2) - 2 * d[8] * d[21] * d[26] + 2 * d[6] * d[23] * d[26] + 2 * d[5] * d[24] * d[26] -
        d[3] * std::pow(d[26], 2);
    coeffs[183] =
        d[12] * std::pow(d[18], 2) + 2 * d[13] * d[18] * d[19] - d[12] * std::pow(d[19], 2) +
        2 * d[14] * d[18] * d[20] - d[12] * std::pow(d[20], 2) + 2 * d[9] * d[18] * d[21] - 2 * d[10] * d[19] * d[21] -
        2 * d[11] * d[20] * d[21] + 3 * d[12] * std::pow(d[21], 2) + 2 * d[10] * d[18] * d[22] +
        2 * d[9] * d[19] * d[22] + 2 * d[13] * d[21] * d[22] + d[12] * std::pow(d[22], 2) + 2 * d[11] * d[18] * d[23] +
        2 * d[9] * d[20] * d[23] + 2 * d[14] * d[21] * d[23] + d[12] * std::pow(d[23], 2) + 2 * d[15] * d[21] * d[24] +
        2 * d[16] * d[22] * d[24] + 2 * d[17] * d[23] * d[24] + d[12] * std::pow(d[24], 2) - 2 * d[16] * d[21] * d[25] +
        2 * d[15] * d[22] * d[25] + 2 * d[13] * d[24] * d[25] - d[12] * std::pow(d[25], 2) - 2 * d[17] * d[21] * d[26] +
        2 * d[15] * d[23] * d[26] + 2 * d[14] * d[24] * d[26] - d[12] * std::pow(d[26], 2);
    coeffs[184] = std::pow(d[18], 2) * d[21] - std::pow(d[19], 2) * d[21] - std::pow(d[20], 2) * d[21] +
                  std::pow(d[21], 3) + 2 * d[18] * d[19] * d[22] + d[21] * std::pow(d[22], 2) +
                  2 * d[18] * d[20] * d[23] + d[21] * std::pow(d[23], 2) + d[21] * std::pow(d[24], 2) +
                  2 * d[22] * d[24] * d[25] - d[21] * std::pow(d[25], 2) + 2 * d[23] * d[24] * d[26] -
                  d[21] * std::pow(d[26], 2);
    coeffs[185] =
        2 * d[0] * d[3] * d[27] + 2 * d[1] * d[4] * d[27] + 2 * d[2] * d[5] * d[27] - 2 * d[1] * d[3] * d[28] +
        2 * d[0] * d[4] * d[28] - 2 * d[2] * d[3] * d[29] + 2 * d[0] * d[5] * d[29] + std::pow(d[0], 2) * d[30] -
        std::pow(d[1], 2) * d[30] - std::pow(d[2], 2) * d[30] + 3 * std::pow(d[3], 2) * d[30] +
        std::pow(d[4], 2) * d[30] + std::pow(d[5], 2) * d[30] + std::pow(d[6], 2) * d[30] - std::pow(d[7], 2) * d[30] -
        std::pow(d[8], 2) * d[30] + 2 * d[0] * d[1] * d[31] + 2 * d[3] * d[4] * d[31] + 2 * d[6] * d[7] * d[31] +
        2 * d[0] * d[2] * d[32] + 2 * d[3] * d[5] * d[32] + 2 * d[6] * d[8] * d[32] + 2 * d[3] * d[6] * d[33] +
        2 * d[4] * d[7] * d[33] + 2 * d[5] * d[8] * d[33] + 2 * d[4] * d[6] * d[34] - 2 * d[3] * d[7] * d[34] +
        2 * d[5] * d[6] * d[35] - 2 * d[3] * d[8] * d[35];
    coeffs[186] =
        2 * d[3] * d[9] * d[27] + 2 * d[4] * d[10] * d[27] + 2 * d[5] * d[11] * d[27] + 2 * d[0] * d[12] * d[27] +
        2 * d[1] * d[13] * d[27] + 2 * d[2] * d[14] * d[27] + 2 * d[4] * d[9] * d[28] - 2 * d[3] * d[10] * d[28] -
        2 * d[1] * d[12] * d[28] + 2 * d[0] * d[13] * d[28] + 2 * d[5] * d[9] * d[29] - 2 * d[3] * d[11] * d[29] -
        2 * d[2] * d[12] * d[29] + 2 * d[0] * d[14] * d[29] + 2 * d[0] * d[9] * d[30] - 2 * d[1] * d[10] * d[30] -
        2 * d[2] * d[11] * d[30] + 6 * d[3] * d[12] * d[30] + 2 * d[4] * d[13] * d[30] + 2 * d[5] * d[14] * d[30] +
        2 * d[6] * d[15] * d[30] - 2 * d[7] * d[16] * d[30] - 2 * d[8] * d[17] * d[30] + 2 * d[1] * d[9] * d[31] +
        2 * d[0] * d[10] * d[31] + 2 * d[4] * d[12] * d[31] + 2 * d[3] * d[13] * d[31] + 2 * d[7] * d[15] * d[31] +
        2 * d[6] * d[16] * d[31] + 2 * d[2] * d[9] * d[32] + 2 * d[0] * d[11] * d[32] + 2 * d[5] * d[12] * d[32] +
        2 * d[3] * d[14] * d[32] + 2 * d[8] * d[15] * d[32] + 2 * d[6] * d[17] * d[32] + 2 * d[6] * d[12] * d[33] +
        2 * d[7] * d[13] * d[33] + 2 * d[8] * d[14] * d[33] + 2 * d[3] * d[15] * d[33] + 2 * d[4] * d[16] * d[33] +
        2 * d[5] * d[17] * d[33] - 2 * d[7] * d[12] * d[34] + 2 * d[6] * d[13] * d[34] + 2 * d[4] * d[15] * d[34] -
        2 * d[3] * d[16] * d[34] - 2 * d[8] * d[12] * d[35] + 2 * d[6] * d[14] * d[35] + 2 * d[5] * d[15] * d[35] -
        2 * d[3] * d[17] * d[35];
    coeffs[187] =
        2 * d[9] * d[12] * d[27] + 2 * d[10] * d[13] * d[27] + 2 * d[11] * d[14] * d[27] - 2 * d[10] * d[12] * d[28] +
        2 * d[9] * d[13] * d[28] - 2 * d[11] * d[12] * d[29] + 2 * d[9] * d[14] * d[29] + std::pow(d[9], 2) * d[30] -
        std::pow(d[10], 2) * d[30] - std::pow(d[11], 2) * d[30] + 3 * std::pow(d[12], 2) * d[30] +
        std::pow(d[13], 2) * d[30] + std::pow(d[14], 2) * d[30] + std::pow(d[15], 2) * d[30] -
        std::pow(d[16], 2) * d[30] - std::pow(d[17], 2) * d[30] + 2 * d[9] * d[10] * d[31] + 2 * d[12] * d[13] * d[31] +
        2 * d[15] * d[16] * d[31] + 2 * d[9] * d[11] * d[32] + 2 * d[12] * d[14] * d[32] + 2 * d[15] * d[17] * d[32] +
        2 * d[12] * d[15] * d[33] + 2 * d[13] * d[16] * d[33] + 2 * d[14] * d[17] * d[33] + 2 * d[13] * d[15] * d[34] -
        2 * d[12] * d[16] * d[34] + 2 * d[14] * d[15] * d[35] - 2 * d[12] * d[17] * d[35];
    coeffs[188] =
        2 * d[3] * d[18] * d[27] + 2 * d[4] * d[19] * d[27] + 2 * d[5] * d[20] * d[27] + 2 * d[0] * d[21] * d[27] +
        2 * d[1] * d[22] * d[27] + 2 * d[2] * d[23] * d[27] + 2 * d[4] * d[18] * d[28] - 2 * d[3] * d[19] * d[28] -
        2 * d[1] * d[21] * d[28] + 2 * d[0] * d[22] * d[28] + 2 * d[5] * d[18] * d[29] - 2 * d[3] * d[20] * d[29] -
        2 * d[2] * d[21] * d[29] + 2 * d[0] * d[23] * d[29] + 2 * d[0] * d[18] * d[30] - 2 * d[1] * d[19] * d[30] -
        2 * d[2] * d[20] * d[30] + 6 * d[3] * d[21] * d[30] + 2 * d[4] * d[22] * d[30] + 2 * d[5] * d[23] * d[30] +
        2 * d[6] * d[24] * d[30] - 2 * d[7] * d[25] * d[30] - 2 * d[8] * d[26] * d[30] + 2 * d[1] * d[18] * d[31] +
        2 * d[0] * d[19] * d[31] + 2 * d[4] * d[21] * d[31] + 2 * d[3] * d[22] * d[31] + 2 * d[7] * d[24] * d[31] +
        2 * d[6] * d[25] * d[31] + 2 * d[2] * d[18] * d[32] + 2 * d[0] * d[20] * d[32] + 2 * d[5] * d[21] * d[32] +
        2 * d[3] * d[23] * d[32] + 2 * d[8] * d[24] * d[32] + 2 * d[6] * d[26] * d[32] + 2 * d[6] * d[21] * d[33] +
        2 * d[7] * d[22] * d[33] + 2 * d[8] * d[23] * d[33] + 2 * d[3] * d[24] * d[33] + 2 * d[4] * d[25] * d[33] +
        2 * d[5] * d[26] * d[33] - 2 * d[7] * d[21] * d[34] + 2 * d[6] * d[22] * d[34] + 2 * d[4] * d[24] * d[34] -
        2 * d[3] * d[25] * d[34] - 2 * d[8] * d[21] * d[35] + 2 * d[6] * d[23] * d[35] + 2 * d[5] * d[24] * d[35] -
        2 * d[3] * d[26] * d[35];
    coeffs[189] =
        2 * d[12] * d[18] * d[27] + 2 * d[13] * d[19] * d[27] + 2 * d[14] * d[20] * d[27] + 2 * d[9] * d[21] * d[27] +
        2 * d[10] * d[22] * d[27] + 2 * d[11] * d[23] * d[27] + 2 * d[13] * d[18] * d[28] - 2 * d[12] * d[19] * d[28] -
        2 * d[10] * d[21] * d[28] + 2 * d[9] * d[22] * d[28] + 2 * d[14] * d[18] * d[29] - 2 * d[12] * d[20] * d[29] -
        2 * d[11] * d[21] * d[29] + 2 * d[9] * d[23] * d[29] + 2 * d[9] * d[18] * d[30] - 2 * d[10] * d[19] * d[30] -
        2 * d[11] * d[20] * d[30] + 6 * d[12] * d[21] * d[30] + 2 * d[13] * d[22] * d[30] + 2 * d[14] * d[23] * d[30] +
        2 * d[15] * d[24] * d[30] - 2 * d[16] * d[25] * d[30] - 2 * d[17] * d[26] * d[30] + 2 * d[10] * d[18] * d[31] +
        2 * d[9] * d[19] * d[31] + 2 * d[13] * d[21] * d[31] + 2 * d[12] * d[22] * d[31] + 2 * d[16] * d[24] * d[31] +
        2 * d[15] * d[25] * d[31] + 2 * d[11] * d[18] * d[32] + 2 * d[9] * d[20] * d[32] + 2 * d[14] * d[21] * d[32] +
        2 * d[12] * d[23] * d[32] + 2 * d[17] * d[24] * d[32] + 2 * d[15] * d[26] * d[32] + 2 * d[15] * d[21] * d[33] +
        2 * d[16] * d[22] * d[33] + 2 * d[17] * d[23] * d[33] + 2 * d[12] * d[24] * d[33] + 2 * d[13] * d[25] * d[33] +
        2 * d[14] * d[26] * d[33] - 2 * d[16] * d[21] * d[34] + 2 * d[15] * d[22] * d[34] + 2 * d[13] * d[24] * d[34] -
        2 * d[12] * d[25] * d[34] - 2 * d[17] * d[21] * d[35] + 2 * d[15] * d[23] * d[35] + 2 * d[14] * d[24] * d[35] -
        2 * d[12] * d[26] * d[35];
    coeffs[190] =
        2 * d[18] * d[21] * d[27] + 2 * d[19] * d[22] * d[27] + 2 * d[20] * d[23] * d[27] - 2 * d[19] * d[21] * d[28] +
        2 * d[18] * d[22] * d[28] - 2 * d[20] * d[21] * d[29] + 2 * d[18] * d[23] * d[29] + std::pow(d[18], 2) * d[30] -
        std::pow(d[19], 2) * d[30] - std::pow(d[20], 2) * d[30] + 3 * std::pow(d[21], 2) * d[30] +
        std::pow(d[22], 2) * d[30] + std::pow(d[23], 2) * d[30] + std::pow(d[24], 2) * d[30] -
        std::pow(d[25], 2) * d[30] - std::pow(d[26], 2) * d[30] + 2 * d[18] * d[19] * d[31] +
        2 * d[21] * d[22] * d[31] + 2 * d[24] * d[25] * d[31] + 2 * d[18] * d[20] * d[32] + 2 * d[21] * d[23] * d[32] +
        2 * d[24] * d[26] * d[32] + 2 * d[21] * d[24] * d[33] + 2 * d[22] * d[25] * d[33] + 2 * d[23] * d[26] * d[33] +
        2 * d[22] * d[24] * d[34] - 2 * d[21] * d[25] * d[34] + 2 * d[23] * d[24] * d[35] - 2 * d[21] * d[26] * d[35];
    coeffs[191] =
        d[3] * std::pow(d[27], 2) + 2 * d[4] * d[27] * d[28] - d[3] * std::pow(d[28], 2) + 2 * d[5] * d[27] * d[29] -
        d[3] * std::pow(d[29], 2) + 2 * d[0] * d[27] * d[30] - 2 * d[1] * d[28] * d[30] - 2 * d[2] * d[29] * d[30] +
        3 * d[3] * std::pow(d[30], 2) + 2 * d[1] * d[27] * d[31] + 2 * d[0] * d[28] * d[31] + 2 * d[4] * d[30] * d[31] +
        d[3] * std::pow(d[31], 2) + 2 * d[2] * d[27] * d[32] + 2 * d[0] * d[29] * d[32] + 2 * d[5] * d[30] * d[32] +
        d[3] * std::pow(d[32], 2) + 2 * d[6] * d[30] * d[33] + 2 * d[7] * d[31] * d[33] + 2 * d[8] * d[32] * d[33] +
        d[3] * std::pow(d[33], 2) - 2 * d[7] * d[30] * d[34] + 2 * d[6] * d[31] * d[34] + 2 * d[4] * d[33] * d[34] -
        d[3] * std::pow(d[34], 2) - 2 * d[8] * d[30] * d[35] + 2 * d[6] * d[32] * d[35] + 2 * d[5] * d[33] * d[35] -
        d[3] * std::pow(d[35], 2);
    coeffs[192] =
        d[12] * std::pow(d[27], 2) + 2 * d[13] * d[27] * d[28] - d[12] * std::pow(d[28], 2) +
        2 * d[14] * d[27] * d[29] - d[12] * std::pow(d[29], 2) + 2 * d[9] * d[27] * d[30] - 2 * d[10] * d[28] * d[30] -
        2 * d[11] * d[29] * d[30] + 3 * d[12] * std::pow(d[30], 2) + 2 * d[10] * d[27] * d[31] +
        2 * d[9] * d[28] * d[31] + 2 * d[13] * d[30] * d[31] + d[12] * std::pow(d[31], 2) + 2 * d[11] * d[27] * d[32] +
        2 * d[9] * d[29] * d[32] + 2 * d[14] * d[30] * d[32] + d[12] * std::pow(d[32], 2) + 2 * d[15] * d[30] * d[33] +
        2 * d[16] * d[31] * d[33] + 2 * d[17] * d[32] * d[33] + d[12] * std::pow(d[33], 2) - 2 * d[16] * d[30] * d[34] +
        2 * d[15] * d[31] * d[34] + 2 * d[13] * d[33] * d[34] - d[12] * std::pow(d[34], 2) - 2 * d[17] * d[30] * d[35] +
        2 * d[15] * d[32] * d[35] + 2 * d[14] * d[33] * d[35] - d[12] * std::pow(d[35], 2);
    coeffs[193] =
        d[21] * std::pow(d[27], 2) + 2 * d[22] * d[27] * d[28] - d[21] * std::pow(d[28], 2) +
        2 * d[23] * d[27] * d[29] - d[21] * std::pow(d[29], 2) + 2 * d[18] * d[27] * d[30] - 2 * d[19] * d[28] * d[30] -
        2 * d[20] * d[29] * d[30] + 3 * d[21] * std::pow(d[30], 2) + 2 * d[19] * d[27] * d[31] +
        2 * d[18] * d[28] * d[31] + 2 * d[22] * d[30] * d[31] + d[21] * std::pow(d[31], 2) + 2 * d[20] * d[27] * d[32] +
        2 * d[18] * d[29] * d[32] + 2 * d[23] * d[30] * d[32] + d[21] * std::pow(d[32], 2) + 2 * d[24] * d[30] * d[33] +
        2 * d[25] * d[31] * d[33] + 2 * d[26] * d[32] * d[33] + d[21] * std::pow(d[33], 2) - 2 * d[25] * d[30] * d[34] +
        2 * d[24] * d[31] * d[34] + 2 * d[22] * d[33] * d[34] - d[21] * std::pow(d[34], 2) - 2 * d[26] * d[30] * d[35] +
        2 * d[24] * d[32] * d[35] + 2 * d[23] * d[33] * d[35] - d[21] * std::pow(d[35], 2);
    coeffs[194] = std::pow(d[27], 2) * d[30] - std::pow(d[28], 2) * d[30] - std::pow(d[29], 2) * d[30] +
                  std::pow(d[30], 3) + 2 * d[27] * d[28] * d[31] + d[30] * std::pow(d[31], 2) +
                  2 * d[27] * d[29] * d[32] + d[30] * std::pow(d[32], 2) + d[30] * std::pow(d[33], 2) +
                  2 * d[31] * d[33] * d[34] - d[30] * std::pow(d[34], 2) + 2 * d[32] * d[33] * d[35] -
                  d[30] * std::pow(d[35], 2);
    coeffs[195] =
        2 * d[0] * d[3] * d[36] + 2 * d[1] * d[4] * d[36] + 2 * d[2] * d[5] * d[36] - 2 * d[1] * d[3] * d[37] +
        2 * d[0] * d[4] * d[37] - 2 * d[2] * d[3] * d[38] + 2 * d[0] * d[5] * d[38] + std::pow(d[0], 2) * d[39] -
        std::pow(d[1], 2) * d[39] - std::pow(d[2], 2) * d[39] + 3 * std::pow(d[3], 2) * d[39] +
        std::pow(d[4], 2) * d[39] + std::pow(d[5], 2) * d[39] + std::pow(d[6], 2) * d[39] - std::pow(d[7], 2) * d[39] -
        std::pow(d[8], 2) * d[39] + 2 * d[0] * d[1] * d[40] + 2 * d[3] * d[4] * d[40] + 2 * d[6] * d[7] * d[40] +
        2 * d[0] * d[2] * d[41] + 2 * d[3] * d[5] * d[41] + 2 * d[6] * d[8] * d[41] + 2 * d[3] * d[6] * d[42] +
        2 * d[4] * d[7] * d[42] + 2 * d[5] * d[8] * d[42] + 2 * d[4] * d[6] * d[43] - 2 * d[3] * d[7] * d[43] +
        2 * d[5] * d[6] * d[44] - 2 * d[3] * d[8] * d[44];
    coeffs[196] =
        2 * d[3] * d[9] * d[36] + 2 * d[4] * d[10] * d[36] + 2 * d[5] * d[11] * d[36] + 2 * d[0] * d[12] * d[36] +
        2 * d[1] * d[13] * d[36] + 2 * d[2] * d[14] * d[36] + 2 * d[4] * d[9] * d[37] - 2 * d[3] * d[10] * d[37] -
        2 * d[1] * d[12] * d[37] + 2 * d[0] * d[13] * d[37] + 2 * d[5] * d[9] * d[38] - 2 * d[3] * d[11] * d[38] -
        2 * d[2] * d[12] * d[38] + 2 * d[0] * d[14] * d[38] + 2 * d[0] * d[9] * d[39] - 2 * d[1] * d[10] * d[39] -
        2 * d[2] * d[11] * d[39] + 6 * d[3] * d[12] * d[39] + 2 * d[4] * d[13] * d[39] + 2 * d[5] * d[14] * d[39] +
        2 * d[6] * d[15] * d[39] - 2 * d[7] * d[16] * d[39] - 2 * d[8] * d[17] * d[39] + 2 * d[1] * d[9] * d[40] +
        2 * d[0] * d[10] * d[40] + 2 * d[4] * d[12] * d[40] + 2 * d[3] * d[13] * d[40] + 2 * d[7] * d[15] * d[40] +
        2 * d[6] * d[16] * d[40] + 2 * d[2] * d[9] * d[41] + 2 * d[0] * d[11] * d[41] + 2 * d[5] * d[12] * d[41] +
        2 * d[3] * d[14] * d[41] + 2 * d[8] * d[15] * d[41] + 2 * d[6] * d[17] * d[41] + 2 * d[6] * d[12] * d[42] +
        2 * d[7] * d[13] * d[42] + 2 * d[8] * d[14] * d[42] + 2 * d[3] * d[15] * d[42] + 2 * d[4] * d[16] * d[42] +
        2 * d[5] * d[17] * d[42] - 2 * d[7] * d[12] * d[43] + 2 * d[6] * d[13] * d[43] + 2 * d[4] * d[15] * d[43] -
        2 * d[3] * d[16] * d[43] - 2 * d[8] * d[12] * d[44] + 2 * d[6] * d[14] * d[44] + 2 * d[5] * d[15] * d[44] -
        2 * d[3] * d[17] * d[44];
    coeffs[197] =
        2 * d[9] * d[12] * d[36] + 2 * d[10] * d[13] * d[36] + 2 * d[11] * d[14] * d[36] - 2 * d[10] * d[12] * d[37] +
        2 * d[9] * d[13] * d[37] - 2 * d[11] * d[12] * d[38] + 2 * d[9] * d[14] * d[38] + std::pow(d[9], 2) * d[39] -
        std::pow(d[10], 2) * d[39] - std::pow(d[11], 2) * d[39] + 3 * std::pow(d[12], 2) * d[39] +
        std::pow(d[13], 2) * d[39] + std::pow(d[14], 2) * d[39] + std::pow(d[15], 2) * d[39] -
        std::pow(d[16], 2) * d[39] - std::pow(d[17], 2) * d[39] + 2 * d[9] * d[10] * d[40] + 2 * d[12] * d[13] * d[40] +
        2 * d[15] * d[16] * d[40] + 2 * d[9] * d[11] * d[41] + 2 * d[12] * d[14] * d[41] + 2 * d[15] * d[17] * d[41] +
        2 * d[12] * d[15] * d[42] + 2 * d[13] * d[16] * d[42] + 2 * d[14] * d[17] * d[42] + 2 * d[13] * d[15] * d[43] -
        2 * d[12] * d[16] * d[43] + 2 * d[14] * d[15] * d[44] - 2 * d[12] * d[17] * d[44];
    coeffs[198] =
        2 * d[3] * d[18] * d[36] + 2 * d[4] * d[19] * d[36] + 2 * d[5] * d[20] * d[36] + 2 * d[0] * d[21] * d[36] +
        2 * d[1] * d[22] * d[36] + 2 * d[2] * d[23] * d[36] + 2 * d[4] * d[18] * d[37] - 2 * d[3] * d[19] * d[37] -
        2 * d[1] * d[21] * d[37] + 2 * d[0] * d[22] * d[37] + 2 * d[5] * d[18] * d[38] - 2 * d[3] * d[20] * d[38] -
        2 * d[2] * d[21] * d[38] + 2 * d[0] * d[23] * d[38] + 2 * d[0] * d[18] * d[39] - 2 * d[1] * d[19] * d[39] -
        2 * d[2] * d[20] * d[39] + 6 * d[3] * d[21] * d[39] + 2 * d[4] * d[22] * d[39] + 2 * d[5] * d[23] * d[39] +
        2 * d[6] * d[24] * d[39] - 2 * d[7] * d[25] * d[39] - 2 * d[8] * d[26] * d[39] + 2 * d[1] * d[18] * d[40] +
        2 * d[0] * d[19] * d[40] + 2 * d[4] * d[21] * d[40] + 2 * d[3] * d[22] * d[40] + 2 * d[7] * d[24] * d[40] +
        2 * d[6] * d[25] * d[40] + 2 * d[2] * d[18] * d[41] + 2 * d[0] * d[20] * d[41] + 2 * d[5] * d[21] * d[41] +
        2 * d[3] * d[23] * d[41] + 2 * d[8] * d[24] * d[41] + 2 * d[6] * d[26] * d[41] + 2 * d[6] * d[21] * d[42] +
        2 * d[7] * d[22] * d[42] + 2 * d[8] * d[23] * d[42] + 2 * d[3] * d[24] * d[42] + 2 * d[4] * d[25] * d[42] +
        2 * d[5] * d[26] * d[42] - 2 * d[7] * d[21] * d[43] + 2 * d[6] * d[22] * d[43] + 2 * d[4] * d[24] * d[43] -
        2 * d[3] * d[25] * d[43] - 2 * d[8] * d[21] * d[44] + 2 * d[6] * d[23] * d[44] + 2 * d[5] * d[24] * d[44] -
        2 * d[3] * d[26] * d[44];
    coeffs[199] =
        2 * d[12] * d[18] * d[36] + 2 * d[13] * d[19] * d[36] + 2 * d[14] * d[20] * d[36] + 2 * d[9] * d[21] * d[36] +
        2 * d[10] * d[22] * d[36] + 2 * d[11] * d[23] * d[36] + 2 * d[13] * d[18] * d[37] - 2 * d[12] * d[19] * d[37] -
        2 * d[10] * d[21] * d[37] + 2 * d[9] * d[22] * d[37] + 2 * d[14] * d[18] * d[38] - 2 * d[12] * d[20] * d[38] -
        2 * d[11] * d[21] * d[38] + 2 * d[9] * d[23] * d[38] + 2 * d[9] * d[18] * d[39] - 2 * d[10] * d[19] * d[39] -
        2 * d[11] * d[20] * d[39] + 6 * d[12] * d[21] * d[39] + 2 * d[13] * d[22] * d[39] + 2 * d[14] * d[23] * d[39] +
        2 * d[15] * d[24] * d[39] - 2 * d[16] * d[25] * d[39] - 2 * d[17] * d[26] * d[39] + 2 * d[10] * d[18] * d[40] +
        2 * d[9] * d[19] * d[40] + 2 * d[13] * d[21] * d[40] + 2 * d[12] * d[22] * d[40] + 2 * d[16] * d[24] * d[40] +
        2 * d[15] * d[25] * d[40] + 2 * d[11] * d[18] * d[41] + 2 * d[9] * d[20] * d[41] + 2 * d[14] * d[21] * d[41] +
        2 * d[12] * d[23] * d[41] + 2 * d[17] * d[24] * d[41] + 2 * d[15] * d[26] * d[41] + 2 * d[15] * d[21] * d[42] +
        2 * d[16] * d[22] * d[42] + 2 * d[17] * d[23] * d[42] + 2 * d[12] * d[24] * d[42] + 2 * d[13] * d[25] * d[42] +
        2 * d[14] * d[26] * d[42] - 2 * d[16] * d[21] * d[43] + 2 * d[15] * d[22] * d[43] + 2 * d[13] * d[24] * d[43] -
        2 * d[12] * d[25] * d[43] - 2 * d[17] * d[21] * d[44] + 2 * d[15] * d[23] * d[44] + 2 * d[14] * d[24] * d[44] -
        2 * d[12] * d[26] * d[44];
    coeffs[200] =
        2 * d[18] * d[21] * d[36] + 2 * d[19] * d[22] * d[36] + 2 * d[20] * d[23] * d[36] - 2 * d[19] * d[21] * d[37] +
        2 * d[18] * d[22] * d[37] - 2 * d[20] * d[21] * d[38] + 2 * d[18] * d[23] * d[38] + std::pow(d[18], 2) * d[39] -
        std::pow(d[19], 2) * d[39] - std::pow(d[20], 2) * d[39] + 3 * std::pow(d[21], 2) * d[39] +
        std::pow(d[22], 2) * d[39] + std::pow(d[23], 2) * d[39] + std::pow(d[24], 2) * d[39] -
        std::pow(d[25], 2) * d[39] - std::pow(d[26], 2) * d[39] + 2 * d[18] * d[19] * d[40] +
        2 * d[21] * d[22] * d[40] + 2 * d[24] * d[25] * d[40] + 2 * d[18] * d[20] * d[41] + 2 * d[21] * d[23] * d[41] +
        2 * d[24] * d[26] * d[41] + 2 * d[21] * d[24] * d[42] + 2 * d[22] * d[25] * d[42] + 2 * d[23] * d[26] * d[42] +
        2 * d[22] * d[24] * d[43] - 2 * d[21] * d[25] * d[43] + 2 * d[23] * d[24] * d[44] - 2 * d[21] * d[26] * d[44];
    coeffs[201] =
        2 * d[3] * d[27] * d[36] + 2 * d[4] * d[28] * d[36] + 2 * d[5] * d[29] * d[36] + 2 * d[0] * d[30] * d[36] +
        2 * d[1] * d[31] * d[36] + 2 * d[2] * d[32] * d[36] + 2 * d[4] * d[27] * d[37] - 2 * d[3] * d[28] * d[37] -
        2 * d[1] * d[30] * d[37] + 2 * d[0] * d[31] * d[37] + 2 * d[5] * d[27] * d[38] - 2 * d[3] * d[29] * d[38] -
        2 * d[2] * d[30] * d[38] + 2 * d[0] * d[32] * d[38] + 2 * d[0] * d[27] * d[39] - 2 * d[1] * d[28] * d[39] -
        2 * d[2] * d[29] * d[39] + 6 * d[3] * d[30] * d[39] + 2 * d[4] * d[31] * d[39] + 2 * d[5] * d[32] * d[39] +
        2 * d[6] * d[33] * d[39] - 2 * d[7] * d[34] * d[39] - 2 * d[8] * d[35] * d[39] + 2 * d[1] * d[27] * d[40] +
        2 * d[0] * d[28] * d[40] + 2 * d[4] * d[30] * d[40] + 2 * d[3] * d[31] * d[40] + 2 * d[7] * d[33] * d[40] +
        2 * d[6] * d[34] * d[40] + 2 * d[2] * d[27] * d[41] + 2 * d[0] * d[29] * d[41] + 2 * d[5] * d[30] * d[41] +
        2 * d[3] * d[32] * d[41] + 2 * d[8] * d[33] * d[41] + 2 * d[6] * d[35] * d[41] + 2 * d[6] * d[30] * d[42] +
        2 * d[7] * d[31] * d[42] + 2 * d[8] * d[32] * d[42] + 2 * d[3] * d[33] * d[42] + 2 * d[4] * d[34] * d[42] +
        2 * d[5] * d[35] * d[42] - 2 * d[7] * d[30] * d[43] + 2 * d[6] * d[31] * d[43] + 2 * d[4] * d[33] * d[43] -
        2 * d[3] * d[34] * d[43] - 2 * d[8] * d[30] * d[44] + 2 * d[6] * d[32] * d[44] + 2 * d[5] * d[33] * d[44] -
        2 * d[3] * d[35] * d[44];
    coeffs[202] =
        2 * d[12] * d[27] * d[36] + 2 * d[13] * d[28] * d[36] + 2 * d[14] * d[29] * d[36] + 2 * d[9] * d[30] * d[36] +
        2 * d[10] * d[31] * d[36] + 2 * d[11] * d[32] * d[36] + 2 * d[13] * d[27] * d[37] - 2 * d[12] * d[28] * d[37] -
        2 * d[10] * d[30] * d[37] + 2 * d[9] * d[31] * d[37] + 2 * d[14] * d[27] * d[38] - 2 * d[12] * d[29] * d[38] -
        2 * d[11] * d[30] * d[38] + 2 * d[9] * d[32] * d[38] + 2 * d[9] * d[27] * d[39] - 2 * d[10] * d[28] * d[39] -
        2 * d[11] * d[29] * d[39] + 6 * d[12] * d[30] * d[39] + 2 * d[13] * d[31] * d[39] + 2 * d[14] * d[32] * d[39] +
        2 * d[15] * d[33] * d[39] - 2 * d[16] * d[34] * d[39] - 2 * d[17] * d[35] * d[39] + 2 * d[10] * d[27] * d[40] +
        2 * d[9] * d[28] * d[40] + 2 * d[13] * d[30] * d[40] + 2 * d[12] * d[31] * d[40] + 2 * d[16] * d[33] * d[40] +
        2 * d[15] * d[34] * d[40] + 2 * d[11] * d[27] * d[41] + 2 * d[9] * d[29] * d[41] + 2 * d[14] * d[30] * d[41] +
        2 * d[12] * d[32] * d[41] + 2 * d[17] * d[33] * d[41] + 2 * d[15] * d[35] * d[41] + 2 * d[15] * d[30] * d[42] +
        2 * d[16] * d[31] * d[42] + 2 * d[17] * d[32] * d[42] + 2 * d[12] * d[33] * d[42] + 2 * d[13] * d[34] * d[42] +
        2 * d[14] * d[35] * d[42] - 2 * d[16] * d[30] * d[43] + 2 * d[15] * d[31] * d[43] + 2 * d[13] * d[33] * d[43] -
        2 * d[12] * d[34] * d[43] - 2 * d[17] * d[30] * d[44] + 2 * d[15] * d[32] * d[44] + 2 * d[14] * d[33] * d[44] -
        2 * d[12] * d[35] * d[44];
    coeffs[203] =
        2 * d[21] * d[27] * d[36] + 2 * d[22] * d[28] * d[36] + 2 * d[23] * d[29] * d[36] + 2 * d[18] * d[30] * d[36] +
        2 * d[19] * d[31] * d[36] + 2 * d[20] * d[32] * d[36] + 2 * d[22] * d[27] * d[37] - 2 * d[21] * d[28] * d[37] -
        2 * d[19] * d[30] * d[37] + 2 * d[18] * d[31] * d[37] + 2 * d[23] * d[27] * d[38] - 2 * d[21] * d[29] * d[38] -
        2 * d[20] * d[30] * d[38] + 2 * d[18] * d[32] * d[38] + 2 * d[18] * d[27] * d[39] - 2 * d[19] * d[28] * d[39] -
        2 * d[20] * d[29] * d[39] + 6 * d[21] * d[30] * d[39] + 2 * d[22] * d[31] * d[39] + 2 * d[23] * d[32] * d[39] +
        2 * d[24] * d[33] * d[39] - 2 * d[25] * d[34] * d[39] - 2 * d[26] * d[35] * d[39] + 2 * d[19] * d[27] * d[40] +
        2 * d[18] * d[28] * d[40] + 2 * d[22] * d[30] * d[40] + 2 * d[21] * d[31] * d[40] + 2 * d[25] * d[33] * d[40] +
        2 * d[24] * d[34] * d[40] + 2 * d[20] * d[27] * d[41] + 2 * d[18] * d[29] * d[41] + 2 * d[23] * d[30] * d[41] +
        2 * d[21] * d[32] * d[41] + 2 * d[26] * d[33] * d[41] + 2 * d[24] * d[35] * d[41] + 2 * d[24] * d[30] * d[42] +
        2 * d[25] * d[31] * d[42] + 2 * d[26] * d[32] * d[42] + 2 * d[21] * d[33] * d[42] + 2 * d[22] * d[34] * d[42] +
        2 * d[23] * d[35] * d[42] - 2 * d[25] * d[30] * d[43] + 2 * d[24] * d[31] * d[43] + 2 * d[22] * d[33] * d[43] -
        2 * d[21] * d[34] * d[43] - 2 * d[26] * d[30] * d[44] + 2 * d[24] * d[32] * d[44] + 2 * d[23] * d[33] * d[44] -
        2 * d[21] * d[35] * d[44];
    coeffs[204] =
        2 * d[27] * d[30] * d[36] + 2 * d[28] * d[31] * d[36] + 2 * d[29] * d[32] * d[36] - 2 * d[28] * d[30] * d[37] +
        2 * d[27] * d[31] * d[37] - 2 * d[29] * d[30] * d[38] + 2 * d[27] * d[32] * d[38] + std::pow(d[27], 2) * d[39] -
        std::pow(d[28], 2) * d[39] - std::pow(d[29], 2) * d[39] + 3 * std::pow(d[30], 2) * d[39] +
        std::pow(d[31], 2) * d[39] + std::pow(d[32], 2) * d[39] + std::pow(d[33], 2) * d[39] -
        std::pow(d[34], 2) * d[39] - std::pow(d[35], 2) * d[39] + 2 * d[27] * d[28] * d[40] +
        2 * d[30] * d[31] * d[40] + 2 * d[33] * d[34] * d[40] + 2 * d[27] * d[29] * d[41] + 2 * d[30] * d[32] * d[41] +
        2 * d[33] * d[35] * d[41] + 2 * d[30] * d[33] * d[42] + 2 * d[31] * d[34] * d[42] + 2 * d[32] * d[35] * d[42] +
        2 * d[31] * d[33] * d[43] - 2 * d[30] * d[34] * d[43] + 2 * d[32] * d[33] * d[44] - 2 * d[30] * d[35] * d[44];
    coeffs[205] =
        d[3] * std::pow(d[36], 2) + 2 * d[4] * d[36] * d[37] - d[3] * std::pow(d[37], 2) + 2 * d[5] * d[36] * d[38] -
        d[3] * std::pow(d[38], 2) + 2 * d[0] * d[36] * d[39] - 2 * d[1] * d[37] * d[39] - 2 * d[2] * d[38] * d[39] +
        3 * d[3] * std::pow(d[39], 2) + 2 * d[1] * d[36] * d[40] + 2 * d[0] * d[37] * d[40] + 2 * d[4] * d[39] * d[40] +
        d[3] * std::pow(d[40], 2) + 2 * d[2] * d[36] * d[41] + 2 * d[0] * d[38] * d[41] + 2 * d[5] * d[39] * d[41] +
        d[3] * std::pow(d[41], 2) + 2 * d[6] * d[39] * d[42] + 2 * d[7] * d[40] * d[42] + 2 * d[8] * d[41] * d[42] +
        d[3] * std::pow(d[42], 2) - 2 * d[7] * d[39] * d[43] + 2 * d[6] * d[40] * d[43] + 2 * d[4] * d[42] * d[43] -
        d[3] * std::pow(d[43], 2) - 2 * d[8] * d[39] * d[44] + 2 * d[6] * d[41] * d[44] + 2 * d[5] * d[42] * d[44] -
        d[3] * std::pow(d[44], 2);
    coeffs[206] =
        d[12] * std::pow(d[36], 2) + 2 * d[13] * d[36] * d[37] - d[12] * std::pow(d[37], 2) +
        2 * d[14] * d[36] * d[38] - d[12] * std::pow(d[38], 2) + 2 * d[9] * d[36] * d[39] - 2 * d[10] * d[37] * d[39] -
        2 * d[11] * d[38] * d[39] + 3 * d[12] * std::pow(d[39], 2) + 2 * d[10] * d[36] * d[40] +
        2 * d[9] * d[37] * d[40] + 2 * d[13] * d[39] * d[40] + d[12] * std::pow(d[40], 2) + 2 * d[11] * d[36] * d[41] +
        2 * d[9] * d[38] * d[41] + 2 * d[14] * d[39] * d[41] + d[12] * std::pow(d[41], 2) + 2 * d[15] * d[39] * d[42] +
        2 * d[16] * d[40] * d[42] + 2 * d[17] * d[41] * d[42] + d[12] * std::pow(d[42], 2) - 2 * d[16] * d[39] * d[43] +
        2 * d[15] * d[40] * d[43] + 2 * d[13] * d[42] * d[43] - d[12] * std::pow(d[43], 2) - 2 * d[17] * d[39] * d[44] +
        2 * d[15] * d[41] * d[44] + 2 * d[14] * d[42] * d[44] - d[12] * std::pow(d[44], 2);
    coeffs[207] =
        d[21] * std::pow(d[36], 2) + 2 * d[22] * d[36] * d[37] - d[21] * std::pow(d[37], 2) +
        2 * d[23] * d[36] * d[38] - d[21] * std::pow(d[38], 2) + 2 * d[18] * d[36] * d[39] - 2 * d[19] * d[37] * d[39] -
        2 * d[20] * d[38] * d[39] + 3 * d[21] * std::pow(d[39], 2) + 2 * d[19] * d[36] * d[40] +
        2 * d[18] * d[37] * d[40] + 2 * d[22] * d[39] * d[40] + d[21] * std::pow(d[40], 2) + 2 * d[20] * d[36] * d[41] +
        2 * d[18] * d[38] * d[41] + 2 * d[23] * d[39] * d[41] + d[21] * std::pow(d[41], 2) + 2 * d[24] * d[39] * d[42] +
        2 * d[25] * d[40] * d[42] + 2 * d[26] * d[41] * d[42] + d[21] * std::pow(d[42], 2) - 2 * d[25] * d[39] * d[43] +
        2 * d[24] * d[40] * d[43] + 2 * d[22] * d[42] * d[43] - d[21] * std::pow(d[43], 2) - 2 * d[26] * d[39] * d[44] +
        2 * d[24] * d[41] * d[44] + 2 * d[23] * d[42] * d[44] - d[21] * std::pow(d[44], 2);
    coeffs[208] =
        d[30] * std::pow(d[36], 2) + 2 * d[31] * d[36] * d[37] - d[30] * std::pow(d[37], 2) +
        2 * d[32] * d[36] * d[38] - d[30] * std::pow(d[38], 2) + 2 * d[27] * d[36] * d[39] - 2 * d[28] * d[37] * d[39] -
        2 * d[29] * d[38] * d[39] + 3 * d[30] * std::pow(d[39], 2) + 2 * d[28] * d[36] * d[40] +
        2 * d[27] * d[37] * d[40] + 2 * d[31] * d[39] * d[40] + d[30] * std::pow(d[40], 2) + 2 * d[29] * d[36] * d[41] +
        2 * d[27] * d[38] * d[41] + 2 * d[32] * d[39] * d[41] + d[30] * std::pow(d[41], 2) + 2 * d[33] * d[39] * d[42] +
        2 * d[34] * d[40] * d[42] + 2 * d[35] * d[41] * d[42] + d[30] * std::pow(d[42], 2) - 2 * d[34] * d[39] * d[43] +
        2 * d[33] * d[40] * d[43] + 2 * d[31] * d[42] * d[43] - d[30] * std::pow(d[43], 2) - 2 * d[35] * d[39] * d[44] +
        2 * d[33] * d[41] * d[44] + 2 * d[32] * d[42] * d[44] - d[30] * std::pow(d[44], 2);
    coeffs[209] = std::pow(d[36], 2) * d[39] - std::pow(d[37], 2) * d[39] - std::pow(d[38], 2) * d[39] +
                  std::pow(d[39], 3) + 2 * d[36] * d[37] * d[40] + d[39] * std::pow(d[40], 2) +
                  2 * d[36] * d[38] * d[41] + d[39] * std::pow(d[41], 2) + d[39] * std::pow(d[42], 2) +
                  2 * d[40] * d[42] * d[43] - d[39] * std::pow(d[43], 2) + 2 * d[41] * d[42] * d[44] -
                  d[39] * std::pow(d[44], 2);
    coeffs[210] = 2 * d[0] * d[1] * d[3] - std::pow(d[0], 2) * d[4] + std::pow(d[1], 2) * d[4] -
                  std::pow(d[2], 2) * d[4] + std::pow(d[3], 2) * d[4] + std::pow(d[4], 3) + 2 * d[1] * d[2] * d[5] +
                  d[4] * std::pow(d[5], 2) - d[4] * std::pow(d[6], 2) + 2 * d[3] * d[6] * d[7] +
                  d[4] * std::pow(d[7], 2) + 2 * d[5] * d[7] * d[8] - d[4] * std::pow(d[8], 2);
    coeffs[211] = 2 * d[1] * d[3] * d[9] - 2 * d[0] * d[4] * d[9] + 2 * d[0] * d[3] * d[10] + 2 * d[1] * d[4] * d[10] +
                  2 * d[2] * d[5] * d[10] - 2 * d[2] * d[4] * d[11] + 2 * d[1] * d[5] * d[11] +
                  2 * d[0] * d[1] * d[12] + 2 * d[3] * d[4] * d[12] + 2 * d[6] * d[7] * d[12] -
                  std::pow(d[0], 2) * d[13] + std::pow(d[1], 2) * d[13] - std::pow(d[2], 2) * d[13] +
                  std::pow(d[3], 2) * d[13] + 3 * std::pow(d[4], 2) * d[13] + std::pow(d[5], 2) * d[13] -
                  std::pow(d[6], 2) * d[13] + std::pow(d[7], 2) * d[13] - std::pow(d[8], 2) * d[13] +
                  2 * d[1] * d[2] * d[14] + 2 * d[4] * d[5] * d[14] + 2 * d[7] * d[8] * d[14] -
                  2 * d[4] * d[6] * d[15] + 2 * d[3] * d[7] * d[15] + 2 * d[3] * d[6] * d[16] +
                  2 * d[4] * d[7] * d[16] + 2 * d[5] * d[8] * d[16] + 2 * d[5] * d[7] * d[17] - 2 * d[4] * d[8] * d[17];
    coeffs[212] =
        -d[4] * std::pow(d[9], 2) + 2 * d[3] * d[9] * d[10] + d[4] * std::pow(d[10], 2) + 2 * d[5] * d[10] * d[11] -
        d[4] * std::pow(d[11], 2) + 2 * d[1] * d[9] * d[12] + 2 * d[0] * d[10] * d[12] + d[4] * std::pow(d[12], 2) -
        2 * d[0] * d[9] * d[13] + 2 * d[1] * d[10] * d[13] - 2 * d[2] * d[11] * d[13] + 2 * d[3] * d[12] * d[13] +
        3 * d[4] * std::pow(d[13], 2) + 2 * d[2] * d[10] * d[14] + 2 * d[1] * d[11] * d[14] + 2 * d[5] * d[13] * d[14] +
        d[4] * std::pow(d[14], 2) + 2 * d[7] * d[12] * d[15] - 2 * d[6] * d[13] * d[15] - d[4] * std::pow(d[15], 2) +
        2 * d[6] * d[12] * d[16] + 2 * d[7] * d[13] * d[16] + 2 * d[8] * d[14] * d[16] + 2 * d[3] * d[15] * d[16] +
        d[4] * std::pow(d[16], 2) - 2 * d[8] * d[13] * d[17] + 2 * d[7] * d[14] * d[17] + 2 * d[5] * d[16] * d[17] -
        d[4] * std::pow(d[17], 2);
    coeffs[213] = 2 * d[9] * d[10] * d[12] - std::pow(d[9], 2) * d[13] + std::pow(d[10], 2) * d[13] -
                  std::pow(d[11], 2) * d[13] + std::pow(d[12], 2) * d[13] + std::pow(d[13], 3) +
                  2 * d[10] * d[11] * d[14] + d[13] * std::pow(d[14], 2) - d[13] * std::pow(d[15], 2) +
                  2 * d[12] * d[15] * d[16] + d[13] * std::pow(d[16], 2) + 2 * d[14] * d[16] * d[17] -
                  d[13] * std::pow(d[17], 2);
    coeffs[214] =
        2 * d[1] * d[3] * d[18] - 2 * d[0] * d[4] * d[18] + 2 * d[0] * d[3] * d[19] + 2 * d[1] * d[4] * d[19] +
        2 * d[2] * d[5] * d[19] - 2 * d[2] * d[4] * d[20] + 2 * d[1] * d[5] * d[20] + 2 * d[0] * d[1] * d[21] +
        2 * d[3] * d[4] * d[21] + 2 * d[6] * d[7] * d[21] - std::pow(d[0], 2) * d[22] + std::pow(d[1], 2) * d[22] -
        std::pow(d[2], 2) * d[22] + std::pow(d[3], 2) * d[22] + 3 * std::pow(d[4], 2) * d[22] +
        std::pow(d[5], 2) * d[22] - std::pow(d[6], 2) * d[22] + std::pow(d[7], 2) * d[22] - std::pow(d[8], 2) * d[22] +
        2 * d[1] * d[2] * d[23] + 2 * d[4] * d[5] * d[23] + 2 * d[7] * d[8] * d[23] - 2 * d[4] * d[6] * d[24] +
        2 * d[3] * d[7] * d[24] + 2 * d[3] * d[6] * d[25] + 2 * d[4] * d[7] * d[25] + 2 * d[5] * d[8] * d[25] +
        2 * d[5] * d[7] * d[26] - 2 * d[4] * d[8] * d[26];
    coeffs[215] =
        -2 * d[4] * d[9] * d[18] + 2 * d[3] * d[10] * d[18] + 2 * d[1] * d[12] * d[18] - 2 * d[0] * d[13] * d[18] +
        2 * d[3] * d[9] * d[19] + 2 * d[4] * d[10] * d[19] + 2 * d[5] * d[11] * d[19] + 2 * d[0] * d[12] * d[19] +
        2 * d[1] * d[13] * d[19] + 2 * d[2] * d[14] * d[19] + 2 * d[5] * d[10] * d[20] - 2 * d[4] * d[11] * d[20] -
        2 * d[2] * d[13] * d[20] + 2 * d[1] * d[14] * d[20] + 2 * d[1] * d[9] * d[21] + 2 * d[0] * d[10] * d[21] +
        2 * d[4] * d[12] * d[21] + 2 * d[3] * d[13] * d[21] + 2 * d[7] * d[15] * d[21] + 2 * d[6] * d[16] * d[21] -
        2 * d[0] * d[9] * d[22] + 2 * d[1] * d[10] * d[22] - 2 * d[2] * d[11] * d[22] + 2 * d[3] * d[12] * d[22] +
        6 * d[4] * d[13] * d[22] + 2 * d[5] * d[14] * d[22] - 2 * d[6] * d[15] * d[22] + 2 * d[7] * d[16] * d[22] -
        2 * d[8] * d[17] * d[22] + 2 * d[2] * d[10] * d[23] + 2 * d[1] * d[11] * d[23] + 2 * d[5] * d[13] * d[23] +
        2 * d[4] * d[14] * d[23] + 2 * d[8] * d[16] * d[23] + 2 * d[7] * d[17] * d[23] + 2 * d[7] * d[12] * d[24] -
        2 * d[6] * d[13] * d[24] - 2 * d[4] * d[15] * d[24] + 2 * d[3] * d[16] * d[24] + 2 * d[6] * d[12] * d[25] +
        2 * d[7] * d[13] * d[25] + 2 * d[8] * d[14] * d[25] + 2 * d[3] * d[15] * d[25] + 2 * d[4] * d[16] * d[25] +
        2 * d[5] * d[17] * d[25] - 2 * d[8] * d[13] * d[26] + 2 * d[7] * d[14] * d[26] + 2 * d[5] * d[16] * d[26] -
        2 * d[4] * d[17] * d[26];
    coeffs[216] =
        2 * d[10] * d[12] * d[18] - 2 * d[9] * d[13] * d[18] + 2 * d[9] * d[12] * d[19] + 2 * d[10] * d[13] * d[19] +
        2 * d[11] * d[14] * d[19] - 2 * d[11] * d[13] * d[20] + 2 * d[10] * d[14] * d[20] + 2 * d[9] * d[10] * d[21] +
        2 * d[12] * d[13] * d[21] + 2 * d[15] * d[16] * d[21] - std::pow(d[9], 2) * d[22] + std::pow(d[10], 2) * d[22] -
        std::pow(d[11], 2) * d[22] + std::pow(d[12], 2) * d[22] + 3 * std::pow(d[13], 2) * d[22] +
        std::pow(d[14], 2) * d[22] - std::pow(d[15], 2) * d[22] + std::pow(d[16], 2) * d[22] -
        std::pow(d[17], 2) * d[22] + 2 * d[10] * d[11] * d[23] + 2 * d[13] * d[14] * d[23] + 2 * d[16] * d[17] * d[23] -
        2 * d[13] * d[15] * d[24] + 2 * d[12] * d[16] * d[24] + 2 * d[12] * d[15] * d[25] + 2 * d[13] * d[16] * d[25] +
        2 * d[14] * d[17] * d[25] + 2 * d[14] * d[16] * d[26] - 2 * d[13] * d[17] * d[26];
    coeffs[217] =
        -d[4] * std::pow(d[18], 2) + 2 * d[3] * d[18] * d[19] + d[4] * std::pow(d[19], 2) + 2 * d[5] * d[19] * d[20] -
        d[4] * std::pow(d[20], 2) + 2 * d[1] * d[18] * d[21] + 2 * d[0] * d[19] * d[21] + d[4] * std::pow(d[21], 2) -
        2 * d[0] * d[18] * d[22] + 2 * d[1] * d[19] * d[22] - 2 * d[2] * d[20] * d[22] + 2 * d[3] * d[21] * d[22] +
        3 * d[4] * std::pow(d[22], 2) + 2 * d[2] * d[19] * d[23] + 2 * d[1] * d[20] * d[23] + 2 * d[5] * d[22] * d[23] +
        d[4] * std::pow(d[23], 2) + 2 * d[7] * d[21] * d[24] - 2 * d[6] * d[22] * d[24] - d[4] * std::pow(d[24], 2) +
        2 * d[6] * d[21] * d[25] + 2 * d[7] * d[22] * d[25] + 2 * d[8] * d[23] * d[25] + 2 * d[3] * d[24] * d[25] +
        d[4] * std::pow(d[25], 2) - 2 * d[8] * d[22] * d[26] + 2 * d[7] * d[23] * d[26] + 2 * d[5] * d[25] * d[26] -
        d[4] * std::pow(d[26], 2);
    coeffs[218] =
        -d[13] * std::pow(d[18], 2) + 2 * d[12] * d[18] * d[19] + d[13] * std::pow(d[19], 2) +
        2 * d[14] * d[19] * d[20] - d[13] * std::pow(d[20], 2) + 2 * d[10] * d[18] * d[21] + 2 * d[9] * d[19] * d[21] +
        d[13] * std::pow(d[21], 2) - 2 * d[9] * d[18] * d[22] + 2 * d[10] * d[19] * d[22] - 2 * d[11] * d[20] * d[22] +
        2 * d[12] * d[21] * d[22] + 3 * d[13] * std::pow(d[22], 2) + 2 * d[11] * d[19] * d[23] +
        2 * d[10] * d[20] * d[23] + 2 * d[14] * d[22] * d[23] + d[13] * std::pow(d[23], 2) + 2 * d[16] * d[21] * d[24] -
        2 * d[15] * d[22] * d[24] - d[13] * std::pow(d[24], 2) + 2 * d[15] * d[21] * d[25] + 2 * d[16] * d[22] * d[25] +
        2 * d[17] * d[23] * d[25] + 2 * d[12] * d[24] * d[25] + d[13] * std::pow(d[25], 2) - 2 * d[17] * d[22] * d[26] +
        2 * d[16] * d[23] * d[26] + 2 * d[14] * d[25] * d[26] - d[13] * std::pow(d[26], 2);
    coeffs[219] = 2 * d[18] * d[19] * d[21] - std::pow(d[18], 2) * d[22] + std::pow(d[19], 2) * d[22] -
                  std::pow(d[20], 2) * d[22] + std::pow(d[21], 2) * d[22] + std::pow(d[22], 3) +
                  2 * d[19] * d[20] * d[23] + d[22] * std::pow(d[23], 2) - d[22] * std::pow(d[24], 2) +
                  2 * d[21] * d[24] * d[25] + d[22] * std::pow(d[25], 2) + 2 * d[23] * d[25] * d[26] -
                  d[22] * std::pow(d[26], 2);
    coeffs[220] =
        2 * d[1] * d[3] * d[27] - 2 * d[0] * d[4] * d[27] + 2 * d[0] * d[3] * d[28] + 2 * d[1] * d[4] * d[28] +
        2 * d[2] * d[5] * d[28] - 2 * d[2] * d[4] * d[29] + 2 * d[1] * d[5] * d[29] + 2 * d[0] * d[1] * d[30] +
        2 * d[3] * d[4] * d[30] + 2 * d[6] * d[7] * d[30] - std::pow(d[0], 2) * d[31] + std::pow(d[1], 2) * d[31] -
        std::pow(d[2], 2) * d[31] + std::pow(d[3], 2) * d[31] + 3 * std::pow(d[4], 2) * d[31] +
        std::pow(d[5], 2) * d[31] - std::pow(d[6], 2) * d[31] + std::pow(d[7], 2) * d[31] - std::pow(d[8], 2) * d[31] +
        2 * d[1] * d[2] * d[32] + 2 * d[4] * d[5] * d[32] + 2 * d[7] * d[8] * d[32] - 2 * d[4] * d[6] * d[33] +
        2 * d[3] * d[7] * d[33] + 2 * d[3] * d[6] * d[34] + 2 * d[4] * d[7] * d[34] + 2 * d[5] * d[8] * d[34] +
        2 * d[5] * d[7] * d[35] - 2 * d[4] * d[8] * d[35];
    coeffs[221] =
        -2 * d[4] * d[9] * d[27] + 2 * d[3] * d[10] * d[27] + 2 * d[1] * d[12] * d[27] - 2 * d[0] * d[13] * d[27] +
        2 * d[3] * d[9] * d[28] + 2 * d[4] * d[10] * d[28] + 2 * d[5] * d[11] * d[28] + 2 * d[0] * d[12] * d[28] +
        2 * d[1] * d[13] * d[28] + 2 * d[2] * d[14] * d[28] + 2 * d[5] * d[10] * d[29] - 2 * d[4] * d[11] * d[29] -
        2 * d[2] * d[13] * d[29] + 2 * d[1] * d[14] * d[29] + 2 * d[1] * d[9] * d[30] + 2 * d[0] * d[10] * d[30] +
        2 * d[4] * d[12] * d[30] + 2 * d[3] * d[13] * d[30] + 2 * d[7] * d[15] * d[30] + 2 * d[6] * d[16] * d[30] -
        2 * d[0] * d[9] * d[31] + 2 * d[1] * d[10] * d[31] - 2 * d[2] * d[11] * d[31] + 2 * d[3] * d[12] * d[31] +
        6 * d[4] * d[13] * d[31] + 2 * d[5] * d[14] * d[31] - 2 * d[6] * d[15] * d[31] + 2 * d[7] * d[16] * d[31] -
        2 * d[8] * d[17] * d[31] + 2 * d[2] * d[10] * d[32] + 2 * d[1] * d[11] * d[32] + 2 * d[5] * d[13] * d[32] +
        2 * d[4] * d[14] * d[32] + 2 * d[8] * d[16] * d[32] + 2 * d[7] * d[17] * d[32] + 2 * d[7] * d[12] * d[33] -
        2 * d[6] * d[13] * d[33] - 2 * d[4] * d[15] * d[33] + 2 * d[3] * d[16] * d[33] + 2 * d[6] * d[12] * d[34] +
        2 * d[7] * d[13] * d[34] + 2 * d[8] * d[14] * d[34] + 2 * d[3] * d[15] * d[34] + 2 * d[4] * d[16] * d[34] +
        2 * d[5] * d[17] * d[34] - 2 * d[8] * d[13] * d[35] + 2 * d[7] * d[14] * d[35] + 2 * d[5] * d[16] * d[35] -
        2 * d[4] * d[17] * d[35];
    coeffs[222] =
        2 * d[10] * d[12] * d[27] - 2 * d[9] * d[13] * d[27] + 2 * d[9] * d[12] * d[28] + 2 * d[10] * d[13] * d[28] +
        2 * d[11] * d[14] * d[28] - 2 * d[11] * d[13] * d[29] + 2 * d[10] * d[14] * d[29] + 2 * d[9] * d[10] * d[30] +
        2 * d[12] * d[13] * d[30] + 2 * d[15] * d[16] * d[30] - std::pow(d[9], 2) * d[31] + std::pow(d[10], 2) * d[31] -
        std::pow(d[11], 2) * d[31] + std::pow(d[12], 2) * d[31] + 3 * std::pow(d[13], 2) * d[31] +
        std::pow(d[14], 2) * d[31] - std::pow(d[15], 2) * d[31] + std::pow(d[16], 2) * d[31] -
        std::pow(d[17], 2) * d[31] + 2 * d[10] * d[11] * d[32] + 2 * d[13] * d[14] * d[32] + 2 * d[16] * d[17] * d[32] -
        2 * d[13] * d[15] * d[33] + 2 * d[12] * d[16] * d[33] + 2 * d[12] * d[15] * d[34] + 2 * d[13] * d[16] * d[34] +
        2 * d[14] * d[17] * d[34] + 2 * d[14] * d[16] * d[35] - 2 * d[13] * d[17] * d[35];
    coeffs[223] =
        -2 * d[4] * d[18] * d[27] + 2 * d[3] * d[19] * d[27] + 2 * d[1] * d[21] * d[27] - 2 * d[0] * d[22] * d[27] +
        2 * d[3] * d[18] * d[28] + 2 * d[4] * d[19] * d[28] + 2 * d[5] * d[20] * d[28] + 2 * d[0] * d[21] * d[28] +
        2 * d[1] * d[22] * d[28] + 2 * d[2] * d[23] * d[28] + 2 * d[5] * d[19] * d[29] - 2 * d[4] * d[20] * d[29] -
        2 * d[2] * d[22] * d[29] + 2 * d[1] * d[23] * d[29] + 2 * d[1] * d[18] * d[30] + 2 * d[0] * d[19] * d[30] +
        2 * d[4] * d[21] * d[30] + 2 * d[3] * d[22] * d[30] + 2 * d[7] * d[24] * d[30] + 2 * d[6] * d[25] * d[30] -
        2 * d[0] * d[18] * d[31] + 2 * d[1] * d[19] * d[31] - 2 * d[2] * d[20] * d[31] + 2 * d[3] * d[21] * d[31] +
        6 * d[4] * d[22] * d[31] + 2 * d[5] * d[23] * d[31] - 2 * d[6] * d[24] * d[31] + 2 * d[7] * d[25] * d[31] -
        2 * d[8] * d[26] * d[31] + 2 * d[2] * d[19] * d[32] + 2 * d[1] * d[20] * d[32] + 2 * d[5] * d[22] * d[32] +
        2 * d[4] * d[23] * d[32] + 2 * d[8] * d[25] * d[32] + 2 * d[7] * d[26] * d[32] + 2 * d[7] * d[21] * d[33] -
        2 * d[6] * d[22] * d[33] - 2 * d[4] * d[24] * d[33] + 2 * d[3] * d[25] * d[33] + 2 * d[6] * d[21] * d[34] +
        2 * d[7] * d[22] * d[34] + 2 * d[8] * d[23] * d[34] + 2 * d[3] * d[24] * d[34] + 2 * d[4] * d[25] * d[34] +
        2 * d[5] * d[26] * d[34] - 2 * d[8] * d[22] * d[35] + 2 * d[7] * d[23] * d[35] + 2 * d[5] * d[25] * d[35] -
        2 * d[4] * d[26] * d[35];
    coeffs[224] =
        -2 * d[13] * d[18] * d[27] + 2 * d[12] * d[19] * d[27] + 2 * d[10] * d[21] * d[27] - 2 * d[9] * d[22] * d[27] +
        2 * d[12] * d[18] * d[28] + 2 * d[13] * d[19] * d[28] + 2 * d[14] * d[20] * d[28] + 2 * d[9] * d[21] * d[28] +
        2 * d[10] * d[22] * d[28] + 2 * d[11] * d[23] * d[28] + 2 * d[14] * d[19] * d[29] - 2 * d[13] * d[20] * d[29] -
        2 * d[11] * d[22] * d[29] + 2 * d[10] * d[23] * d[29] + 2 * d[10] * d[18] * d[30] + 2 * d[9] * d[19] * d[30] +
        2 * d[13] * d[21] * d[30] + 2 * d[12] * d[22] * d[30] + 2 * d[16] * d[24] * d[30] + 2 * d[15] * d[25] * d[30] -
        2 * d[9] * d[18] * d[31] + 2 * d[10] * d[19] * d[31] - 2 * d[11] * d[20] * d[31] + 2 * d[12] * d[21] * d[31] +
        6 * d[13] * d[22] * d[31] + 2 * d[14] * d[23] * d[31] - 2 * d[15] * d[24] * d[31] + 2 * d[16] * d[25] * d[31] -
        2 * d[17] * d[26] * d[31] + 2 * d[11] * d[19] * d[32] + 2 * d[10] * d[20] * d[32] + 2 * d[14] * d[22] * d[32] +
        2 * d[13] * d[23] * d[32] + 2 * d[17] * d[25] * d[32] + 2 * d[16] * d[26] * d[32] + 2 * d[16] * d[21] * d[33] -
        2 * d[15] * d[22] * d[33] - 2 * d[13] * d[24] * d[33] + 2 * d[12] * d[25] * d[33] + 2 * d[15] * d[21] * d[34] +
        2 * d[16] * d[22] * d[34] + 2 * d[17] * d[23] * d[34] + 2 * d[12] * d[24] * d[34] + 2 * d[13] * d[25] * d[34] +
        2 * d[14] * d[26] * d[34] - 2 * d[17] * d[22] * d[35] + 2 * d[16] * d[23] * d[35] + 2 * d[14] * d[25] * d[35] -
        2 * d[13] * d[26] * d[35];
    coeffs[225] = 2 * d[19] * d[21] * d[27] - 2 * d[18] * d[22] * d[27] + 2 * d[18] * d[21] * d[28] +
                  2 * d[19] * d[22] * d[28] + 2 * d[20] * d[23] * d[28] - 2 * d[20] * d[22] * d[29] +
                  2 * d[19] * d[23] * d[29] + 2 * d[18] * d[19] * d[30] + 2 * d[21] * d[22] * d[30] +
                  2 * d[24] * d[25] * d[30] - std::pow(d[18], 2) * d[31] + std::pow(d[19], 2) * d[31] -
                  std::pow(d[20], 2) * d[31] + std::pow(d[21], 2) * d[31] + 3 * std::pow(d[22], 2) * d[31] +
                  std::pow(d[23], 2) * d[31] - std::pow(d[24], 2) * d[31] + std::pow(d[25], 2) * d[31] -
                  std::pow(d[26], 2) * d[31] + 2 * d[19] * d[20] * d[32] + 2 * d[22] * d[23] * d[32] +
                  2 * d[25] * d[26] * d[32] - 2 * d[22] * d[24] * d[33] + 2 * d[21] * d[25] * d[33] +
                  2 * d[21] * d[24] * d[34] + 2 * d[22] * d[25] * d[34] + 2 * d[23] * d[26] * d[34] +
                  2 * d[23] * d[25] * d[35] - 2 * d[22] * d[26] * d[35];
    coeffs[226] =
        -d[4] * std::pow(d[27], 2) + 2 * d[3] * d[27] * d[28] + d[4] * std::pow(d[28], 2) + 2 * d[5] * d[28] * d[29] -
        d[4] * std::pow(d[29], 2) + 2 * d[1] * d[27] * d[30] + 2 * d[0] * d[28] * d[30] + d[4] * std::pow(d[30], 2) -
        2 * d[0] * d[27] * d[31] + 2 * d[1] * d[28] * d[31] - 2 * d[2] * d[29] * d[31] + 2 * d[3] * d[30] * d[31] +
        3 * d[4] * std::pow(d[31], 2) + 2 * d[2] * d[28] * d[32] + 2 * d[1] * d[29] * d[32] + 2 * d[5] * d[31] * d[32] +
        d[4] * std::pow(d[32], 2) + 2 * d[7] * d[30] * d[33] - 2 * d[6] * d[31] * d[33] - d[4] * std::pow(d[33], 2) +
        2 * d[6] * d[30] * d[34] + 2 * d[7] * d[31] * d[34] + 2 * d[8] * d[32] * d[34] + 2 * d[3] * d[33] * d[34] +
        d[4] * std::pow(d[34], 2) - 2 * d[8] * d[31] * d[35] + 2 * d[7] * d[32] * d[35] + 2 * d[5] * d[34] * d[35] -
        d[4] * std::pow(d[35], 2);
    coeffs[227] =
        -d[13] * std::pow(d[27], 2) + 2 * d[12] * d[27] * d[28] + d[13] * std::pow(d[28], 2) +
        2 * d[14] * d[28] * d[29] - d[13] * std::pow(d[29], 2) + 2 * d[10] * d[27] * d[30] + 2 * d[9] * d[28] * d[30] +
        d[13] * std::pow(d[30], 2) - 2 * d[9] * d[27] * d[31] + 2 * d[10] * d[28] * d[31] - 2 * d[11] * d[29] * d[31] +
        2 * d[12] * d[30] * d[31] + 3 * d[13] * std::pow(d[31], 2) + 2 * d[11] * d[28] * d[32] +
        2 * d[10] * d[29] * d[32] + 2 * d[14] * d[31] * d[32] + d[13] * std::pow(d[32], 2) + 2 * d[16] * d[30] * d[33] -
        2 * d[15] * d[31] * d[33] - d[13] * std::pow(d[33], 2) + 2 * d[15] * d[30] * d[34] + 2 * d[16] * d[31] * d[34] +
        2 * d[17] * d[32] * d[34] + 2 * d[12] * d[33] * d[34] + d[13] * std::pow(d[34], 2) - 2 * d[17] * d[31] * d[35] +
        2 * d[16] * d[32] * d[35] + 2 * d[14] * d[34] * d[35] - d[13] * std::pow(d[35], 2);
    coeffs[228] =
        -d[22] * std::pow(d[27], 2) + 2 * d[21] * d[27] * d[28] + d[22] * std::pow(d[28], 2) +
        2 * d[23] * d[28] * d[29] - d[22] * std::pow(d[29], 2) + 2 * d[19] * d[27] * d[30] + 2 * d[18] * d[28] * d[30] +
        d[22] * std::pow(d[30], 2) - 2 * d[18] * d[27] * d[31] + 2 * d[19] * d[28] * d[31] - 2 * d[20] * d[29] * d[31] +
        2 * d[21] * d[30] * d[31] + 3 * d[22] * std::pow(d[31], 2) + 2 * d[20] * d[28] * d[32] +
        2 * d[19] * d[29] * d[32] + 2 * d[23] * d[31] * d[32] + d[22] * std::pow(d[32], 2) + 2 * d[25] * d[30] * d[33] -
        2 * d[24] * d[31] * d[33] - d[22] * std::pow(d[33], 2) + 2 * d[24] * d[30] * d[34] + 2 * d[25] * d[31] * d[34] +
        2 * d[26] * d[32] * d[34] + 2 * d[21] * d[33] * d[34] + d[22] * std::pow(d[34], 2) - 2 * d[26] * d[31] * d[35] +
        2 * d[25] * d[32] * d[35] + 2 * d[23] * d[34] * d[35] - d[22] * std::pow(d[35], 2);
    coeffs[229] = 2 * d[27] * d[28] * d[30] - std::pow(d[27], 2) * d[31] + std::pow(d[28], 2) * d[31] -
                  std::pow(d[29], 2) * d[31] + std::pow(d[30], 2) * d[31] + std::pow(d[31], 3) +
                  2 * d[28] * d[29] * d[32] + d[31] * std::pow(d[32], 2) - d[31] * std::pow(d[33], 2) +
                  2 * d[30] * d[33] * d[34] + d[31] * std::pow(d[34], 2) + 2 * d[32] * d[34] * d[35] -
                  d[31] * std::pow(d[35], 2);
    coeffs[230] =
        2 * d[1] * d[3] * d[36] - 2 * d[0] * d[4] * d[36] + 2 * d[0] * d[3] * d[37] + 2 * d[1] * d[4] * d[37] +
        2 * d[2] * d[5] * d[37] - 2 * d[2] * d[4] * d[38] + 2 * d[1] * d[5] * d[38] + 2 * d[0] * d[1] * d[39] +
        2 * d[3] * d[4] * d[39] + 2 * d[6] * d[7] * d[39] - std::pow(d[0], 2) * d[40] + std::pow(d[1], 2) * d[40] -
        std::pow(d[2], 2) * d[40] + std::pow(d[3], 2) * d[40] + 3 * std::pow(d[4], 2) * d[40] +
        std::pow(d[5], 2) * d[40] - std::pow(d[6], 2) * d[40] + std::pow(d[7], 2) * d[40] - std::pow(d[8], 2) * d[40] +
        2 * d[1] * d[2] * d[41] + 2 * d[4] * d[5] * d[41] + 2 * d[7] * d[8] * d[41] - 2 * d[4] * d[6] * d[42] +
        2 * d[3] * d[7] * d[42] + 2 * d[3] * d[6] * d[43] + 2 * d[4] * d[7] * d[43] + 2 * d[5] * d[8] * d[43] +
        2 * d[5] * d[7] * d[44] - 2 * d[4] * d[8] * d[44];
    coeffs[231] =
        -2 * d[4] * d[9] * d[36] + 2 * d[3] * d[10] * d[36] + 2 * d[1] * d[12] * d[36] - 2 * d[0] * d[13] * d[36] +
        2 * d[3] * d[9] * d[37] + 2 * d[4] * d[10] * d[37] + 2 * d[5] * d[11] * d[37] + 2 * d[0] * d[12] * d[37] +
        2 * d[1] * d[13] * d[37] + 2 * d[2] * d[14] * d[37] + 2 * d[5] * d[10] * d[38] - 2 * d[4] * d[11] * d[38] -
        2 * d[2] * d[13] * d[38] + 2 * d[1] * d[14] * d[38] + 2 * d[1] * d[9] * d[39] + 2 * d[0] * d[10] * d[39] +
        2 * d[4] * d[12] * d[39] + 2 * d[3] * d[13] * d[39] + 2 * d[7] * d[15] * d[39] + 2 * d[6] * d[16] * d[39] -
        2 * d[0] * d[9] * d[40] + 2 * d[1] * d[10] * d[40] - 2 * d[2] * d[11] * d[40] + 2 * d[3] * d[12] * d[40] +
        6 * d[4] * d[13] * d[40] + 2 * d[5] * d[14] * d[40] - 2 * d[6] * d[15] * d[40] + 2 * d[7] * d[16] * d[40] -
        2 * d[8] * d[17] * d[40] + 2 * d[2] * d[10] * d[41] + 2 * d[1] * d[11] * d[41] + 2 * d[5] * d[13] * d[41] +
        2 * d[4] * d[14] * d[41] + 2 * d[8] * d[16] * d[41] + 2 * d[7] * d[17] * d[41] + 2 * d[7] * d[12] * d[42] -
        2 * d[6] * d[13] * d[42] - 2 * d[4] * d[15] * d[42] + 2 * d[3] * d[16] * d[42] + 2 * d[6] * d[12] * d[43] +
        2 * d[7] * d[13] * d[43] + 2 * d[8] * d[14] * d[43] + 2 * d[3] * d[15] * d[43] + 2 * d[4] * d[16] * d[43] +
        2 * d[5] * d[17] * d[43] - 2 * d[8] * d[13] * d[44] + 2 * d[7] * d[14] * d[44] + 2 * d[5] * d[16] * d[44] -
        2 * d[4] * d[17] * d[44];
    coeffs[232] =
        2 * d[10] * d[12] * d[36] - 2 * d[9] * d[13] * d[36] + 2 * d[9] * d[12] * d[37] + 2 * d[10] * d[13] * d[37] +
        2 * d[11] * d[14] * d[37] - 2 * d[11] * d[13] * d[38] + 2 * d[10] * d[14] * d[38] + 2 * d[9] * d[10] * d[39] +
        2 * d[12] * d[13] * d[39] + 2 * d[15] * d[16] * d[39] - std::pow(d[9], 2) * d[40] + std::pow(d[10], 2) * d[40] -
        std::pow(d[11], 2) * d[40] + std::pow(d[12], 2) * d[40] + 3 * std::pow(d[13], 2) * d[40] +
        std::pow(d[14], 2) * d[40] - std::pow(d[15], 2) * d[40] + std::pow(d[16], 2) * d[40] -
        std::pow(d[17], 2) * d[40] + 2 * d[10] * d[11] * d[41] + 2 * d[13] * d[14] * d[41] + 2 * d[16] * d[17] * d[41] -
        2 * d[13] * d[15] * d[42] + 2 * d[12] * d[16] * d[42] + 2 * d[12] * d[15] * d[43] + 2 * d[13] * d[16] * d[43] +
        2 * d[14] * d[17] * d[43] + 2 * d[14] * d[16] * d[44] - 2 * d[13] * d[17] * d[44];
    coeffs[233] =
        -2 * d[4] * d[18] * d[36] + 2 * d[3] * d[19] * d[36] + 2 * d[1] * d[21] * d[36] - 2 * d[0] * d[22] * d[36] +
        2 * d[3] * d[18] * d[37] + 2 * d[4] * d[19] * d[37] + 2 * d[5] * d[20] * d[37] + 2 * d[0] * d[21] * d[37] +
        2 * d[1] * d[22] * d[37] + 2 * d[2] * d[23] * d[37] + 2 * d[5] * d[19] * d[38] - 2 * d[4] * d[20] * d[38] -
        2 * d[2] * d[22] * d[38] + 2 * d[1] * d[23] * d[38] + 2 * d[1] * d[18] * d[39] + 2 * d[0] * d[19] * d[39] +
        2 * d[4] * d[21] * d[39] + 2 * d[3] * d[22] * d[39] + 2 * d[7] * d[24] * d[39] + 2 * d[6] * d[25] * d[39] -
        2 * d[0] * d[18] * d[40] + 2 * d[1] * d[19] * d[40] - 2 * d[2] * d[20] * d[40] + 2 * d[3] * d[21] * d[40] +
        6 * d[4] * d[22] * d[40] + 2 * d[5] * d[23] * d[40] - 2 * d[6] * d[24] * d[40] + 2 * d[7] * d[25] * d[40] -
        2 * d[8] * d[26] * d[40] + 2 * d[2] * d[19] * d[41] + 2 * d[1] * d[20] * d[41] + 2 * d[5] * d[22] * d[41] +
        2 * d[4] * d[23] * d[41] + 2 * d[8] * d[25] * d[41] + 2 * d[7] * d[26] * d[41] + 2 * d[7] * d[21] * d[42] -
        2 * d[6] * d[22] * d[42] - 2 * d[4] * d[24] * d[42] + 2 * d[3] * d[25] * d[42] + 2 * d[6] * d[21] * d[43] +
        2 * d[7] * d[22] * d[43] + 2 * d[8] * d[23] * d[43] + 2 * d[3] * d[24] * d[43] + 2 * d[4] * d[25] * d[43] +
        2 * d[5] * d[26] * d[43] - 2 * d[8] * d[22] * d[44] + 2 * d[7] * d[23] * d[44] + 2 * d[5] * d[25] * d[44] -
        2 * d[4] * d[26] * d[44];
    coeffs[234] =
        -2 * d[13] * d[18] * d[36] + 2 * d[12] * d[19] * d[36] + 2 * d[10] * d[21] * d[36] - 2 * d[9] * d[22] * d[36] +
        2 * d[12] * d[18] * d[37] + 2 * d[13] * d[19] * d[37] + 2 * d[14] * d[20] * d[37] + 2 * d[9] * d[21] * d[37] +
        2 * d[10] * d[22] * d[37] + 2 * d[11] * d[23] * d[37] + 2 * d[14] * d[19] * d[38] - 2 * d[13] * d[20] * d[38] -
        2 * d[11] * d[22] * d[38] + 2 * d[10] * d[23] * d[38] + 2 * d[10] * d[18] * d[39] + 2 * d[9] * d[19] * d[39] +
        2 * d[13] * d[21] * d[39] + 2 * d[12] * d[22] * d[39] + 2 * d[16] * d[24] * d[39] + 2 * d[15] * d[25] * d[39] -
        2 * d[9] * d[18] * d[40] + 2 * d[10] * d[19] * d[40] - 2 * d[11] * d[20] * d[40] + 2 * d[12] * d[21] * d[40] +
        6 * d[13] * d[22] * d[40] + 2 * d[14] * d[23] * d[40] - 2 * d[15] * d[24] * d[40] + 2 * d[16] * d[25] * d[40] -
        2 * d[17] * d[26] * d[40] + 2 * d[11] * d[19] * d[41] + 2 * d[10] * d[20] * d[41] + 2 * d[14] * d[22] * d[41] +
        2 * d[13] * d[23] * d[41] + 2 * d[17] * d[25] * d[41] + 2 * d[16] * d[26] * d[41] + 2 * d[16] * d[21] * d[42] -
        2 * d[15] * d[22] * d[42] - 2 * d[13] * d[24] * d[42] + 2 * d[12] * d[25] * d[42] + 2 * d[15] * d[21] * d[43] +
        2 * d[16] * d[22] * d[43] + 2 * d[17] * d[23] * d[43] + 2 * d[12] * d[24] * d[43] + 2 * d[13] * d[25] * d[43] +
        2 * d[14] * d[26] * d[43] - 2 * d[17] * d[22] * d[44] + 2 * d[16] * d[23] * d[44] + 2 * d[14] * d[25] * d[44] -
        2 * d[13] * d[26] * d[44];
    coeffs[235] = 2 * d[19] * d[21] * d[36] - 2 * d[18] * d[22] * d[36] + 2 * d[18] * d[21] * d[37] +
                  2 * d[19] * d[22] * d[37] + 2 * d[20] * d[23] * d[37] - 2 * d[20] * d[22] * d[38] +
                  2 * d[19] * d[23] * d[38] + 2 * d[18] * d[19] * d[39] + 2 * d[21] * d[22] * d[39] +
                  2 * d[24] * d[25] * d[39] - std::pow(d[18], 2) * d[40] + std::pow(d[19], 2) * d[40] -
                  std::pow(d[20], 2) * d[40] + std::pow(d[21], 2) * d[40] + 3 * std::pow(d[22], 2) * d[40] +
                  std::pow(d[23], 2) * d[40] - std::pow(d[24], 2) * d[40] + std::pow(d[25], 2) * d[40] -
                  std::pow(d[26], 2) * d[40] + 2 * d[19] * d[20] * d[41] + 2 * d[22] * d[23] * d[41] +
                  2 * d[25] * d[26] * d[41] - 2 * d[22] * d[24] * d[42] + 2 * d[21] * d[25] * d[42] +
                  2 * d[21] * d[24] * d[43] + 2 * d[22] * d[25] * d[43] + 2 * d[23] * d[26] * d[43] +
                  2 * d[23] * d[25] * d[44] - 2 * d[22] * d[26] * d[44];
    coeffs[236] =
        -2 * d[4] * d[27] * d[36] + 2 * d[3] * d[28] * d[36] + 2 * d[1] * d[30] * d[36] - 2 * d[0] * d[31] * d[36] +
        2 * d[3] * d[27] * d[37] + 2 * d[4] * d[28] * d[37] + 2 * d[5] * d[29] * d[37] + 2 * d[0] * d[30] * d[37] +
        2 * d[1] * d[31] * d[37] + 2 * d[2] * d[32] * d[37] + 2 * d[5] * d[28] * d[38] - 2 * d[4] * d[29] * d[38] -
        2 * d[2] * d[31] * d[38] + 2 * d[1] * d[32] * d[38] + 2 * d[1] * d[27] * d[39] + 2 * d[0] * d[28] * d[39] +
        2 * d[4] * d[30] * d[39] + 2 * d[3] * d[31] * d[39] + 2 * d[7] * d[33] * d[39] + 2 * d[6] * d[34] * d[39] -
        2 * d[0] * d[27] * d[40] + 2 * d[1] * d[28] * d[40] - 2 * d[2] * d[29] * d[40] + 2 * d[3] * d[30] * d[40] +
        6 * d[4] * d[31] * d[40] + 2 * d[5] * d[32] * d[40] - 2 * d[6] * d[33] * d[40] + 2 * d[7] * d[34] * d[40] -
        2 * d[8] * d[35] * d[40] + 2 * d[2] * d[28] * d[41] + 2 * d[1] * d[29] * d[41] + 2 * d[5] * d[31] * d[41] +
        2 * d[4] * d[32] * d[41] + 2 * d[8] * d[34] * d[41] + 2 * d[7] * d[35] * d[41] + 2 * d[7] * d[30] * d[42] -
        2 * d[6] * d[31] * d[42] - 2 * d[4] * d[33] * d[42] + 2 * d[3] * d[34] * d[42] + 2 * d[6] * d[30] * d[43] +
        2 * d[7] * d[31] * d[43] + 2 * d[8] * d[32] * d[43] + 2 * d[3] * d[33] * d[43] + 2 * d[4] * d[34] * d[43] +
        2 * d[5] * d[35] * d[43] - 2 * d[8] * d[31] * d[44] + 2 * d[7] * d[32] * d[44] + 2 * d[5] * d[34] * d[44] -
        2 * d[4] * d[35] * d[44];
    coeffs[237] =
        -2 * d[13] * d[27] * d[36] + 2 * d[12] * d[28] * d[36] + 2 * d[10] * d[30] * d[36] - 2 * d[9] * d[31] * d[36] +
        2 * d[12] * d[27] * d[37] + 2 * d[13] * d[28] * d[37] + 2 * d[14] * d[29] * d[37] + 2 * d[9] * d[30] * d[37] +
        2 * d[10] * d[31] * d[37] + 2 * d[11] * d[32] * d[37] + 2 * d[14] * d[28] * d[38] - 2 * d[13] * d[29] * d[38] -
        2 * d[11] * d[31] * d[38] + 2 * d[10] * d[32] * d[38] + 2 * d[10] * d[27] * d[39] + 2 * d[9] * d[28] * d[39] +
        2 * d[13] * d[30] * d[39] + 2 * d[12] * d[31] * d[39] + 2 * d[16] * d[33] * d[39] + 2 * d[15] * d[34] * d[39] -
        2 * d[9] * d[27] * d[40] + 2 * d[10] * d[28] * d[40] - 2 * d[11] * d[29] * d[40] + 2 * d[12] * d[30] * d[40] +
        6 * d[13] * d[31] * d[40] + 2 * d[14] * d[32] * d[40] - 2 * d[15] * d[33] * d[40] + 2 * d[16] * d[34] * d[40] -
        2 * d[17] * d[35] * d[40] + 2 * d[11] * d[28] * d[41] + 2 * d[10] * d[29] * d[41] + 2 * d[14] * d[31] * d[41] +
        2 * d[13] * d[32] * d[41] + 2 * d[17] * d[34] * d[41] + 2 * d[16] * d[35] * d[41] + 2 * d[16] * d[30] * d[42] -
        2 * d[15] * d[31] * d[42] - 2 * d[13] * d[33] * d[42] + 2 * d[12] * d[34] * d[42] + 2 * d[15] * d[30] * d[43] +
        2 * d[16] * d[31] * d[43] + 2 * d[17] * d[32] * d[43] + 2 * d[12] * d[33] * d[43] + 2 * d[13] * d[34] * d[43] +
        2 * d[14] * d[35] * d[43] - 2 * d[17] * d[31] * d[44] + 2 * d[16] * d[32] * d[44] + 2 * d[14] * d[34] * d[44] -
        2 * d[13] * d[35] * d[44];
    coeffs[238] =
        -2 * d[22] * d[27] * d[36] + 2 * d[21] * d[28] * d[36] + 2 * d[19] * d[30] * d[36] - 2 * d[18] * d[31] * d[36] +
        2 * d[21] * d[27] * d[37] + 2 * d[22] * d[28] * d[37] + 2 * d[23] * d[29] * d[37] + 2 * d[18] * d[30] * d[37] +
        2 * d[19] * d[31] * d[37] + 2 * d[20] * d[32] * d[37] + 2 * d[23] * d[28] * d[38] - 2 * d[22] * d[29] * d[38] -
        2 * d[20] * d[31] * d[38] + 2 * d[19] * d[32] * d[38] + 2 * d[19] * d[27] * d[39] + 2 * d[18] * d[28] * d[39] +
        2 * d[22] * d[30] * d[39] + 2 * d[21] * d[31] * d[39] + 2 * d[25] * d[33] * d[39] + 2 * d[24] * d[34] * d[39] -
        2 * d[18] * d[27] * d[40] + 2 * d[19] * d[28] * d[40] - 2 * d[20] * d[29] * d[40] + 2 * d[21] * d[30] * d[40] +
        6 * d[22] * d[31] * d[40] + 2 * d[23] * d[32] * d[40] - 2 * d[24] * d[33] * d[40] + 2 * d[25] * d[34] * d[40] -
        2 * d[26] * d[35] * d[40] + 2 * d[20] * d[28] * d[41] + 2 * d[19] * d[29] * d[41] + 2 * d[23] * d[31] * d[41] +
        2 * d[22] * d[32] * d[41] + 2 * d[26] * d[34] * d[41] + 2 * d[25] * d[35] * d[41] + 2 * d[25] * d[30] * d[42] -
        2 * d[24] * d[31] * d[42] - 2 * d[22] * d[33] * d[42] + 2 * d[21] * d[34] * d[42] + 2 * d[24] * d[30] * d[43] +
        2 * d[25] * d[31] * d[43] + 2 * d[26] * d[32] * d[43] + 2 * d[21] * d[33] * d[43] + 2 * d[22] * d[34] * d[43] +
        2 * d[23] * d[35] * d[43] - 2 * d[26] * d[31] * d[44] + 2 * d[25] * d[32] * d[44] + 2 * d[23] * d[34] * d[44] -
        2 * d[22] * d[35] * d[44];
    coeffs[239] = 2 * d[28] * d[30] * d[36] - 2 * d[27] * d[31] * d[36] + 2 * d[27] * d[30] * d[37] +
                  2 * d[28] * d[31] * d[37] + 2 * d[29] * d[32] * d[37] - 2 * d[29] * d[31] * d[38] +
                  2 * d[28] * d[32] * d[38] + 2 * d[27] * d[28] * d[39] + 2 * d[30] * d[31] * d[39] +
                  2 * d[33] * d[34] * d[39] - std::pow(d[27], 2) * d[40] + std::pow(d[28], 2) * d[40] -
                  std::pow(d[29], 2) * d[40] + std::pow(d[30], 2) * d[40] + 3 * std::pow(d[31], 2) * d[40] +
                  std::pow(d[32], 2) * d[40] - std::pow(d[33], 2) * d[40] + std::pow(d[34], 2) * d[40] -
                  std::pow(d[35], 2) * d[40] + 2 * d[28] * d[29] * d[41] + 2 * d[31] * d[32] * d[41] +
                  2 * d[34] * d[35] * d[41] - 2 * d[31] * d[33] * d[42] + 2 * d[30] * d[34] * d[42] +
                  2 * d[30] * d[33] * d[43] + 2 * d[31] * d[34] * d[43] + 2 * d[32] * d[35] * d[43] +
                  2 * d[32] * d[34] * d[44] - 2 * d[31] * d[35] * d[44];
    coeffs[240] =
        -d[4] * std::pow(d[36], 2) + 2 * d[3] * d[36] * d[37] + d[4] * std::pow(d[37], 2) + 2 * d[5] * d[37] * d[38] -
        d[4] * std::pow(d[38], 2) + 2 * d[1] * d[36] * d[39] + 2 * d[0] * d[37] * d[39] + d[4] * std::pow(d[39], 2) -
        2 * d[0] * d[36] * d[40] + 2 * d[1] * d[37] * d[40] - 2 * d[2] * d[38] * d[40] + 2 * d[3] * d[39] * d[40] +
        3 * d[4] * std::pow(d[40], 2) + 2 * d[2] * d[37] * d[41] + 2 * d[1] * d[38] * d[41] + 2 * d[5] * d[40] * d[41] +
        d[4] * std::pow(d[41], 2) + 2 * d[7] * d[39] * d[42] - 2 * d[6] * d[40] * d[42] - d[4] * std::pow(d[42], 2) +
        2 * d[6] * d[39] * d[43] + 2 * d[7] * d[40] * d[43] + 2 * d[8] * d[41] * d[43] + 2 * d[3] * d[42] * d[43] +
        d[4] * std::pow(d[43], 2) - 2 * d[8] * d[40] * d[44] + 2 * d[7] * d[41] * d[44] + 2 * d[5] * d[43] * d[44] -
        d[4] * std::pow(d[44], 2);
    coeffs[241] =
        -d[13] * std::pow(d[36], 2) + 2 * d[12] * d[36] * d[37] + d[13] * std::pow(d[37], 2) +
        2 * d[14] * d[37] * d[38] - d[13] * std::pow(d[38], 2) + 2 * d[10] * d[36] * d[39] + 2 * d[9] * d[37] * d[39] +
        d[13] * std::pow(d[39], 2) - 2 * d[9] * d[36] * d[40] + 2 * d[10] * d[37] * d[40] - 2 * d[11] * d[38] * d[40] +
        2 * d[12] * d[39] * d[40] + 3 * d[13] * std::pow(d[40], 2) + 2 * d[11] * d[37] * d[41] +
        2 * d[10] * d[38] * d[41] + 2 * d[14] * d[40] * d[41] + d[13] * std::pow(d[41], 2) + 2 * d[16] * d[39] * d[42] -
        2 * d[15] * d[40] * d[42] - d[13] * std::pow(d[42], 2) + 2 * d[15] * d[39] * d[43] + 2 * d[16] * d[40] * d[43] +
        2 * d[17] * d[41] * d[43] + 2 * d[12] * d[42] * d[43] + d[13] * std::pow(d[43], 2) - 2 * d[17] * d[40] * d[44] +
        2 * d[16] * d[41] * d[44] + 2 * d[14] * d[43] * d[44] - d[13] * std::pow(d[44], 2);
    coeffs[242] =
        -d[22] * std::pow(d[36], 2) + 2 * d[21] * d[36] * d[37] + d[22] * std::pow(d[37], 2) +
        2 * d[23] * d[37] * d[38] - d[22] * std::pow(d[38], 2) + 2 * d[19] * d[36] * d[39] + 2 * d[18] * d[37] * d[39] +
        d[22] * std::pow(d[39], 2) - 2 * d[18] * d[36] * d[40] + 2 * d[19] * d[37] * d[40] - 2 * d[20] * d[38] * d[40] +
        2 * d[21] * d[39] * d[40] + 3 * d[22] * std::pow(d[40], 2) + 2 * d[20] * d[37] * d[41] +
        2 * d[19] * d[38] * d[41] + 2 * d[23] * d[40] * d[41] + d[22] * std::pow(d[41], 2) + 2 * d[25] * d[39] * d[42] -
        2 * d[24] * d[40] * d[42] - d[22] * std::pow(d[42], 2) + 2 * d[24] * d[39] * d[43] + 2 * d[25] * d[40] * d[43] +
        2 * d[26] * d[41] * d[43] + 2 * d[21] * d[42] * d[43] + d[22] * std::pow(d[43], 2) - 2 * d[26] * d[40] * d[44] +
        2 * d[25] * d[41] * d[44] + 2 * d[23] * d[43] * d[44] - d[22] * std::pow(d[44], 2);
    coeffs[243] =
        -d[31] * std::pow(d[36], 2) + 2 * d[30] * d[36] * d[37] + d[31] * std::pow(d[37], 2) +
        2 * d[32] * d[37] * d[38] - d[31] * std::pow(d[38], 2) + 2 * d[28] * d[36] * d[39] + 2 * d[27] * d[37] * d[39] +
        d[31] * std::pow(d[39], 2) - 2 * d[27] * d[36] * d[40] + 2 * d[28] * d[37] * d[40] - 2 * d[29] * d[38] * d[40] +
        2 * d[30] * d[39] * d[40] + 3 * d[31] * std::pow(d[40], 2) + 2 * d[29] * d[37] * d[41] +
        2 * d[28] * d[38] * d[41] + 2 * d[32] * d[40] * d[41] + d[31] * std::pow(d[41], 2) + 2 * d[34] * d[39] * d[42] -
        2 * d[33] * d[40] * d[42] - d[31] * std::pow(d[42], 2) + 2 * d[33] * d[39] * d[43] + 2 * d[34] * d[40] * d[43] +
        2 * d[35] * d[41] * d[43] + 2 * d[30] * d[42] * d[43] + d[31] * std::pow(d[43], 2) - 2 * d[35] * d[40] * d[44] +
        2 * d[34] * d[41] * d[44] + 2 * d[32] * d[43] * d[44] - d[31] * std::pow(d[44], 2);
    coeffs[244] = 2 * d[36] * d[37] * d[39] - std::pow(d[36], 2) * d[40] + std::pow(d[37], 2) * d[40] -
                  std::pow(d[38], 2) * d[40] + std::pow(d[39], 2) * d[40] + std::pow(d[40], 3) +
                  2 * d[37] * d[38] * d[41] + d[40] * std::pow(d[41], 2) - d[40] * std::pow(d[42], 2) +
                  2 * d[39] * d[42] * d[43] + d[40] * std::pow(d[43], 2) + 2 * d[41] * d[43] * d[44] -
                  d[40] * std::pow(d[44], 2);
    coeffs[245] = 2 * d[0] * d[2] * d[3] + 2 * d[1] * d[2] * d[4] - std::pow(d[0], 2) * d[5] -
                  std::pow(d[1], 2) * d[5] + std::pow(d[2], 2) * d[5] + std::pow(d[3], 2) * d[5] +
                  std::pow(d[4], 2) * d[5] + std::pow(d[5], 3) - d[5] * std::pow(d[6], 2) - d[5] * std::pow(d[7], 2) +
                  2 * d[3] * d[6] * d[8] + 2 * d[4] * d[7] * d[8] + d[5] * std::pow(d[8], 2);
    coeffs[246] = 2 * d[2] * d[3] * d[9] - 2 * d[0] * d[5] * d[9] + 2 * d[2] * d[4] * d[10] - 2 * d[1] * d[5] * d[10] +
                  2 * d[0] * d[3] * d[11] + 2 * d[1] * d[4] * d[11] + 2 * d[2] * d[5] * d[11] +
                  2 * d[0] * d[2] * d[12] + 2 * d[3] * d[5] * d[12] + 2 * d[6] * d[8] * d[12] +
                  2 * d[1] * d[2] * d[13] + 2 * d[4] * d[5] * d[13] + 2 * d[7] * d[8] * d[13] -
                  std::pow(d[0], 2) * d[14] - std::pow(d[1], 2) * d[14] + std::pow(d[2], 2) * d[14] +
                  std::pow(d[3], 2) * d[14] + std::pow(d[4], 2) * d[14] + 3 * std::pow(d[5], 2) * d[14] -
                  std::pow(d[6], 2) * d[14] - std::pow(d[7], 2) * d[14] + std::pow(d[8], 2) * d[14] -
                  2 * d[5] * d[6] * d[15] + 2 * d[3] * d[8] * d[15] - 2 * d[5] * d[7] * d[16] +
                  2 * d[4] * d[8] * d[16] + 2 * d[3] * d[6] * d[17] + 2 * d[4] * d[7] * d[17] + 2 * d[5] * d[8] * d[17];
    coeffs[247] =
        -d[5] * std::pow(d[9], 2) - d[5] * std::pow(d[10], 2) + 2 * d[3] * d[9] * d[11] + 2 * d[4] * d[10] * d[11] +
        d[5] * std::pow(d[11], 2) + 2 * d[2] * d[9] * d[12] + 2 * d[0] * d[11] * d[12] + d[5] * std::pow(d[12], 2) +
        2 * d[2] * d[10] * d[13] + 2 * d[1] * d[11] * d[13] + d[5] * std::pow(d[13], 2) - 2 * d[0] * d[9] * d[14] -
        2 * d[1] * d[10] * d[14] + 2 * d[2] * d[11] * d[14] + 2 * d[3] * d[12] * d[14] + 2 * d[4] * d[13] * d[14] +
        3 * d[5] * std::pow(d[14], 2) + 2 * d[8] * d[12] * d[15] - 2 * d[6] * d[14] * d[15] -
        d[5] * std::pow(d[15], 2) + 2 * d[8] * d[13] * d[16] - 2 * d[7] * d[14] * d[16] - d[5] * std::pow(d[16], 2) +
        2 * d[6] * d[12] * d[17] + 2 * d[7] * d[13] * d[17] + 2 * d[8] * d[14] * d[17] + 2 * d[3] * d[15] * d[17] +
        2 * d[4] * d[16] * d[17] + d[5] * std::pow(d[17], 2);
    coeffs[248] = 2 * d[9] * d[11] * d[12] + 2 * d[10] * d[11] * d[13] - std::pow(d[9], 2) * d[14] -
                  std::pow(d[10], 2) * d[14] + std::pow(d[11], 2) * d[14] + std::pow(d[12], 2) * d[14] +
                  std::pow(d[13], 2) * d[14] + std::pow(d[14], 3) - d[14] * std::pow(d[15], 2) -
                  d[14] * std::pow(d[16], 2) + 2 * d[12] * d[15] * d[17] + 2 * d[13] * d[16] * d[17] +
                  d[14] * std::pow(d[17], 2);
    coeffs[249] =
        2 * d[2] * d[3] * d[18] - 2 * d[0] * d[5] * d[18] + 2 * d[2] * d[4] * d[19] - 2 * d[1] * d[5] * d[19] +
        2 * d[0] * d[3] * d[20] + 2 * d[1] * d[4] * d[20] + 2 * d[2] * d[5] * d[20] + 2 * d[0] * d[2] * d[21] +
        2 * d[3] * d[5] * d[21] + 2 * d[6] * d[8] * d[21] + 2 * d[1] * d[2] * d[22] + 2 * d[4] * d[5] * d[22] +
        2 * d[7] * d[8] * d[22] - std::pow(d[0], 2) * d[23] - std::pow(d[1], 2) * d[23] + std::pow(d[2], 2) * d[23] +
        std::pow(d[3], 2) * d[23] + std::pow(d[4], 2) * d[23] + 3 * std::pow(d[5], 2) * d[23] -
        std::pow(d[6], 2) * d[23] - std::pow(d[7], 2) * d[23] + std::pow(d[8], 2) * d[23] - 2 * d[5] * d[6] * d[24] +
        2 * d[3] * d[8] * d[24] - 2 * d[5] * d[7] * d[25] + 2 * d[4] * d[8] * d[25] + 2 * d[3] * d[6] * d[26] +
        2 * d[4] * d[7] * d[26] + 2 * d[5] * d[8] * d[26];
    coeffs[250] =
        -2 * d[5] * d[9] * d[18] + 2 * d[3] * d[11] * d[18] + 2 * d[2] * d[12] * d[18] - 2 * d[0] * d[14] * d[18] -
        2 * d[5] * d[10] * d[19] + 2 * d[4] * d[11] * d[19] + 2 * d[2] * d[13] * d[19] - 2 * d[1] * d[14] * d[19] +
        2 * d[3] * d[9] * d[20] + 2 * d[4] * d[10] * d[20] + 2 * d[5] * d[11] * d[20] + 2 * d[0] * d[12] * d[20] +
        2 * d[1] * d[13] * d[20] + 2 * d[2] * d[14] * d[20] + 2 * d[2] * d[9] * d[21] + 2 * d[0] * d[11] * d[21] +
        2 * d[5] * d[12] * d[21] + 2 * d[3] * d[14] * d[21] + 2 * d[8] * d[15] * d[21] + 2 * d[6] * d[17] * d[21] +
        2 * d[2] * d[10] * d[22] + 2 * d[1] * d[11] * d[22] + 2 * d[5] * d[13] * d[22] + 2 * d[4] * d[14] * d[22] +
        2 * d[8] * d[16] * d[22] + 2 * d[7] * d[17] * d[22] - 2 * d[0] * d[9] * d[23] - 2 * d[1] * d[10] * d[23] +
        2 * d[2] * d[11] * d[23] + 2 * d[3] * d[12] * d[23] + 2 * d[4] * d[13] * d[23] + 6 * d[5] * d[14] * d[23] -
        2 * d[6] * d[15] * d[23] - 2 * d[7] * d[16] * d[23] + 2 * d[8] * d[17] * d[23] + 2 * d[8] * d[12] * d[24] -
        2 * d[6] * d[14] * d[24] - 2 * d[5] * d[15] * d[24] + 2 * d[3] * d[17] * d[24] + 2 * d[8] * d[13] * d[25] -
        2 * d[7] * d[14] * d[25] - 2 * d[5] * d[16] * d[25] + 2 * d[4] * d[17] * d[25] + 2 * d[6] * d[12] * d[26] +
        2 * d[7] * d[13] * d[26] + 2 * d[8] * d[14] * d[26] + 2 * d[3] * d[15] * d[26] + 2 * d[4] * d[16] * d[26] +
        2 * d[5] * d[17] * d[26];
    coeffs[251] =
        2 * d[11] * d[12] * d[18] - 2 * d[9] * d[14] * d[18] + 2 * d[11] * d[13] * d[19] - 2 * d[10] * d[14] * d[19] +
        2 * d[9] * d[12] * d[20] + 2 * d[10] * d[13] * d[20] + 2 * d[11] * d[14] * d[20] + 2 * d[9] * d[11] * d[21] +
        2 * d[12] * d[14] * d[21] + 2 * d[15] * d[17] * d[21] + 2 * d[10] * d[11] * d[22] + 2 * d[13] * d[14] * d[22] +
        2 * d[16] * d[17] * d[22] - std::pow(d[9], 2) * d[23] - std::pow(d[10], 2) * d[23] +
        std::pow(d[11], 2) * d[23] + std::pow(d[12], 2) * d[23] + std::pow(d[13], 2) * d[23] +
        3 * std::pow(d[14], 2) * d[23] - std::pow(d[15], 2) * d[23] - std::pow(d[16], 2) * d[23] +
        std::pow(d[17], 2) * d[23] - 2 * d[14] * d[15] * d[24] + 2 * d[12] * d[17] * d[24] - 2 * d[14] * d[16] * d[25] +
        2 * d[13] * d[17] * d[25] + 2 * d[12] * d[15] * d[26] + 2 * d[13] * d[16] * d[26] + 2 * d[14] * d[17] * d[26];
    coeffs[252] =
        -d[5] * std::pow(d[18], 2) - d[5] * std::pow(d[19], 2) + 2 * d[3] * d[18] * d[20] + 2 * d[4] * d[19] * d[20] +
        d[5] * std::pow(d[20], 2) + 2 * d[2] * d[18] * d[21] + 2 * d[0] * d[20] * d[21] + d[5] * std::pow(d[21], 2) +
        2 * d[2] * d[19] * d[22] + 2 * d[1] * d[20] * d[22] + d[5] * std::pow(d[22], 2) - 2 * d[0] * d[18] * d[23] -
        2 * d[1] * d[19] * d[23] + 2 * d[2] * d[20] * d[23] + 2 * d[3] * d[21] * d[23] + 2 * d[4] * d[22] * d[23] +
        3 * d[5] * std::pow(d[23], 2) + 2 * d[8] * d[21] * d[24] - 2 * d[6] * d[23] * d[24] -
        d[5] * std::pow(d[24], 2) + 2 * d[8] * d[22] * d[25] - 2 * d[7] * d[23] * d[25] - d[5] * std::pow(d[25], 2) +
        2 * d[6] * d[21] * d[26] + 2 * d[7] * d[22] * d[26] + 2 * d[8] * d[23] * d[26] + 2 * d[3] * d[24] * d[26] +
        2 * d[4] * d[25] * d[26] + d[5] * std::pow(d[26], 2);
    coeffs[253] =
        -d[14] * std::pow(d[18], 2) - d[14] * std::pow(d[19], 2) + 2 * d[12] * d[18] * d[20] +
        2 * d[13] * d[19] * d[20] + d[14] * std::pow(d[20], 2) + 2 * d[11] * d[18] * d[21] + 2 * d[9] * d[20] * d[21] +
        d[14] * std::pow(d[21], 2) + 2 * d[11] * d[19] * d[22] + 2 * d[10] * d[20] * d[22] +
        d[14] * std::pow(d[22], 2) - 2 * d[9] * d[18] * d[23] - 2 * d[10] * d[19] * d[23] + 2 * d[11] * d[20] * d[23] +
        2 * d[12] * d[21] * d[23] + 2 * d[13] * d[22] * d[23] + 3 * d[14] * std::pow(d[23], 2) +
        2 * d[17] * d[21] * d[24] - 2 * d[15] * d[23] * d[24] - d[14] * std::pow(d[24], 2) + 2 * d[17] * d[22] * d[25] -
        2 * d[16] * d[23] * d[25] - d[14] * std::pow(d[25], 2) + 2 * d[15] * d[21] * d[26] + 2 * d[16] * d[22] * d[26] +
        2 * d[17] * d[23] * d[26] + 2 * d[12] * d[24] * d[26] + 2 * d[13] * d[25] * d[26] + d[14] * std::pow(d[26], 2);
    coeffs[254] = 2 * d[18] * d[20] * d[21] + 2 * d[19] * d[20] * d[22] - std::pow(d[18], 2) * d[23] -
                  std::pow(d[19], 2) * d[23] + std::pow(d[20], 2) * d[23] + std::pow(d[21], 2) * d[23] +
                  std::pow(d[22], 2) * d[23] + std::pow(d[23], 3) - d[23] * std::pow(d[24], 2) -
                  d[23] * std::pow(d[25], 2) + 2 * d[21] * d[24] * d[26] + 2 * d[22] * d[25] * d[26] +
                  d[23] * std::pow(d[26], 2);
    coeffs[255] =
        2 * d[2] * d[3] * d[27] - 2 * d[0] * d[5] * d[27] + 2 * d[2] * d[4] * d[28] - 2 * d[1] * d[5] * d[28] +
        2 * d[0] * d[3] * d[29] + 2 * d[1] * d[4] * d[29] + 2 * d[2] * d[5] * d[29] + 2 * d[0] * d[2] * d[30] +
        2 * d[3] * d[5] * d[30] + 2 * d[6] * d[8] * d[30] + 2 * d[1] * d[2] * d[31] + 2 * d[4] * d[5] * d[31] +
        2 * d[7] * d[8] * d[31] - std::pow(d[0], 2) * d[32] - std::pow(d[1], 2) * d[32] + std::pow(d[2], 2) * d[32] +
        std::pow(d[3], 2) * d[32] + std::pow(d[4], 2) * d[32] + 3 * std::pow(d[5], 2) * d[32] -
        std::pow(d[6], 2) * d[32] - std::pow(d[7], 2) * d[32] + std::pow(d[8], 2) * d[32] - 2 * d[5] * d[6] * d[33] +
        2 * d[3] * d[8] * d[33] - 2 * d[5] * d[7] * d[34] + 2 * d[4] * d[8] * d[34] + 2 * d[3] * d[6] * d[35] +
        2 * d[4] * d[7] * d[35] + 2 * d[5] * d[8] * d[35];
    coeffs[256] =
        -2 * d[5] * d[9] * d[27] + 2 * d[3] * d[11] * d[27] + 2 * d[2] * d[12] * d[27] - 2 * d[0] * d[14] * d[27] -
        2 * d[5] * d[10] * d[28] + 2 * d[4] * d[11] * d[28] + 2 * d[2] * d[13] * d[28] - 2 * d[1] * d[14] * d[28] +
        2 * d[3] * d[9] * d[29] + 2 * d[4] * d[10] * d[29] + 2 * d[5] * d[11] * d[29] + 2 * d[0] * d[12] * d[29] +
        2 * d[1] * d[13] * d[29] + 2 * d[2] * d[14] * d[29] + 2 * d[2] * d[9] * d[30] + 2 * d[0] * d[11] * d[30] +
        2 * d[5] * d[12] * d[30] + 2 * d[3] * d[14] * d[30] + 2 * d[8] * d[15] * d[30] + 2 * d[6] * d[17] * d[30] +
        2 * d[2] * d[10] * d[31] + 2 * d[1] * d[11] * d[31] + 2 * d[5] * d[13] * d[31] + 2 * d[4] * d[14] * d[31] +
        2 * d[8] * d[16] * d[31] + 2 * d[7] * d[17] * d[31] - 2 * d[0] * d[9] * d[32] - 2 * d[1] * d[10] * d[32] +
        2 * d[2] * d[11] * d[32] + 2 * d[3] * d[12] * d[32] + 2 * d[4] * d[13] * d[32] + 6 * d[5] * d[14] * d[32] -
        2 * d[6] * d[15] * d[32] - 2 * d[7] * d[16] * d[32] + 2 * d[8] * d[17] * d[32] + 2 * d[8] * d[12] * d[33] -
        2 * d[6] * d[14] * d[33] - 2 * d[5] * d[15] * d[33] + 2 * d[3] * d[17] * d[33] + 2 * d[8] * d[13] * d[34] -
        2 * d[7] * d[14] * d[34] - 2 * d[5] * d[16] * d[34] + 2 * d[4] * d[17] * d[34] + 2 * d[6] * d[12] * d[35] +
        2 * d[7] * d[13] * d[35] + 2 * d[8] * d[14] * d[35] + 2 * d[3] * d[15] * d[35] + 2 * d[4] * d[16] * d[35] +
        2 * d[5] * d[17] * d[35];
    coeffs[257] =
        2 * d[11] * d[12] * d[27] - 2 * d[9] * d[14] * d[27] + 2 * d[11] * d[13] * d[28] - 2 * d[10] * d[14] * d[28] +
        2 * d[9] * d[12] * d[29] + 2 * d[10] * d[13] * d[29] + 2 * d[11] * d[14] * d[29] + 2 * d[9] * d[11] * d[30] +
        2 * d[12] * d[14] * d[30] + 2 * d[15] * d[17] * d[30] + 2 * d[10] * d[11] * d[31] + 2 * d[13] * d[14] * d[31] +
        2 * d[16] * d[17] * d[31] - std::pow(d[9], 2) * d[32] - std::pow(d[10], 2) * d[32] +
        std::pow(d[11], 2) * d[32] + std::pow(d[12], 2) * d[32] + std::pow(d[13], 2) * d[32] +
        3 * std::pow(d[14], 2) * d[32] - std::pow(d[15], 2) * d[32] - std::pow(d[16], 2) * d[32] +
        std::pow(d[17], 2) * d[32] - 2 * d[14] * d[15] * d[33] + 2 * d[12] * d[17] * d[33] - 2 * d[14] * d[16] * d[34] +
        2 * d[13] * d[17] * d[34] + 2 * d[12] * d[15] * d[35] + 2 * d[13] * d[16] * d[35] + 2 * d[14] * d[17] * d[35];
    coeffs[258] =
        -2 * d[5] * d[18] * d[27] + 2 * d[3] * d[20] * d[27] + 2 * d[2] * d[21] * d[27] - 2 * d[0] * d[23] * d[27] -
        2 * d[5] * d[19] * d[28] + 2 * d[4] * d[20] * d[28] + 2 * d[2] * d[22] * d[28] - 2 * d[1] * d[23] * d[28] +
        2 * d[3] * d[18] * d[29] + 2 * d[4] * d[19] * d[29] + 2 * d[5] * d[20] * d[29] + 2 * d[0] * d[21] * d[29] +
        2 * d[1] * d[22] * d[29] + 2 * d[2] * d[23] * d[29] + 2 * d[2] * d[18] * d[30] + 2 * d[0] * d[20] * d[30] +
        2 * d[5] * d[21] * d[30] + 2 * d[3] * d[23] * d[30] + 2 * d[8] * d[24] * d[30] + 2 * d[6] * d[26] * d[30] +
        2 * d[2] * d[19] * d[31] + 2 * d[1] * d[20] * d[31] + 2 * d[5] * d[22] * d[31] + 2 * d[4] * d[23] * d[31] +
        2 * d[8] * d[25] * d[31] + 2 * d[7] * d[26] * d[31] - 2 * d[0] * d[18] * d[32] - 2 * d[1] * d[19] * d[32] +
        2 * d[2] * d[20] * d[32] + 2 * d[3] * d[21] * d[32] + 2 * d[4] * d[22] * d[32] + 6 * d[5] * d[23] * d[32] -
        2 * d[6] * d[24] * d[32] - 2 * d[7] * d[25] * d[32] + 2 * d[8] * d[26] * d[32] + 2 * d[8] * d[21] * d[33] -
        2 * d[6] * d[23] * d[33] - 2 * d[5] * d[24] * d[33] + 2 * d[3] * d[26] * d[33] + 2 * d[8] * d[22] * d[34] -
        2 * d[7] * d[23] * d[34] - 2 * d[5] * d[25] * d[34] + 2 * d[4] * d[26] * d[34] + 2 * d[6] * d[21] * d[35] +
        2 * d[7] * d[22] * d[35] + 2 * d[8] * d[23] * d[35] + 2 * d[3] * d[24] * d[35] + 2 * d[4] * d[25] * d[35] +
        2 * d[5] * d[26] * d[35];
    coeffs[259] =
        -2 * d[14] * d[18] * d[27] + 2 * d[12] * d[20] * d[27] + 2 * d[11] * d[21] * d[27] - 2 * d[9] * d[23] * d[27] -
        2 * d[14] * d[19] * d[28] + 2 * d[13] * d[20] * d[28] + 2 * d[11] * d[22] * d[28] - 2 * d[10] * d[23] * d[28] +
        2 * d[12] * d[18] * d[29] + 2 * d[13] * d[19] * d[29] + 2 * d[14] * d[20] * d[29] + 2 * d[9] * d[21] * d[29] +
        2 * d[10] * d[22] * d[29] + 2 * d[11] * d[23] * d[29] + 2 * d[11] * d[18] * d[30] + 2 * d[9] * d[20] * d[30] +
        2 * d[14] * d[21] * d[30] + 2 * d[12] * d[23] * d[30] + 2 * d[17] * d[24] * d[30] + 2 * d[15] * d[26] * d[30] +
        2 * d[11] * d[19] * d[31] + 2 * d[10] * d[20] * d[31] + 2 * d[14] * d[22] * d[31] + 2 * d[13] * d[23] * d[31] +
        2 * d[17] * d[25] * d[31] + 2 * d[16] * d[26] * d[31] - 2 * d[9] * d[18] * d[32] - 2 * d[10] * d[19] * d[32] +
        2 * d[11] * d[20] * d[32] + 2 * d[12] * d[21] * d[32] + 2 * d[13] * d[22] * d[32] + 6 * d[14] * d[23] * d[32] -
        2 * d[15] * d[24] * d[32] - 2 * d[16] * d[25] * d[32] + 2 * d[17] * d[26] * d[32] + 2 * d[17] * d[21] * d[33] -
        2 * d[15] * d[23] * d[33] - 2 * d[14] * d[24] * d[33] + 2 * d[12] * d[26] * d[33] + 2 * d[17] * d[22] * d[34] -
        2 * d[16] * d[23] * d[34] - 2 * d[14] * d[25] * d[34] + 2 * d[13] * d[26] * d[34] + 2 * d[15] * d[21] * d[35] +
        2 * d[16] * d[22] * d[35] + 2 * d[17] * d[23] * d[35] + 2 * d[12] * d[24] * d[35] + 2 * d[13] * d[25] * d[35] +
        2 * d[14] * d[26] * d[35];
    coeffs[260] =
        2 * d[20] * d[21] * d[27] - 2 * d[18] * d[23] * d[27] + 2 * d[20] * d[22] * d[28] - 2 * d[19] * d[23] * d[28] +
        2 * d[18] * d[21] * d[29] + 2 * d[19] * d[22] * d[29] + 2 * d[20] * d[23] * d[29] + 2 * d[18] * d[20] * d[30] +
        2 * d[21] * d[23] * d[30] + 2 * d[24] * d[26] * d[30] + 2 * d[19] * d[20] * d[31] + 2 * d[22] * d[23] * d[31] +
        2 * d[25] * d[26] * d[31] - std::pow(d[18], 2) * d[32] - std::pow(d[19], 2) * d[32] +
        std::pow(d[20], 2) * d[32] + std::pow(d[21], 2) * d[32] + std::pow(d[22], 2) * d[32] +
        3 * std::pow(d[23], 2) * d[32] - std::pow(d[24], 2) * d[32] - std::pow(d[25], 2) * d[32] +
        std::pow(d[26], 2) * d[32] - 2 * d[23] * d[24] * d[33] + 2 * d[21] * d[26] * d[33] - 2 * d[23] * d[25] * d[34] +
        2 * d[22] * d[26] * d[34] + 2 * d[21] * d[24] * d[35] + 2 * d[22] * d[25] * d[35] + 2 * d[23] * d[26] * d[35];
    coeffs[261] =
        -d[5] * std::pow(d[27], 2) - d[5] * std::pow(d[28], 2) + 2 * d[3] * d[27] * d[29] + 2 * d[4] * d[28] * d[29] +
        d[5] * std::pow(d[29], 2) + 2 * d[2] * d[27] * d[30] + 2 * d[0] * d[29] * d[30] + d[5] * std::pow(d[30], 2) +
        2 * d[2] * d[28] * d[31] + 2 * d[1] * d[29] * d[31] + d[5] * std::pow(d[31], 2) - 2 * d[0] * d[27] * d[32] -
        2 * d[1] * d[28] * d[32] + 2 * d[2] * d[29] * d[32] + 2 * d[3] * d[30] * d[32] + 2 * d[4] * d[31] * d[32] +
        3 * d[5] * std::pow(d[32], 2) + 2 * d[8] * d[30] * d[33] - 2 * d[6] * d[32] * d[33] -
        d[5] * std::pow(d[33], 2) + 2 * d[8] * d[31] * d[34] - 2 * d[7] * d[32] * d[34] - d[5] * std::pow(d[34], 2) +
        2 * d[6] * d[30] * d[35] + 2 * d[7] * d[31] * d[35] + 2 * d[8] * d[32] * d[35] + 2 * d[3] * d[33] * d[35] +
        2 * d[4] * d[34] * d[35] + d[5] * std::pow(d[35], 2);
    coeffs[262] =
        -d[14] * std::pow(d[27], 2) - d[14] * std::pow(d[28], 2) + 2 * d[12] * d[27] * d[29] +
        2 * d[13] * d[28] * d[29] + d[14] * std::pow(d[29], 2) + 2 * d[11] * d[27] * d[30] + 2 * d[9] * d[29] * d[30] +
        d[14] * std::pow(d[30], 2) + 2 * d[11] * d[28] * d[31] + 2 * d[10] * d[29] * d[31] +
        d[14] * std::pow(d[31], 2) - 2 * d[9] * d[27] * d[32] - 2 * d[10] * d[28] * d[32] + 2 * d[11] * d[29] * d[32] +
        2 * d[12] * d[30] * d[32] + 2 * d[13] * d[31] * d[32] + 3 * d[14] * std::pow(d[32], 2) +
        2 * d[17] * d[30] * d[33] - 2 * d[15] * d[32] * d[33] - d[14] * std::pow(d[33], 2) + 2 * d[17] * d[31] * d[34] -
        2 * d[16] * d[32] * d[34] - d[14] * std::pow(d[34], 2) + 2 * d[15] * d[30] * d[35] + 2 * d[16] * d[31] * d[35] +
        2 * d[17] * d[32] * d[35] + 2 * d[12] * d[33] * d[35] + 2 * d[13] * d[34] * d[35] + d[14] * std::pow(d[35], 2);
    coeffs[263] =
        -d[23] * std::pow(d[27], 2) - d[23] * std::pow(d[28], 2) + 2 * d[21] * d[27] * d[29] +
        2 * d[22] * d[28] * d[29] + d[23] * std::pow(d[29], 2) + 2 * d[20] * d[27] * d[30] + 2 * d[18] * d[29] * d[30] +
        d[23] * std::pow(d[30], 2) + 2 * d[20] * d[28] * d[31] + 2 * d[19] * d[29] * d[31] +
        d[23] * std::pow(d[31], 2) - 2 * d[18] * d[27] * d[32] - 2 * d[19] * d[28] * d[32] + 2 * d[20] * d[29] * d[32] +
        2 * d[21] * d[30] * d[32] + 2 * d[22] * d[31] * d[32] + 3 * d[23] * std::pow(d[32], 2) +
        2 * d[26] * d[30] * d[33] - 2 * d[24] * d[32] * d[33] - d[23] * std::pow(d[33], 2) + 2 * d[26] * d[31] * d[34] -
        2 * d[25] * d[32] * d[34] - d[23] * std::pow(d[34], 2) + 2 * d[24] * d[30] * d[35] + 2 * d[25] * d[31] * d[35] +
        2 * d[26] * d[32] * d[35] + 2 * d[21] * d[33] * d[35] + 2 * d[22] * d[34] * d[35] + d[23] * std::pow(d[35], 2);
    coeffs[264] = 2 * d[27] * d[29] * d[30] + 2 * d[28] * d[29] * d[31] - std::pow(d[27], 2) * d[32] -
                  std::pow(d[28], 2) * d[32] + std::pow(d[29], 2) * d[32] + std::pow(d[30], 2) * d[32] +
                  std::pow(d[31], 2) * d[32] + std::pow(d[32], 3) - d[32] * std::pow(d[33], 2) -
                  d[32] * std::pow(d[34], 2) + 2 * d[30] * d[33] * d[35] + 2 * d[31] * d[34] * d[35] +
                  d[32] * std::pow(d[35], 2);
    coeffs[265] =
        2 * d[2] * d[3] * d[36] - 2 * d[0] * d[5] * d[36] + 2 * d[2] * d[4] * d[37] - 2 * d[1] * d[5] * d[37] +
        2 * d[0] * d[3] * d[38] + 2 * d[1] * d[4] * d[38] + 2 * d[2] * d[5] * d[38] + 2 * d[0] * d[2] * d[39] +
        2 * d[3] * d[5] * d[39] + 2 * d[6] * d[8] * d[39] + 2 * d[1] * d[2] * d[40] + 2 * d[4] * d[5] * d[40] +
        2 * d[7] * d[8] * d[40] - std::pow(d[0], 2) * d[41] - std::pow(d[1], 2) * d[41] + std::pow(d[2], 2) * d[41] +
        std::pow(d[3], 2) * d[41] + std::pow(d[4], 2) * d[41] + 3 * std::pow(d[5], 2) * d[41] -
        std::pow(d[6], 2) * d[41] - std::pow(d[7], 2) * d[41] + std::pow(d[8], 2) * d[41] - 2 * d[5] * d[6] * d[42] +
        2 * d[3] * d[8] * d[42] - 2 * d[5] * d[7] * d[43] + 2 * d[4] * d[8] * d[43] + 2 * d[3] * d[6] * d[44] +
        2 * d[4] * d[7] * d[44] + 2 * d[5] * d[8] * d[44];
    coeffs[266] =
        -2 * d[5] * d[9] * d[36] + 2 * d[3] * d[11] * d[36] + 2 * d[2] * d[12] * d[36] - 2 * d[0] * d[14] * d[36] -
        2 * d[5] * d[10] * d[37] + 2 * d[4] * d[11] * d[37] + 2 * d[2] * d[13] * d[37] - 2 * d[1] * d[14] * d[37] +
        2 * d[3] * d[9] * d[38] + 2 * d[4] * d[10] * d[38] + 2 * d[5] * d[11] * d[38] + 2 * d[0] * d[12] * d[38] +
        2 * d[1] * d[13] * d[38] + 2 * d[2] * d[14] * d[38] + 2 * d[2] * d[9] * d[39] + 2 * d[0] * d[11] * d[39] +
        2 * d[5] * d[12] * d[39] + 2 * d[3] * d[14] * d[39] + 2 * d[8] * d[15] * d[39] + 2 * d[6] * d[17] * d[39] +
        2 * d[2] * d[10] * d[40] + 2 * d[1] * d[11] * d[40] + 2 * d[5] * d[13] * d[40] + 2 * d[4] * d[14] * d[40] +
        2 * d[8] * d[16] * d[40] + 2 * d[7] * d[17] * d[40] - 2 * d[0] * d[9] * d[41] - 2 * d[1] * d[10] * d[41] +
        2 * d[2] * d[11] * d[41] + 2 * d[3] * d[12] * d[41] + 2 * d[4] * d[13] * d[41] + 6 * d[5] * d[14] * d[41] -
        2 * d[6] * d[15] * d[41] - 2 * d[7] * d[16] * d[41] + 2 * d[8] * d[17] * d[41] + 2 * d[8] * d[12] * d[42] -
        2 * d[6] * d[14] * d[42] - 2 * d[5] * d[15] * d[42] + 2 * d[3] * d[17] * d[42] + 2 * d[8] * d[13] * d[43] -
        2 * d[7] * d[14] * d[43] - 2 * d[5] * d[16] * d[43] + 2 * d[4] * d[17] * d[43] + 2 * d[6] * d[12] * d[44] +
        2 * d[7] * d[13] * d[44] + 2 * d[8] * d[14] * d[44] + 2 * d[3] * d[15] * d[44] + 2 * d[4] * d[16] * d[44] +
        2 * d[5] * d[17] * d[44];
    coeffs[267] =
        2 * d[11] * d[12] * d[36] - 2 * d[9] * d[14] * d[36] + 2 * d[11] * d[13] * d[37] - 2 * d[10] * d[14] * d[37] +
        2 * d[9] * d[12] * d[38] + 2 * d[10] * d[13] * d[38] + 2 * d[11] * d[14] * d[38] + 2 * d[9] * d[11] * d[39] +
        2 * d[12] * d[14] * d[39] + 2 * d[15] * d[17] * d[39] + 2 * d[10] * d[11] * d[40] + 2 * d[13] * d[14] * d[40] +
        2 * d[16] * d[17] * d[40] - std::pow(d[9], 2) * d[41] - std::pow(d[10], 2) * d[41] +
        std::pow(d[11], 2) * d[41] + std::pow(d[12], 2) * d[41] + std::pow(d[13], 2) * d[41] +
        3 * std::pow(d[14], 2) * d[41] - std::pow(d[15], 2) * d[41] - std::pow(d[16], 2) * d[41] +
        std::pow(d[17], 2) * d[41] - 2 * d[14] * d[15] * d[42] + 2 * d[12] * d[17] * d[42] - 2 * d[14] * d[16] * d[43] +
        2 * d[13] * d[17] * d[43] + 2 * d[12] * d[15] * d[44] + 2 * d[13] * d[16] * d[44] + 2 * d[14] * d[17] * d[44];
    coeffs[268] =
        -2 * d[5] * d[18] * d[36] + 2 * d[3] * d[20] * d[36] + 2 * d[2] * d[21] * d[36] - 2 * d[0] * d[23] * d[36] -
        2 * d[5] * d[19] * d[37] + 2 * d[4] * d[20] * d[37] + 2 * d[2] * d[22] * d[37] - 2 * d[1] * d[23] * d[37] +
        2 * d[3] * d[18] * d[38] + 2 * d[4] * d[19] * d[38] + 2 * d[5] * d[20] * d[38] + 2 * d[0] * d[21] * d[38] +
        2 * d[1] * d[22] * d[38] + 2 * d[2] * d[23] * d[38] + 2 * d[2] * d[18] * d[39] + 2 * d[0] * d[20] * d[39] +
        2 * d[5] * d[21] * d[39] + 2 * d[3] * d[23] * d[39] + 2 * d[8] * d[24] * d[39] + 2 * d[6] * d[26] * d[39] +
        2 * d[2] * d[19] * d[40] + 2 * d[1] * d[20] * d[40] + 2 * d[5] * d[22] * d[40] + 2 * d[4] * d[23] * d[40] +
        2 * d[8] * d[25] * d[40] + 2 * d[7] * d[26] * d[40] - 2 * d[0] * d[18] * d[41] - 2 * d[1] * d[19] * d[41] +
        2 * d[2] * d[20] * d[41] + 2 * d[3] * d[21] * d[41] + 2 * d[4] * d[22] * d[41] + 6 * d[5] * d[23] * d[41] -
        2 * d[6] * d[24] * d[41] - 2 * d[7] * d[25] * d[41] + 2 * d[8] * d[26] * d[41] + 2 * d[8] * d[21] * d[42] -
        2 * d[6] * d[23] * d[42] - 2 * d[5] * d[24] * d[42] + 2 * d[3] * d[26] * d[42] + 2 * d[8] * d[22] * d[43] -
        2 * d[7] * d[23] * d[43] - 2 * d[5] * d[25] * d[43] + 2 * d[4] * d[26] * d[43] + 2 * d[6] * d[21] * d[44] +
        2 * d[7] * d[22] * d[44] + 2 * d[8] * d[23] * d[44] + 2 * d[3] * d[24] * d[44] + 2 * d[4] * d[25] * d[44] +
        2 * d[5] * d[26] * d[44];
    coeffs[269] =
        -2 * d[14] * d[18] * d[36] + 2 * d[12] * d[20] * d[36] + 2 * d[11] * d[21] * d[36] - 2 * d[9] * d[23] * d[36] -
        2 * d[14] * d[19] * d[37] + 2 * d[13] * d[20] * d[37] + 2 * d[11] * d[22] * d[37] - 2 * d[10] * d[23] * d[37] +
        2 * d[12] * d[18] * d[38] + 2 * d[13] * d[19] * d[38] + 2 * d[14] * d[20] * d[38] + 2 * d[9] * d[21] * d[38] +
        2 * d[10] * d[22] * d[38] + 2 * d[11] * d[23] * d[38] + 2 * d[11] * d[18] * d[39] + 2 * d[9] * d[20] * d[39] +
        2 * d[14] * d[21] * d[39] + 2 * d[12] * d[23] * d[39] + 2 * d[17] * d[24] * d[39] + 2 * d[15] * d[26] * d[39] +
        2 * d[11] * d[19] * d[40] + 2 * d[10] * d[20] * d[40] + 2 * d[14] * d[22] * d[40] + 2 * d[13] * d[23] * d[40] +
        2 * d[17] * d[25] * d[40] + 2 * d[16] * d[26] * d[40] - 2 * d[9] * d[18] * d[41] - 2 * d[10] * d[19] * d[41] +
        2 * d[11] * d[20] * d[41] + 2 * d[12] * d[21] * d[41] + 2 * d[13] * d[22] * d[41] + 6 * d[14] * d[23] * d[41] -
        2 * d[15] * d[24] * d[41] - 2 * d[16] * d[25] * d[41] + 2 * d[17] * d[26] * d[41] + 2 * d[17] * d[21] * d[42] -
        2 * d[15] * d[23] * d[42] - 2 * d[14] * d[24] * d[42] + 2 * d[12] * d[26] * d[42] + 2 * d[17] * d[22] * d[43] -
        2 * d[16] * d[23] * d[43] - 2 * d[14] * d[25] * d[43] + 2 * d[13] * d[26] * d[43] + 2 * d[15] * d[21] * d[44] +
        2 * d[16] * d[22] * d[44] + 2 * d[17] * d[23] * d[44] + 2 * d[12] * d[24] * d[44] + 2 * d[13] * d[25] * d[44] +
        2 * d[14] * d[26] * d[44];
    coeffs[270] =
        2 * d[20] * d[21] * d[36] - 2 * d[18] * d[23] * d[36] + 2 * d[20] * d[22] * d[37] - 2 * d[19] * d[23] * d[37] +
        2 * d[18] * d[21] * d[38] + 2 * d[19] * d[22] * d[38] + 2 * d[20] * d[23] * d[38] + 2 * d[18] * d[20] * d[39] +
        2 * d[21] * d[23] * d[39] + 2 * d[24] * d[26] * d[39] + 2 * d[19] * d[20] * d[40] + 2 * d[22] * d[23] * d[40] +
        2 * d[25] * d[26] * d[40] - std::pow(d[18], 2) * d[41] - std::pow(d[19], 2) * d[41] +
        std::pow(d[20], 2) * d[41] + std::pow(d[21], 2) * d[41] + std::pow(d[22], 2) * d[41] +
        3 * std::pow(d[23], 2) * d[41] - std::pow(d[24], 2) * d[41] - std::pow(d[25], 2) * d[41] +
        std::pow(d[26], 2) * d[41] - 2 * d[23] * d[24] * d[42] + 2 * d[21] * d[26] * d[42] - 2 * d[23] * d[25] * d[43] +
        2 * d[22] * d[26] * d[43] + 2 * d[21] * d[24] * d[44] + 2 * d[22] * d[25] * d[44] + 2 * d[23] * d[26] * d[44];
    coeffs[271] =
        -2 * d[5] * d[27] * d[36] + 2 * d[3] * d[29] * d[36] + 2 * d[2] * d[30] * d[36] - 2 * d[0] * d[32] * d[36] -
        2 * d[5] * d[28] * d[37] + 2 * d[4] * d[29] * d[37] + 2 * d[2] * d[31] * d[37] - 2 * d[1] * d[32] * d[37] +
        2 * d[3] * d[27] * d[38] + 2 * d[4] * d[28] * d[38] + 2 * d[5] * d[29] * d[38] + 2 * d[0] * d[30] * d[38] +
        2 * d[1] * d[31] * d[38] + 2 * d[2] * d[32] * d[38] + 2 * d[2] * d[27] * d[39] + 2 * d[0] * d[29] * d[39] +
        2 * d[5] * d[30] * d[39] + 2 * d[3] * d[32] * d[39] + 2 * d[8] * d[33] * d[39] + 2 * d[6] * d[35] * d[39] +
        2 * d[2] * d[28] * d[40] + 2 * d[1] * d[29] * d[40] + 2 * d[5] * d[31] * d[40] + 2 * d[4] * d[32] * d[40] +
        2 * d[8] * d[34] * d[40] + 2 * d[7] * d[35] * d[40] - 2 * d[0] * d[27] * d[41] - 2 * d[1] * d[28] * d[41] +
        2 * d[2] * d[29] * d[41] + 2 * d[3] * d[30] * d[41] + 2 * d[4] * d[31] * d[41] + 6 * d[5] * d[32] * d[41] -
        2 * d[6] * d[33] * d[41] - 2 * d[7] * d[34] * d[41] + 2 * d[8] * d[35] * d[41] + 2 * d[8] * d[30] * d[42] -
        2 * d[6] * d[32] * d[42] - 2 * d[5] * d[33] * d[42] + 2 * d[3] * d[35] * d[42] + 2 * d[8] * d[31] * d[43] -
        2 * d[7] * d[32] * d[43] - 2 * d[5] * d[34] * d[43] + 2 * d[4] * d[35] * d[43] + 2 * d[6] * d[30] * d[44] +
        2 * d[7] * d[31] * d[44] + 2 * d[8] * d[32] * d[44] + 2 * d[3] * d[33] * d[44] + 2 * d[4] * d[34] * d[44] +
        2 * d[5] * d[35] * d[44];
    coeffs[272] =
        -2 * d[14] * d[27] * d[36] + 2 * d[12] * d[29] * d[36] + 2 * d[11] * d[30] * d[36] - 2 * d[9] * d[32] * d[36] -
        2 * d[14] * d[28] * d[37] + 2 * d[13] * d[29] * d[37] + 2 * d[11] * d[31] * d[37] - 2 * d[10] * d[32] * d[37] +
        2 * d[12] * d[27] * d[38] + 2 * d[13] * d[28] * d[38] + 2 * d[14] * d[29] * d[38] + 2 * d[9] * d[30] * d[38] +
        2 * d[10] * d[31] * d[38] + 2 * d[11] * d[32] * d[38] + 2 * d[11] * d[27] * d[39] + 2 * d[9] * d[29] * d[39] +
        2 * d[14] * d[30] * d[39] + 2 * d[12] * d[32] * d[39] + 2 * d[17] * d[33] * d[39] + 2 * d[15] * d[35] * d[39] +
        2 * d[11] * d[28] * d[40] + 2 * d[10] * d[29] * d[40] + 2 * d[14] * d[31] * d[40] + 2 * d[13] * d[32] * d[40] +
        2 * d[17] * d[34] * d[40] + 2 * d[16] * d[35] * d[40] - 2 * d[9] * d[27] * d[41] - 2 * d[10] * d[28] * d[41] +
        2 * d[11] * d[29] * d[41] + 2 * d[12] * d[30] * d[41] + 2 * d[13] * d[31] * d[41] + 6 * d[14] * d[32] * d[41] -
        2 * d[15] * d[33] * d[41] - 2 * d[16] * d[34] * d[41] + 2 * d[17] * d[35] * d[41] + 2 * d[17] * d[30] * d[42] -
        2 * d[15] * d[32] * d[42] - 2 * d[14] * d[33] * d[42] + 2 * d[12] * d[35] * d[42] + 2 * d[17] * d[31] * d[43] -
        2 * d[16] * d[32] * d[43] - 2 * d[14] * d[34] * d[43] + 2 * d[13] * d[35] * d[43] + 2 * d[15] * d[30] * d[44] +
        2 * d[16] * d[31] * d[44] + 2 * d[17] * d[32] * d[44] + 2 * d[12] * d[33] * d[44] + 2 * d[13] * d[34] * d[44] +
        2 * d[14] * d[35] * d[44];
    coeffs[273] =
        -2 * d[23] * d[27] * d[36] + 2 * d[21] * d[29] * d[36] + 2 * d[20] * d[30] * d[36] - 2 * d[18] * d[32] * d[36] -
        2 * d[23] * d[28] * d[37] + 2 * d[22] * d[29] * d[37] + 2 * d[20] * d[31] * d[37] - 2 * d[19] * d[32] * d[37] +
        2 * d[21] * d[27] * d[38] + 2 * d[22] * d[28] * d[38] + 2 * d[23] * d[29] * d[38] + 2 * d[18] * d[30] * d[38] +
        2 * d[19] * d[31] * d[38] + 2 * d[20] * d[32] * d[38] + 2 * d[20] * d[27] * d[39] + 2 * d[18] * d[29] * d[39] +
        2 * d[23] * d[30] * d[39] + 2 * d[21] * d[32] * d[39] + 2 * d[26] * d[33] * d[39] + 2 * d[24] * d[35] * d[39] +
        2 * d[20] * d[28] * d[40] + 2 * d[19] * d[29] * d[40] + 2 * d[23] * d[31] * d[40] + 2 * d[22] * d[32] * d[40] +
        2 * d[26] * d[34] * d[40] + 2 * d[25] * d[35] * d[40] - 2 * d[18] * d[27] * d[41] - 2 * d[19] * d[28] * d[41] +
        2 * d[20] * d[29] * d[41] + 2 * d[21] * d[30] * d[41] + 2 * d[22] * d[31] * d[41] + 6 * d[23] * d[32] * d[41] -
        2 * d[24] * d[33] * d[41] - 2 * d[25] * d[34] * d[41] + 2 * d[26] * d[35] * d[41] + 2 * d[26] * d[30] * d[42] -
        2 * d[24] * d[32] * d[42] - 2 * d[23] * d[33] * d[42] + 2 * d[21] * d[35] * d[42] + 2 * d[26] * d[31] * d[43] -
        2 * d[25] * d[32] * d[43] - 2 * d[23] * d[34] * d[43] + 2 * d[22] * d[35] * d[43] + 2 * d[24] * d[30] * d[44] +
        2 * d[25] * d[31] * d[44] + 2 * d[26] * d[32] * d[44] + 2 * d[21] * d[33] * d[44] + 2 * d[22] * d[34] * d[44] +
        2 * d[23] * d[35] * d[44];
    coeffs[274] =
        2 * d[29] * d[30] * d[36] - 2 * d[27] * d[32] * d[36] + 2 * d[29] * d[31] * d[37] - 2 * d[28] * d[32] * d[37] +
        2 * d[27] * d[30] * d[38] + 2 * d[28] * d[31] * d[38] + 2 * d[29] * d[32] * d[38] + 2 * d[27] * d[29] * d[39] +
        2 * d[30] * d[32] * d[39] + 2 * d[33] * d[35] * d[39] + 2 * d[28] * d[29] * d[40] + 2 * d[31] * d[32] * d[40] +
        2 * d[34] * d[35] * d[40] - std::pow(d[27], 2) * d[41] - std::pow(d[28], 2) * d[41] +
        std::pow(d[29], 2) * d[41] + std::pow(d[30], 2) * d[41] + std::pow(d[31], 2) * d[41] +
        3 * std::pow(d[32], 2) * d[41] - std::pow(d[33], 2) * d[41] - std::pow(d[34], 2) * d[41] +
        std::pow(d[35], 2) * d[41] - 2 * d[32] * d[33] * d[42] + 2 * d[30] * d[35] * d[42] - 2 * d[32] * d[34] * d[43] +
        2 * d[31] * d[35] * d[43] + 2 * d[30] * d[33] * d[44] + 2 * d[31] * d[34] * d[44] + 2 * d[32] * d[35] * d[44];
    coeffs[275] =
        -d[5] * std::pow(d[36], 2) - d[5] * std::pow(d[37], 2) + 2 * d[3] * d[36] * d[38] + 2 * d[4] * d[37] * d[38] +
        d[5] * std::pow(d[38], 2) + 2 * d[2] * d[36] * d[39] + 2 * d[0] * d[38] * d[39] + d[5] * std::pow(d[39], 2) +
        2 * d[2] * d[37] * d[40] + 2 * d[1] * d[38] * d[40] + d[5] * std::pow(d[40], 2) - 2 * d[0] * d[36] * d[41] -
        2 * d[1] * d[37] * d[41] + 2 * d[2] * d[38] * d[41] + 2 * d[3] * d[39] * d[41] + 2 * d[4] * d[40] * d[41] +
        3 * d[5] * std::pow(d[41], 2) + 2 * d[8] * d[39] * d[42] - 2 * d[6] * d[41] * d[42] -
        d[5] * std::pow(d[42], 2) + 2 * d[8] * d[40] * d[43] - 2 * d[7] * d[41] * d[43] - d[5] * std::pow(d[43], 2) +
        2 * d[6] * d[39] * d[44] + 2 * d[7] * d[40] * d[44] + 2 * d[8] * d[41] * d[44] + 2 * d[3] * d[42] * d[44] +
        2 * d[4] * d[43] * d[44] + d[5] * std::pow(d[44], 2);
    coeffs[276] =
        -d[14] * std::pow(d[36], 2) - d[14] * std::pow(d[37], 2) + 2 * d[12] * d[36] * d[38] +
        2 * d[13] * d[37] * d[38] + d[14] * std::pow(d[38], 2) + 2 * d[11] * d[36] * d[39] + 2 * d[9] * d[38] * d[39] +
        d[14] * std::pow(d[39], 2) + 2 * d[11] * d[37] * d[40] + 2 * d[10] * d[38] * d[40] +
        d[14] * std::pow(d[40], 2) - 2 * d[9] * d[36] * d[41] - 2 * d[10] * d[37] * d[41] + 2 * d[11] * d[38] * d[41] +
        2 * d[12] * d[39] * d[41] + 2 * d[13] * d[40] * d[41] + 3 * d[14] * std::pow(d[41], 2) +
        2 * d[17] * d[39] * d[42] - 2 * d[15] * d[41] * d[42] - d[14] * std::pow(d[42], 2) + 2 * d[17] * d[40] * d[43] -
        2 * d[16] * d[41] * d[43] - d[14] * std::pow(d[43], 2) + 2 * d[15] * d[39] * d[44] + 2 * d[16] * d[40] * d[44] +
        2 * d[17] * d[41] * d[44] + 2 * d[12] * d[42] * d[44] + 2 * d[13] * d[43] * d[44] + d[14] * std::pow(d[44], 2);
    coeffs[277] =
        -d[23] * std::pow(d[36], 2) - d[23] * std::pow(d[37], 2) + 2 * d[21] * d[36] * d[38] +
        2 * d[22] * d[37] * d[38] + d[23] * std::pow(d[38], 2) + 2 * d[20] * d[36] * d[39] + 2 * d[18] * d[38] * d[39] +
        d[23] * std::pow(d[39], 2) + 2 * d[20] * d[37] * d[40] + 2 * d[19] * d[38] * d[40] +
        d[23] * std::pow(d[40], 2) - 2 * d[18] * d[36] * d[41] - 2 * d[19] * d[37] * d[41] + 2 * d[20] * d[38] * d[41] +
        2 * d[21] * d[39] * d[41] + 2 * d[22] * d[40] * d[41] + 3 * d[23] * std::pow(d[41], 2) +
        2 * d[26] * d[39] * d[42] - 2 * d[24] * d[41] * d[42] - d[23] * std::pow(d[42], 2) + 2 * d[26] * d[40] * d[43] -
        2 * d[25] * d[41] * d[43] - d[23] * std::pow(d[43], 2) + 2 * d[24] * d[39] * d[44] + 2 * d[25] * d[40] * d[44] +
        2 * d[26] * d[41] * d[44] + 2 * d[21] * d[42] * d[44] + 2 * d[22] * d[43] * d[44] + d[23] * std::pow(d[44], 2);
    coeffs[278] =
        -d[32] * std::pow(d[36], 2) - d[32] * std::pow(d[37], 2) + 2 * d[30] * d[36] * d[38] +
        2 * d[31] * d[37] * d[38] + d[32] * std::pow(d[38], 2) + 2 * d[29] * d[36] * d[39] + 2 * d[27] * d[38] * d[39] +
        d[32] * std::pow(d[39], 2) + 2 * d[29] * d[37] * d[40] + 2 * d[28] * d[38] * d[40] +
        d[32] * std::pow(d[40], 2) - 2 * d[27] * d[36] * d[41] - 2 * d[28] * d[37] * d[41] + 2 * d[29] * d[38] * d[41] +
        2 * d[30] * d[39] * d[41] + 2 * d[31] * d[40] * d[41] + 3 * d[32] * std::pow(d[41], 2) +
        2 * d[35] * d[39] * d[42] - 2 * d[33] * d[41] * d[42] - d[32] * std::pow(d[42], 2) + 2 * d[35] * d[40] * d[43] -
        2 * d[34] * d[41] * d[43] - d[32] * std::pow(d[43], 2) + 2 * d[33] * d[39] * d[44] + 2 * d[34] * d[40] * d[44] +
        2 * d[35] * d[41] * d[44] + 2 * d[30] * d[42] * d[44] + 2 * d[31] * d[43] * d[44] + d[32] * std::pow(d[44], 2);
    coeffs[279] = 2 * d[36] * d[38] * d[39] + 2 * d[37] * d[38] * d[40] - std::pow(d[36], 2) * d[41] -
                  std::pow(d[37], 2) * d[41] + std::pow(d[38], 2) * d[41] + std::pow(d[39], 2) * d[41] +
                  std::pow(d[40], 2) * d[41] + std::pow(d[41], 3) - d[41] * std::pow(d[42], 2) -
                  d[41] * std::pow(d[43], 2) + 2 * d[39] * d[42] * d[44] + 2 * d[40] * d[43] * d[44] +
                  d[41] * std::pow(d[44], 2);
    coeffs[280] = std::pow(d[0], 2) * d[6] - std::pow(d[1], 2) * d[6] - std::pow(d[2], 2) * d[6] +
                  std::pow(d[3], 2) * d[6] - std::pow(d[4], 2) * d[6] - std::pow(d[5], 2) * d[6] + std::pow(d[6], 3) +
                  2 * d[0] * d[1] * d[7] + 2 * d[3] * d[4] * d[7] + d[6] * std::pow(d[7], 2) + 2 * d[0] * d[2] * d[8] +
                  2 * d[3] * d[5] * d[8] + d[6] * std::pow(d[8], 2);
    coeffs[281] = 2 * d[0] * d[6] * d[9] + 2 * d[1] * d[7] * d[9] + 2 * d[2] * d[8] * d[9] - 2 * d[1] * d[6] * d[10] +
                  2 * d[0] * d[7] * d[10] - 2 * d[2] * d[6] * d[11] + 2 * d[0] * d[8] * d[11] +
                  2 * d[3] * d[6] * d[12] + 2 * d[4] * d[7] * d[12] + 2 * d[5] * d[8] * d[12] -
                  2 * d[4] * d[6] * d[13] + 2 * d[3] * d[7] * d[13] - 2 * d[5] * d[6] * d[14] +
                  2 * d[3] * d[8] * d[14] + std::pow(d[0], 2) * d[15] - std::pow(d[1], 2) * d[15] -
                  std::pow(d[2], 2) * d[15] + std::pow(d[3], 2) * d[15] - std::pow(d[4], 2) * d[15] -
                  std::pow(d[5], 2) * d[15] + 3 * std::pow(d[6], 2) * d[15] + std::pow(d[7], 2) * d[15] +
                  std::pow(d[8], 2) * d[15] + 2 * d[0] * d[1] * d[16] + 2 * d[3] * d[4] * d[16] +
                  2 * d[6] * d[7] * d[16] + 2 * d[0] * d[2] * d[17] + 2 * d[3] * d[5] * d[17] + 2 * d[6] * d[8] * d[17];
    coeffs[282] =
        d[6] * std::pow(d[9], 2) + 2 * d[7] * d[9] * d[10] - d[6] * std::pow(d[10], 2) + 2 * d[8] * d[9] * d[11] -
        d[6] * std::pow(d[11], 2) + d[6] * std::pow(d[12], 2) + 2 * d[7] * d[12] * d[13] - d[6] * std::pow(d[13], 2) +
        2 * d[8] * d[12] * d[14] - d[6] * std::pow(d[14], 2) + 2 * d[0] * d[9] * d[15] - 2 * d[1] * d[10] * d[15] -
        2 * d[2] * d[11] * d[15] + 2 * d[3] * d[12] * d[15] - 2 * d[4] * d[13] * d[15] - 2 * d[5] * d[14] * d[15] +
        3 * d[6] * std::pow(d[15], 2) + 2 * d[1] * d[9] * d[16] + 2 * d[0] * d[10] * d[16] + 2 * d[4] * d[12] * d[16] +
        2 * d[3] * d[13] * d[16] + 2 * d[7] * d[15] * d[16] + d[6] * std::pow(d[16], 2) + 2 * d[2] * d[9] * d[17] +
        2 * d[0] * d[11] * d[17] + 2 * d[5] * d[12] * d[17] + 2 * d[3] * d[14] * d[17] + 2 * d[8] * d[15] * d[17] +
        d[6] * std::pow(d[17], 2);
    coeffs[283] = std::pow(d[9], 2) * d[15] - std::pow(d[10], 2) * d[15] - std::pow(d[11], 2) * d[15] +
                  std::pow(d[12], 2) * d[15] - std::pow(d[13], 2) * d[15] - std::pow(d[14], 2) * d[15] +
                  std::pow(d[15], 3) + 2 * d[9] * d[10] * d[16] + 2 * d[12] * d[13] * d[16] +
                  d[15] * std::pow(d[16], 2) + 2 * d[9] * d[11] * d[17] + 2 * d[12] * d[14] * d[17] +
                  d[15] * std::pow(d[17], 2);
    coeffs[284] =
        2 * d[0] * d[6] * d[18] + 2 * d[1] * d[7] * d[18] + 2 * d[2] * d[8] * d[18] - 2 * d[1] * d[6] * d[19] +
        2 * d[0] * d[7] * d[19] - 2 * d[2] * d[6] * d[20] + 2 * d[0] * d[8] * d[20] + 2 * d[3] * d[6] * d[21] +
        2 * d[4] * d[7] * d[21] + 2 * d[5] * d[8] * d[21] - 2 * d[4] * d[6] * d[22] + 2 * d[3] * d[7] * d[22] -
        2 * d[5] * d[6] * d[23] + 2 * d[3] * d[8] * d[23] + std::pow(d[0], 2) * d[24] - std::pow(d[1], 2) * d[24] -
        std::pow(d[2], 2) * d[24] + std::pow(d[3], 2) * d[24] - std::pow(d[4], 2) * d[24] - std::pow(d[5], 2) * d[24] +
        3 * std::pow(d[6], 2) * d[24] + std::pow(d[7], 2) * d[24] + std::pow(d[8], 2) * d[24] +
        2 * d[0] * d[1] * d[25] + 2 * d[3] * d[4] * d[25] + 2 * d[6] * d[7] * d[25] + 2 * d[0] * d[2] * d[26] +
        2 * d[3] * d[5] * d[26] + 2 * d[6] * d[8] * d[26];
    coeffs[285] =
        2 * d[6] * d[9] * d[18] + 2 * d[7] * d[10] * d[18] + 2 * d[8] * d[11] * d[18] + 2 * d[0] * d[15] * d[18] +
        2 * d[1] * d[16] * d[18] + 2 * d[2] * d[17] * d[18] + 2 * d[7] * d[9] * d[19] - 2 * d[6] * d[10] * d[19] -
        2 * d[1] * d[15] * d[19] + 2 * d[0] * d[16] * d[19] + 2 * d[8] * d[9] * d[20] - 2 * d[6] * d[11] * d[20] -
        2 * d[2] * d[15] * d[20] + 2 * d[0] * d[17] * d[20] + 2 * d[6] * d[12] * d[21] + 2 * d[7] * d[13] * d[21] +
        2 * d[8] * d[14] * d[21] + 2 * d[3] * d[15] * d[21] + 2 * d[4] * d[16] * d[21] + 2 * d[5] * d[17] * d[21] +
        2 * d[7] * d[12] * d[22] - 2 * d[6] * d[13] * d[22] - 2 * d[4] * d[15] * d[22] + 2 * d[3] * d[16] * d[22] +
        2 * d[8] * d[12] * d[23] - 2 * d[6] * d[14] * d[23] - 2 * d[5] * d[15] * d[23] + 2 * d[3] * d[17] * d[23] +
        2 * d[0] * d[9] * d[24] - 2 * d[1] * d[10] * d[24] - 2 * d[2] * d[11] * d[24] + 2 * d[3] * d[12] * d[24] -
        2 * d[4] * d[13] * d[24] - 2 * d[5] * d[14] * d[24] + 6 * d[6] * d[15] * d[24] + 2 * d[7] * d[16] * d[24] +
        2 * d[8] * d[17] * d[24] + 2 * d[1] * d[9] * d[25] + 2 * d[0] * d[10] * d[25] + 2 * d[4] * d[12] * d[25] +
        2 * d[3] * d[13] * d[25] + 2 * d[7] * d[15] * d[25] + 2 * d[6] * d[16] * d[25] + 2 * d[2] * d[9] * d[26] +
        2 * d[0] * d[11] * d[26] + 2 * d[5] * d[12] * d[26] + 2 * d[3] * d[14] * d[26] + 2 * d[8] * d[15] * d[26] +
        2 * d[6] * d[17] * d[26];
    coeffs[286] =
        2 * d[9] * d[15] * d[18] + 2 * d[10] * d[16] * d[18] + 2 * d[11] * d[17] * d[18] - 2 * d[10] * d[15] * d[19] +
        2 * d[9] * d[16] * d[19] - 2 * d[11] * d[15] * d[20] + 2 * d[9] * d[17] * d[20] + 2 * d[12] * d[15] * d[21] +
        2 * d[13] * d[16] * d[21] + 2 * d[14] * d[17] * d[21] - 2 * d[13] * d[15] * d[22] + 2 * d[12] * d[16] * d[22] -
        2 * d[14] * d[15] * d[23] + 2 * d[12] * d[17] * d[23] + std::pow(d[9], 2) * d[24] - std::pow(d[10], 2) * d[24] -
        std::pow(d[11], 2) * d[24] + std::pow(d[12], 2) * d[24] - std::pow(d[13], 2) * d[24] -
        std::pow(d[14], 2) * d[24] + 3 * std::pow(d[15], 2) * d[24] + std::pow(d[16], 2) * d[24] +
        std::pow(d[17], 2) * d[24] + 2 * d[9] * d[10] * d[25] + 2 * d[12] * d[13] * d[25] + 2 * d[15] * d[16] * d[25] +
        2 * d[9] * d[11] * d[26] + 2 * d[12] * d[14] * d[26] + 2 * d[15] * d[17] * d[26];
    coeffs[287] =
        d[6] * std::pow(d[18], 2) + 2 * d[7] * d[18] * d[19] - d[6] * std::pow(d[19], 2) + 2 * d[8] * d[18] * d[20] -
        d[6] * std::pow(d[20], 2) + d[6] * std::pow(d[21], 2) + 2 * d[7] * d[21] * d[22] - d[6] * std::pow(d[22], 2) +
        2 * d[8] * d[21] * d[23] - d[6] * std::pow(d[23], 2) + 2 * d[0] * d[18] * d[24] - 2 * d[1] * d[19] * d[24] -
        2 * d[2] * d[20] * d[24] + 2 * d[3] * d[21] * d[24] - 2 * d[4] * d[22] * d[24] - 2 * d[5] * d[23] * d[24] +
        3 * d[6] * std::pow(d[24], 2) + 2 * d[1] * d[18] * d[25] + 2 * d[0] * d[19] * d[25] + 2 * d[4] * d[21] * d[25] +
        2 * d[3] * d[22] * d[25] + 2 * d[7] * d[24] * d[25] + d[6] * std::pow(d[25], 2) + 2 * d[2] * d[18] * d[26] +
        2 * d[0] * d[20] * d[26] + 2 * d[5] * d[21] * d[26] + 2 * d[3] * d[23] * d[26] + 2 * d[8] * d[24] * d[26] +
        d[6] * std::pow(d[26], 2);
    coeffs[288] = d[15] * std::pow(d[18], 2) + 2 * d[16] * d[18] * d[19] - d[15] * std::pow(d[19], 2) +
                  2 * d[17] * d[18] * d[20] - d[15] * std::pow(d[20], 2) + d[15] * std::pow(d[21], 2) +
                  2 * d[16] * d[21] * d[22] - d[15] * std::pow(d[22], 2) + 2 * d[17] * d[21] * d[23] -
                  d[15] * std::pow(d[23], 2) + 2 * d[9] * d[18] * d[24] - 2 * d[10] * d[19] * d[24] -
                  2 * d[11] * d[20] * d[24] + 2 * d[12] * d[21] * d[24] - 2 * d[13] * d[22] * d[24] -
                  2 * d[14] * d[23] * d[24] + 3 * d[15] * std::pow(d[24], 2) + 2 * d[10] * d[18] * d[25] +
                  2 * d[9] * d[19] * d[25] + 2 * d[13] * d[21] * d[25] + 2 * d[12] * d[22] * d[25] +
                  2 * d[16] * d[24] * d[25] + d[15] * std::pow(d[25], 2) + 2 * d[11] * d[18] * d[26] +
                  2 * d[9] * d[20] * d[26] + 2 * d[14] * d[21] * d[26] + 2 * d[12] * d[23] * d[26] +
                  2 * d[17] * d[24] * d[26] + d[15] * std::pow(d[26], 2);
    coeffs[289] = std::pow(d[18], 2) * d[24] - std::pow(d[19], 2) * d[24] - std::pow(d[20], 2) * d[24] +
                  std::pow(d[21], 2) * d[24] - std::pow(d[22], 2) * d[24] - std::pow(d[23], 2) * d[24] +
                  std::pow(d[24], 3) + 2 * d[18] * d[19] * d[25] + 2 * d[21] * d[22] * d[25] +
                  d[24] * std::pow(d[25], 2) + 2 * d[18] * d[20] * d[26] + 2 * d[21] * d[23] * d[26] +
                  d[24] * std::pow(d[26], 2);
    coeffs[290] =
        2 * d[0] * d[6] * d[27] + 2 * d[1] * d[7] * d[27] + 2 * d[2] * d[8] * d[27] - 2 * d[1] * d[6] * d[28] +
        2 * d[0] * d[7] * d[28] - 2 * d[2] * d[6] * d[29] + 2 * d[0] * d[8] * d[29] + 2 * d[3] * d[6] * d[30] +
        2 * d[4] * d[7] * d[30] + 2 * d[5] * d[8] * d[30] - 2 * d[4] * d[6] * d[31] + 2 * d[3] * d[7] * d[31] -
        2 * d[5] * d[6] * d[32] + 2 * d[3] * d[8] * d[32] + std::pow(d[0], 2) * d[33] - std::pow(d[1], 2) * d[33] -
        std::pow(d[2], 2) * d[33] + std::pow(d[3], 2) * d[33] - std::pow(d[4], 2) * d[33] - std::pow(d[5], 2) * d[33] +
        3 * std::pow(d[6], 2) * d[33] + std::pow(d[7], 2) * d[33] + std::pow(d[8], 2) * d[33] +
        2 * d[0] * d[1] * d[34] + 2 * d[3] * d[4] * d[34] + 2 * d[6] * d[7] * d[34] + 2 * d[0] * d[2] * d[35] +
        2 * d[3] * d[5] * d[35] + 2 * d[6] * d[8] * d[35];
    coeffs[291] =
        2 * d[6] * d[9] * d[27] + 2 * d[7] * d[10] * d[27] + 2 * d[8] * d[11] * d[27] + 2 * d[0] * d[15] * d[27] +
        2 * d[1] * d[16] * d[27] + 2 * d[2] * d[17] * d[27] + 2 * d[7] * d[9] * d[28] - 2 * d[6] * d[10] * d[28] -
        2 * d[1] * d[15] * d[28] + 2 * d[0] * d[16] * d[28] + 2 * d[8] * d[9] * d[29] - 2 * d[6] * d[11] * d[29] -
        2 * d[2] * d[15] * d[29] + 2 * d[0] * d[17] * d[29] + 2 * d[6] * d[12] * d[30] + 2 * d[7] * d[13] * d[30] +
        2 * d[8] * d[14] * d[30] + 2 * d[3] * d[15] * d[30] + 2 * d[4] * d[16] * d[30] + 2 * d[5] * d[17] * d[30] +
        2 * d[7] * d[12] * d[31] - 2 * d[6] * d[13] * d[31] - 2 * d[4] * d[15] * d[31] + 2 * d[3] * d[16] * d[31] +
        2 * d[8] * d[12] * d[32] - 2 * d[6] * d[14] * d[32] - 2 * d[5] * d[15] * d[32] + 2 * d[3] * d[17] * d[32] +
        2 * d[0] * d[9] * d[33] - 2 * d[1] * d[10] * d[33] - 2 * d[2] * d[11] * d[33] + 2 * d[3] * d[12] * d[33] -
        2 * d[4] * d[13] * d[33] - 2 * d[5] * d[14] * d[33] + 6 * d[6] * d[15] * d[33] + 2 * d[7] * d[16] * d[33] +
        2 * d[8] * d[17] * d[33] + 2 * d[1] * d[9] * d[34] + 2 * d[0] * d[10] * d[34] + 2 * d[4] * d[12] * d[34] +
        2 * d[3] * d[13] * d[34] + 2 * d[7] * d[15] * d[34] + 2 * d[6] * d[16] * d[34] + 2 * d[2] * d[9] * d[35] +
        2 * d[0] * d[11] * d[35] + 2 * d[5] * d[12] * d[35] + 2 * d[3] * d[14] * d[35] + 2 * d[8] * d[15] * d[35] +
        2 * d[6] * d[17] * d[35];
    coeffs[292] =
        2 * d[9] * d[15] * d[27] + 2 * d[10] * d[16] * d[27] + 2 * d[11] * d[17] * d[27] - 2 * d[10] * d[15] * d[28] +
        2 * d[9] * d[16] * d[28] - 2 * d[11] * d[15] * d[29] + 2 * d[9] * d[17] * d[29] + 2 * d[12] * d[15] * d[30] +
        2 * d[13] * d[16] * d[30] + 2 * d[14] * d[17] * d[30] - 2 * d[13] * d[15] * d[31] + 2 * d[12] * d[16] * d[31] -
        2 * d[14] * d[15] * d[32] + 2 * d[12] * d[17] * d[32] + std::pow(d[9], 2) * d[33] - std::pow(d[10], 2) * d[33] -
        std::pow(d[11], 2) * d[33] + std::pow(d[12], 2) * d[33] - std::pow(d[13], 2) * d[33] -
        std::pow(d[14], 2) * d[33] + 3 * std::pow(d[15], 2) * d[33] + std::pow(d[16], 2) * d[33] +
        std::pow(d[17], 2) * d[33] + 2 * d[9] * d[10] * d[34] + 2 * d[12] * d[13] * d[34] + 2 * d[15] * d[16] * d[34] +
        2 * d[9] * d[11] * d[35] + 2 * d[12] * d[14] * d[35] + 2 * d[15] * d[17] * d[35];
    coeffs[293] =
        2 * d[6] * d[18] * d[27] + 2 * d[7] * d[19] * d[27] + 2 * d[8] * d[20] * d[27] + 2 * d[0] * d[24] * d[27] +
        2 * d[1] * d[25] * d[27] + 2 * d[2] * d[26] * d[27] + 2 * d[7] * d[18] * d[28] - 2 * d[6] * d[19] * d[28] -
        2 * d[1] * d[24] * d[28] + 2 * d[0] * d[25] * d[28] + 2 * d[8] * d[18] * d[29] - 2 * d[6] * d[20] * d[29] -
        2 * d[2] * d[24] * d[29] + 2 * d[0] * d[26] * d[29] + 2 * d[6] * d[21] * d[30] + 2 * d[7] * d[22] * d[30] +
        2 * d[8] * d[23] * d[30] + 2 * d[3] * d[24] * d[30] + 2 * d[4] * d[25] * d[30] + 2 * d[5] * d[26] * d[30] +
        2 * d[7] * d[21] * d[31] - 2 * d[6] * d[22] * d[31] - 2 * d[4] * d[24] * d[31] + 2 * d[3] * d[25] * d[31] +
        2 * d[8] * d[21] * d[32] - 2 * d[6] * d[23] * d[32] - 2 * d[5] * d[24] * d[32] + 2 * d[3] * d[26] * d[32] +
        2 * d[0] * d[18] * d[33] - 2 * d[1] * d[19] * d[33] - 2 * d[2] * d[20] * d[33] + 2 * d[3] * d[21] * d[33] -
        2 * d[4] * d[22] * d[33] - 2 * d[5] * d[23] * d[33] + 6 * d[6] * d[24] * d[33] + 2 * d[7] * d[25] * d[33] +
        2 * d[8] * d[26] * d[33] + 2 * d[1] * d[18] * d[34] + 2 * d[0] * d[19] * d[34] + 2 * d[4] * d[21] * d[34] +
        2 * d[3] * d[22] * d[34] + 2 * d[7] * d[24] * d[34] + 2 * d[6] * d[25] * d[34] + 2 * d[2] * d[18] * d[35] +
        2 * d[0] * d[20] * d[35] + 2 * d[5] * d[21] * d[35] + 2 * d[3] * d[23] * d[35] + 2 * d[8] * d[24] * d[35] +
        2 * d[6] * d[26] * d[35];
    coeffs[294] =
        2 * d[15] * d[18] * d[27] + 2 * d[16] * d[19] * d[27] + 2 * d[17] * d[20] * d[27] + 2 * d[9] * d[24] * d[27] +
        2 * d[10] * d[25] * d[27] + 2 * d[11] * d[26] * d[27] + 2 * d[16] * d[18] * d[28] - 2 * d[15] * d[19] * d[28] -
        2 * d[10] * d[24] * d[28] + 2 * d[9] * d[25] * d[28] + 2 * d[17] * d[18] * d[29] - 2 * d[15] * d[20] * d[29] -
        2 * d[11] * d[24] * d[29] + 2 * d[9] * d[26] * d[29] + 2 * d[15] * d[21] * d[30] + 2 * d[16] * d[22] * d[30] +
        2 * d[17] * d[23] * d[30] + 2 * d[12] * d[24] * d[30] + 2 * d[13] * d[25] * d[30] + 2 * d[14] * d[26] * d[30] +
        2 * d[16] * d[21] * d[31] - 2 * d[15] * d[22] * d[31] - 2 * d[13] * d[24] * d[31] + 2 * d[12] * d[25] * d[31] +
        2 * d[17] * d[21] * d[32] - 2 * d[15] * d[23] * d[32] - 2 * d[14] * d[24] * d[32] + 2 * d[12] * d[26] * d[32] +
        2 * d[9] * d[18] * d[33] - 2 * d[10] * d[19] * d[33] - 2 * d[11] * d[20] * d[33] + 2 * d[12] * d[21] * d[33] -
        2 * d[13] * d[22] * d[33] - 2 * d[14] * d[23] * d[33] + 6 * d[15] * d[24] * d[33] + 2 * d[16] * d[25] * d[33] +
        2 * d[17] * d[26] * d[33] + 2 * d[10] * d[18] * d[34] + 2 * d[9] * d[19] * d[34] + 2 * d[13] * d[21] * d[34] +
        2 * d[12] * d[22] * d[34] + 2 * d[16] * d[24] * d[34] + 2 * d[15] * d[25] * d[34] + 2 * d[11] * d[18] * d[35] +
        2 * d[9] * d[20] * d[35] + 2 * d[14] * d[21] * d[35] + 2 * d[12] * d[23] * d[35] + 2 * d[17] * d[24] * d[35] +
        2 * d[15] * d[26] * d[35];
    coeffs[295] = 2 * d[18] * d[24] * d[27] + 2 * d[19] * d[25] * d[27] + 2 * d[20] * d[26] * d[27] -
                  2 * d[19] * d[24] * d[28] + 2 * d[18] * d[25] * d[28] - 2 * d[20] * d[24] * d[29] +
                  2 * d[18] * d[26] * d[29] + 2 * d[21] * d[24] * d[30] + 2 * d[22] * d[25] * d[30] +
                  2 * d[23] * d[26] * d[30] - 2 * d[22] * d[24] * d[31] + 2 * d[21] * d[25] * d[31] -
                  2 * d[23] * d[24] * d[32] + 2 * d[21] * d[26] * d[32] + std::pow(d[18], 2) * d[33] -
                  std::pow(d[19], 2) * d[33] - std::pow(d[20], 2) * d[33] + std::pow(d[21], 2) * d[33] -
                  std::pow(d[22], 2) * d[33] - std::pow(d[23], 2) * d[33] + 3 * std::pow(d[24], 2) * d[33] +
                  std::pow(d[25], 2) * d[33] + std::pow(d[26], 2) * d[33] + 2 * d[18] * d[19] * d[34] +
                  2 * d[21] * d[22] * d[34] + 2 * d[24] * d[25] * d[34] + 2 * d[18] * d[20] * d[35] +
                  2 * d[21] * d[23] * d[35] + 2 * d[24] * d[26] * d[35];
    coeffs[296] =
        d[6] * std::pow(d[27], 2) + 2 * d[7] * d[27] * d[28] - d[6] * std::pow(d[28], 2) + 2 * d[8] * d[27] * d[29] -
        d[6] * std::pow(d[29], 2) + d[6] * std::pow(d[30], 2) + 2 * d[7] * d[30] * d[31] - d[6] * std::pow(d[31], 2) +
        2 * d[8] * d[30] * d[32] - d[6] * std::pow(d[32], 2) + 2 * d[0] * d[27] * d[33] - 2 * d[1] * d[28] * d[33] -
        2 * d[2] * d[29] * d[33] + 2 * d[3] * d[30] * d[33] - 2 * d[4] * d[31] * d[33] - 2 * d[5] * d[32] * d[33] +
        3 * d[6] * std::pow(d[33], 2) + 2 * d[1] * d[27] * d[34] + 2 * d[0] * d[28] * d[34] + 2 * d[4] * d[30] * d[34] +
        2 * d[3] * d[31] * d[34] + 2 * d[7] * d[33] * d[34] + d[6] * std::pow(d[34], 2) + 2 * d[2] * d[27] * d[35] +
        2 * d[0] * d[29] * d[35] + 2 * d[5] * d[30] * d[35] + 2 * d[3] * d[32] * d[35] + 2 * d[8] * d[33] * d[35] +
        d[6] * std::pow(d[35], 2);
    coeffs[297] = d[15] * std::pow(d[27], 2) + 2 * d[16] * d[27] * d[28] - d[15] * std::pow(d[28], 2) +
                  2 * d[17] * d[27] * d[29] - d[15] * std::pow(d[29], 2) + d[15] * std::pow(d[30], 2) +
                  2 * d[16] * d[30] * d[31] - d[15] * std::pow(d[31], 2) + 2 * d[17] * d[30] * d[32] -
                  d[15] * std::pow(d[32], 2) + 2 * d[9] * d[27] * d[33] - 2 * d[10] * d[28] * d[33] -
                  2 * d[11] * d[29] * d[33] + 2 * d[12] * d[30] * d[33] - 2 * d[13] * d[31] * d[33] -
                  2 * d[14] * d[32] * d[33] + 3 * d[15] * std::pow(d[33], 2) + 2 * d[10] * d[27] * d[34] +
                  2 * d[9] * d[28] * d[34] + 2 * d[13] * d[30] * d[34] + 2 * d[12] * d[31] * d[34] +
                  2 * d[16] * d[33] * d[34] + d[15] * std::pow(d[34], 2) + 2 * d[11] * d[27] * d[35] +
                  2 * d[9] * d[29] * d[35] + 2 * d[14] * d[30] * d[35] + 2 * d[12] * d[32] * d[35] +
                  2 * d[17] * d[33] * d[35] + d[15] * std::pow(d[35], 2);
    coeffs[298] = d[24] * std::pow(d[27], 2) + 2 * d[25] * d[27] * d[28] - d[24] * std::pow(d[28], 2) +
                  2 * d[26] * d[27] * d[29] - d[24] * std::pow(d[29], 2) + d[24] * std::pow(d[30], 2) +
                  2 * d[25] * d[30] * d[31] - d[24] * std::pow(d[31], 2) + 2 * d[26] * d[30] * d[32] -
                  d[24] * std::pow(d[32], 2) + 2 * d[18] * d[27] * d[33] - 2 * d[19] * d[28] * d[33] -
                  2 * d[20] * d[29] * d[33] + 2 * d[21] * d[30] * d[33] - 2 * d[22] * d[31] * d[33] -
                  2 * d[23] * d[32] * d[33] + 3 * d[24] * std::pow(d[33], 2) + 2 * d[19] * d[27] * d[34] +
                  2 * d[18] * d[28] * d[34] + 2 * d[22] * d[30] * d[34] + 2 * d[21] * d[31] * d[34] +
                  2 * d[25] * d[33] * d[34] + d[24] * std::pow(d[34], 2) + 2 * d[20] * d[27] * d[35] +
                  2 * d[18] * d[29] * d[35] + 2 * d[23] * d[30] * d[35] + 2 * d[21] * d[32] * d[35] +
                  2 * d[26] * d[33] * d[35] + d[24] * std::pow(d[35], 2);
    coeffs[299] = std::pow(d[27], 2) * d[33] - std::pow(d[28], 2) * d[33] - std::pow(d[29], 2) * d[33] +
                  std::pow(d[30], 2) * d[33] - std::pow(d[31], 2) * d[33] - std::pow(d[32], 2) * d[33] +
                  std::pow(d[33], 3) + 2 * d[27] * d[28] * d[34] + 2 * d[30] * d[31] * d[34] +
                  d[33] * std::pow(d[34], 2) + 2 * d[27] * d[29] * d[35] + 2 * d[30] * d[32] * d[35] +
                  d[33] * std::pow(d[35], 2);
    coeffs[300] =
        2 * d[0] * d[6] * d[36] + 2 * d[1] * d[7] * d[36] + 2 * d[2] * d[8] * d[36] - 2 * d[1] * d[6] * d[37] +
        2 * d[0] * d[7] * d[37] - 2 * d[2] * d[6] * d[38] + 2 * d[0] * d[8] * d[38] + 2 * d[3] * d[6] * d[39] +
        2 * d[4] * d[7] * d[39] + 2 * d[5] * d[8] * d[39] - 2 * d[4] * d[6] * d[40] + 2 * d[3] * d[7] * d[40] -
        2 * d[5] * d[6] * d[41] + 2 * d[3] * d[8] * d[41] + std::pow(d[0], 2) * d[42] - std::pow(d[1], 2) * d[42] -
        std::pow(d[2], 2) * d[42] + std::pow(d[3], 2) * d[42] - std::pow(d[4], 2) * d[42] - std::pow(d[5], 2) * d[42] +
        3 * std::pow(d[6], 2) * d[42] + std::pow(d[7], 2) * d[42] + std::pow(d[8], 2) * d[42] +
        2 * d[0] * d[1] * d[43] + 2 * d[3] * d[4] * d[43] + 2 * d[6] * d[7] * d[43] + 2 * d[0] * d[2] * d[44] +
        2 * d[3] * d[5] * d[44] + 2 * d[6] * d[8] * d[44];
    coeffs[301] =
        2 * d[6] * d[9] * d[36] + 2 * d[7] * d[10] * d[36] + 2 * d[8] * d[11] * d[36] + 2 * d[0] * d[15] * d[36] +
        2 * d[1] * d[16] * d[36] + 2 * d[2] * d[17] * d[36] + 2 * d[7] * d[9] * d[37] - 2 * d[6] * d[10] * d[37] -
        2 * d[1] * d[15] * d[37] + 2 * d[0] * d[16] * d[37] + 2 * d[8] * d[9] * d[38] - 2 * d[6] * d[11] * d[38] -
        2 * d[2] * d[15] * d[38] + 2 * d[0] * d[17] * d[38] + 2 * d[6] * d[12] * d[39] + 2 * d[7] * d[13] * d[39] +
        2 * d[8] * d[14] * d[39] + 2 * d[3] * d[15] * d[39] + 2 * d[4] * d[16] * d[39] + 2 * d[5] * d[17] * d[39] +
        2 * d[7] * d[12] * d[40] - 2 * d[6] * d[13] * d[40] - 2 * d[4] * d[15] * d[40] + 2 * d[3] * d[16] * d[40] +
        2 * d[8] * d[12] * d[41] - 2 * d[6] * d[14] * d[41] - 2 * d[5] * d[15] * d[41] + 2 * d[3] * d[17] * d[41] +
        2 * d[0] * d[9] * d[42] - 2 * d[1] * d[10] * d[42] - 2 * d[2] * d[11] * d[42] + 2 * d[3] * d[12] * d[42] -
        2 * d[4] * d[13] * d[42] - 2 * d[5] * d[14] * d[42] + 6 * d[6] * d[15] * d[42] + 2 * d[7] * d[16] * d[42] +
        2 * d[8] * d[17] * d[42] + 2 * d[1] * d[9] * d[43] + 2 * d[0] * d[10] * d[43] + 2 * d[4] * d[12] * d[43] +
        2 * d[3] * d[13] * d[43] + 2 * d[7] * d[15] * d[43] + 2 * d[6] * d[16] * d[43] + 2 * d[2] * d[9] * d[44] +
        2 * d[0] * d[11] * d[44] + 2 * d[5] * d[12] * d[44] + 2 * d[3] * d[14] * d[44] + 2 * d[8] * d[15] * d[44] +
        2 * d[6] * d[17] * d[44];
    coeffs[302] =
        2 * d[9] * d[15] * d[36] + 2 * d[10] * d[16] * d[36] + 2 * d[11] * d[17] * d[36] - 2 * d[10] * d[15] * d[37] +
        2 * d[9] * d[16] * d[37] - 2 * d[11] * d[15] * d[38] + 2 * d[9] * d[17] * d[38] + 2 * d[12] * d[15] * d[39] +
        2 * d[13] * d[16] * d[39] + 2 * d[14] * d[17] * d[39] - 2 * d[13] * d[15] * d[40] + 2 * d[12] * d[16] * d[40] -
        2 * d[14] * d[15] * d[41] + 2 * d[12] * d[17] * d[41] + std::pow(d[9], 2) * d[42] - std::pow(d[10], 2) * d[42] -
        std::pow(d[11], 2) * d[42] + std::pow(d[12], 2) * d[42] - std::pow(d[13], 2) * d[42] -
        std::pow(d[14], 2) * d[42] + 3 * std::pow(d[15], 2) * d[42] + std::pow(d[16], 2) * d[42] +
        std::pow(d[17], 2) * d[42] + 2 * d[9] * d[10] * d[43] + 2 * d[12] * d[13] * d[43] + 2 * d[15] * d[16] * d[43] +
        2 * d[9] * d[11] * d[44] + 2 * d[12] * d[14] * d[44] + 2 * d[15] * d[17] * d[44];
    coeffs[303] =
        2 * d[6] * d[18] * d[36] + 2 * d[7] * d[19] * d[36] + 2 * d[8] * d[20] * d[36] + 2 * d[0] * d[24] * d[36] +
        2 * d[1] * d[25] * d[36] + 2 * d[2] * d[26] * d[36] + 2 * d[7] * d[18] * d[37] - 2 * d[6] * d[19] * d[37] -
        2 * d[1] * d[24] * d[37] + 2 * d[0] * d[25] * d[37] + 2 * d[8] * d[18] * d[38] - 2 * d[6] * d[20] * d[38] -
        2 * d[2] * d[24] * d[38] + 2 * d[0] * d[26] * d[38] + 2 * d[6] * d[21] * d[39] + 2 * d[7] * d[22] * d[39] +
        2 * d[8] * d[23] * d[39] + 2 * d[3] * d[24] * d[39] + 2 * d[4] * d[25] * d[39] + 2 * d[5] * d[26] * d[39] +
        2 * d[7] * d[21] * d[40] - 2 * d[6] * d[22] * d[40] - 2 * d[4] * d[24] * d[40] + 2 * d[3] * d[25] * d[40] +
        2 * d[8] * d[21] * d[41] - 2 * d[6] * d[23] * d[41] - 2 * d[5] * d[24] * d[41] + 2 * d[3] * d[26] * d[41] +
        2 * d[0] * d[18] * d[42] - 2 * d[1] * d[19] * d[42] - 2 * d[2] * d[20] * d[42] + 2 * d[3] * d[21] * d[42] -
        2 * d[4] * d[22] * d[42] - 2 * d[5] * d[23] * d[42] + 6 * d[6] * d[24] * d[42] + 2 * d[7] * d[25] * d[42] +
        2 * d[8] * d[26] * d[42] + 2 * d[1] * d[18] * d[43] + 2 * d[0] * d[19] * d[43] + 2 * d[4] * d[21] * d[43] +
        2 * d[3] * d[22] * d[43] + 2 * d[7] * d[24] * d[43] + 2 * d[6] * d[25] * d[43] + 2 * d[2] * d[18] * d[44] +
        2 * d[0] * d[20] * d[44] + 2 * d[5] * d[21] * d[44] + 2 * d[3] * d[23] * d[44] + 2 * d[8] * d[24] * d[44] +
        2 * d[6] * d[26] * d[44];
    coeffs[304] =
        2 * d[15] * d[18] * d[36] + 2 * d[16] * d[19] * d[36] + 2 * d[17] * d[20] * d[36] + 2 * d[9] * d[24] * d[36] +
        2 * d[10] * d[25] * d[36] + 2 * d[11] * d[26] * d[36] + 2 * d[16] * d[18] * d[37] - 2 * d[15] * d[19] * d[37] -
        2 * d[10] * d[24] * d[37] + 2 * d[9] * d[25] * d[37] + 2 * d[17] * d[18] * d[38] - 2 * d[15] * d[20] * d[38] -
        2 * d[11] * d[24] * d[38] + 2 * d[9] * d[26] * d[38] + 2 * d[15] * d[21] * d[39] + 2 * d[16] * d[22] * d[39] +
        2 * d[17] * d[23] * d[39] + 2 * d[12] * d[24] * d[39] + 2 * d[13] * d[25] * d[39] + 2 * d[14] * d[26] * d[39] +
        2 * d[16] * d[21] * d[40] - 2 * d[15] * d[22] * d[40] - 2 * d[13] * d[24] * d[40] + 2 * d[12] * d[25] * d[40] +
        2 * d[17] * d[21] * d[41] - 2 * d[15] * d[23] * d[41] - 2 * d[14] * d[24] * d[41] + 2 * d[12] * d[26] * d[41] +
        2 * d[9] * d[18] * d[42] - 2 * d[10] * d[19] * d[42] - 2 * d[11] * d[20] * d[42] + 2 * d[12] * d[21] * d[42] -
        2 * d[13] * d[22] * d[42] - 2 * d[14] * d[23] * d[42] + 6 * d[15] * d[24] * d[42] + 2 * d[16] * d[25] * d[42] +
        2 * d[17] * d[26] * d[42] + 2 * d[10] * d[18] * d[43] + 2 * d[9] * d[19] * d[43] + 2 * d[13] * d[21] * d[43] +
        2 * d[12] * d[22] * d[43] + 2 * d[16] * d[24] * d[43] + 2 * d[15] * d[25] * d[43] + 2 * d[11] * d[18] * d[44] +
        2 * d[9] * d[20] * d[44] + 2 * d[14] * d[21] * d[44] + 2 * d[12] * d[23] * d[44] + 2 * d[17] * d[24] * d[44] +
        2 * d[15] * d[26] * d[44];
    coeffs[305] = 2 * d[18] * d[24] * d[36] + 2 * d[19] * d[25] * d[36] + 2 * d[20] * d[26] * d[36] -
                  2 * d[19] * d[24] * d[37] + 2 * d[18] * d[25] * d[37] - 2 * d[20] * d[24] * d[38] +
                  2 * d[18] * d[26] * d[38] + 2 * d[21] * d[24] * d[39] + 2 * d[22] * d[25] * d[39] +
                  2 * d[23] * d[26] * d[39] - 2 * d[22] * d[24] * d[40] + 2 * d[21] * d[25] * d[40] -
                  2 * d[23] * d[24] * d[41] + 2 * d[21] * d[26] * d[41] + std::pow(d[18], 2) * d[42] -
                  std::pow(d[19], 2) * d[42] - std::pow(d[20], 2) * d[42] + std::pow(d[21], 2) * d[42] -
                  std::pow(d[22], 2) * d[42] - std::pow(d[23], 2) * d[42] + 3 * std::pow(d[24], 2) * d[42] +
                  std::pow(d[25], 2) * d[42] + std::pow(d[26], 2) * d[42] + 2 * d[18] * d[19] * d[43] +
                  2 * d[21] * d[22] * d[43] + 2 * d[24] * d[25] * d[43] + 2 * d[18] * d[20] * d[44] +
                  2 * d[21] * d[23] * d[44] + 2 * d[24] * d[26] * d[44];
    coeffs[306] =
        2 * d[6] * d[27] * d[36] + 2 * d[7] * d[28] * d[36] + 2 * d[8] * d[29] * d[36] + 2 * d[0] * d[33] * d[36] +
        2 * d[1] * d[34] * d[36] + 2 * d[2] * d[35] * d[36] + 2 * d[7] * d[27] * d[37] - 2 * d[6] * d[28] * d[37] -
        2 * d[1] * d[33] * d[37] + 2 * d[0] * d[34] * d[37] + 2 * d[8] * d[27] * d[38] - 2 * d[6] * d[29] * d[38] -
        2 * d[2] * d[33] * d[38] + 2 * d[0] * d[35] * d[38] + 2 * d[6] * d[30] * d[39] + 2 * d[7] * d[31] * d[39] +
        2 * d[8] * d[32] * d[39] + 2 * d[3] * d[33] * d[39] + 2 * d[4] * d[34] * d[39] + 2 * d[5] * d[35] * d[39] +
        2 * d[7] * d[30] * d[40] - 2 * d[6] * d[31] * d[40] - 2 * d[4] * d[33] * d[40] + 2 * d[3] * d[34] * d[40] +
        2 * d[8] * d[30] * d[41] - 2 * d[6] * d[32] * d[41] - 2 * d[5] * d[33] * d[41] + 2 * d[3] * d[35] * d[41] +
        2 * d[0] * d[27] * d[42] - 2 * d[1] * d[28] * d[42] - 2 * d[2] * d[29] * d[42] + 2 * d[3] * d[30] * d[42] -
        2 * d[4] * d[31] * d[42] - 2 * d[5] * d[32] * d[42] + 6 * d[6] * d[33] * d[42] + 2 * d[7] * d[34] * d[42] +
        2 * d[8] * d[35] * d[42] + 2 * d[1] * d[27] * d[43] + 2 * d[0] * d[28] * d[43] + 2 * d[4] * d[30] * d[43] +
        2 * d[3] * d[31] * d[43] + 2 * d[7] * d[33] * d[43] + 2 * d[6] * d[34] * d[43] + 2 * d[2] * d[27] * d[44] +
        2 * d[0] * d[29] * d[44] + 2 * d[5] * d[30] * d[44] + 2 * d[3] * d[32] * d[44] + 2 * d[8] * d[33] * d[44] +
        2 * d[6] * d[35] * d[44];
    coeffs[307] =
        2 * d[15] * d[27] * d[36] + 2 * d[16] * d[28] * d[36] + 2 * d[17] * d[29] * d[36] + 2 * d[9] * d[33] * d[36] +
        2 * d[10] * d[34] * d[36] + 2 * d[11] * d[35] * d[36] + 2 * d[16] * d[27] * d[37] - 2 * d[15] * d[28] * d[37] -
        2 * d[10] * d[33] * d[37] + 2 * d[9] * d[34] * d[37] + 2 * d[17] * d[27] * d[38] - 2 * d[15] * d[29] * d[38] -
        2 * d[11] * d[33] * d[38] + 2 * d[9] * d[35] * d[38] + 2 * d[15] * d[30] * d[39] + 2 * d[16] * d[31] * d[39] +
        2 * d[17] * d[32] * d[39] + 2 * d[12] * d[33] * d[39] + 2 * d[13] * d[34] * d[39] + 2 * d[14] * d[35] * d[39] +
        2 * d[16] * d[30] * d[40] - 2 * d[15] * d[31] * d[40] - 2 * d[13] * d[33] * d[40] + 2 * d[12] * d[34] * d[40] +
        2 * d[17] * d[30] * d[41] - 2 * d[15] * d[32] * d[41] - 2 * d[14] * d[33] * d[41] + 2 * d[12] * d[35] * d[41] +
        2 * d[9] * d[27] * d[42] - 2 * d[10] * d[28] * d[42] - 2 * d[11] * d[29] * d[42] + 2 * d[12] * d[30] * d[42] -
        2 * d[13] * d[31] * d[42] - 2 * d[14] * d[32] * d[42] + 6 * d[15] * d[33] * d[42] + 2 * d[16] * d[34] * d[42] +
        2 * d[17] * d[35] * d[42] + 2 * d[10] * d[27] * d[43] + 2 * d[9] * d[28] * d[43] + 2 * d[13] * d[30] * d[43] +
        2 * d[12] * d[31] * d[43] + 2 * d[16] * d[33] * d[43] + 2 * d[15] * d[34] * d[43] + 2 * d[11] * d[27] * d[44] +
        2 * d[9] * d[29] * d[44] + 2 * d[14] * d[30] * d[44] + 2 * d[12] * d[32] * d[44] + 2 * d[17] * d[33] * d[44] +
        2 * d[15] * d[35] * d[44];
    coeffs[308] =
        2 * d[24] * d[27] * d[36] + 2 * d[25] * d[28] * d[36] + 2 * d[26] * d[29] * d[36] + 2 * d[18] * d[33] * d[36] +
        2 * d[19] * d[34] * d[36] + 2 * d[20] * d[35] * d[36] + 2 * d[25] * d[27] * d[37] - 2 * d[24] * d[28] * d[37] -
        2 * d[19] * d[33] * d[37] + 2 * d[18] * d[34] * d[37] + 2 * d[26] * d[27] * d[38] - 2 * d[24] * d[29] * d[38] -
        2 * d[20] * d[33] * d[38] + 2 * d[18] * d[35] * d[38] + 2 * d[24] * d[30] * d[39] + 2 * d[25] * d[31] * d[39] +
        2 * d[26] * d[32] * d[39] + 2 * d[21] * d[33] * d[39] + 2 * d[22] * d[34] * d[39] + 2 * d[23] * d[35] * d[39] +
        2 * d[25] * d[30] * d[40] - 2 * d[24] * d[31] * d[40] - 2 * d[22] * d[33] * d[40] + 2 * d[21] * d[34] * d[40] +
        2 * d[26] * d[30] * d[41] - 2 * d[24] * d[32] * d[41] - 2 * d[23] * d[33] * d[41] + 2 * d[21] * d[35] * d[41] +
        2 * d[18] * d[27] * d[42] - 2 * d[19] * d[28] * d[42] - 2 * d[20] * d[29] * d[42] + 2 * d[21] * d[30] * d[42] -
        2 * d[22] * d[31] * d[42] - 2 * d[23] * d[32] * d[42] + 6 * d[24] * d[33] * d[42] + 2 * d[25] * d[34] * d[42] +
        2 * d[26] * d[35] * d[42] + 2 * d[19] * d[27] * d[43] + 2 * d[18] * d[28] * d[43] + 2 * d[22] * d[30] * d[43] +
        2 * d[21] * d[31] * d[43] + 2 * d[25] * d[33] * d[43] + 2 * d[24] * d[34] * d[43] + 2 * d[20] * d[27] * d[44] +
        2 * d[18] * d[29] * d[44] + 2 * d[23] * d[30] * d[44] + 2 * d[21] * d[32] * d[44] + 2 * d[26] * d[33] * d[44] +
        2 * d[24] * d[35] * d[44];
    coeffs[309] = 2 * d[27] * d[33] * d[36] + 2 * d[28] * d[34] * d[36] + 2 * d[29] * d[35] * d[36] -
                  2 * d[28] * d[33] * d[37] + 2 * d[27] * d[34] * d[37] - 2 * d[29] * d[33] * d[38] +
                  2 * d[27] * d[35] * d[38] + 2 * d[30] * d[33] * d[39] + 2 * d[31] * d[34] * d[39] +
                  2 * d[32] * d[35] * d[39] - 2 * d[31] * d[33] * d[40] + 2 * d[30] * d[34] * d[40] -
                  2 * d[32] * d[33] * d[41] + 2 * d[30] * d[35] * d[41] + std::pow(d[27], 2) * d[42] -
                  std::pow(d[28], 2) * d[42] - std::pow(d[29], 2) * d[42] + std::pow(d[30], 2) * d[42] -
                  std::pow(d[31], 2) * d[42] - std::pow(d[32], 2) * d[42] + 3 * std::pow(d[33], 2) * d[42] +
                  std::pow(d[34], 2) * d[42] + std::pow(d[35], 2) * d[42] + 2 * d[27] * d[28] * d[43] +
                  2 * d[30] * d[31] * d[43] + 2 * d[33] * d[34] * d[43] + 2 * d[27] * d[29] * d[44] +
                  2 * d[30] * d[32] * d[44] + 2 * d[33] * d[35] * d[44];
    coeffs[310] =
        d[6] * std::pow(d[36], 2) + 2 * d[7] * d[36] * d[37] - d[6] * std::pow(d[37], 2) + 2 * d[8] * d[36] * d[38] -
        d[6] * std::pow(d[38], 2) + d[6] * std::pow(d[39], 2) + 2 * d[7] * d[39] * d[40] - d[6] * std::pow(d[40], 2) +
        2 * d[8] * d[39] * d[41] - d[6] * std::pow(d[41], 2) + 2 * d[0] * d[36] * d[42] - 2 * d[1] * d[37] * d[42] -
        2 * d[2] * d[38] * d[42] + 2 * d[3] * d[39] * d[42] - 2 * d[4] * d[40] * d[42] - 2 * d[5] * d[41] * d[42] +
        3 * d[6] * std::pow(d[42], 2) + 2 * d[1] * d[36] * d[43] + 2 * d[0] * d[37] * d[43] + 2 * d[4] * d[39] * d[43] +
        2 * d[3] * d[40] * d[43] + 2 * d[7] * d[42] * d[43] + d[6] * std::pow(d[43], 2) + 2 * d[2] * d[36] * d[44] +
        2 * d[0] * d[38] * d[44] + 2 * d[5] * d[39] * d[44] + 2 * d[3] * d[41] * d[44] + 2 * d[8] * d[42] * d[44] +
        d[6] * std::pow(d[44], 2);
    coeffs[311] = d[15] * std::pow(d[36], 2) + 2 * d[16] * d[36] * d[37] - d[15] * std::pow(d[37], 2) +
                  2 * d[17] * d[36] * d[38] - d[15] * std::pow(d[38], 2) + d[15] * std::pow(d[39], 2) +
                  2 * d[16] * d[39] * d[40] - d[15] * std::pow(d[40], 2) + 2 * d[17] * d[39] * d[41] -
                  d[15] * std::pow(d[41], 2) + 2 * d[9] * d[36] * d[42] - 2 * d[10] * d[37] * d[42] -
                  2 * d[11] * d[38] * d[42] + 2 * d[12] * d[39] * d[42] - 2 * d[13] * d[40] * d[42] -
                  2 * d[14] * d[41] * d[42] + 3 * d[15] * std::pow(d[42], 2) + 2 * d[10] * d[36] * d[43] +
                  2 * d[9] * d[37] * d[43] + 2 * d[13] * d[39] * d[43] + 2 * d[12] * d[40] * d[43] +
                  2 * d[16] * d[42] * d[43] + d[15] * std::pow(d[43], 2) + 2 * d[11] * d[36] * d[44] +
                  2 * d[9] * d[38] * d[44] + 2 * d[14] * d[39] * d[44] + 2 * d[12] * d[41] * d[44] +
                  2 * d[17] * d[42] * d[44] + d[15] * std::pow(d[44], 2);
    coeffs[312] = d[24] * std::pow(d[36], 2) + 2 * d[25] * d[36] * d[37] - d[24] * std::pow(d[37], 2) +
                  2 * d[26] * d[36] * d[38] - d[24] * std::pow(d[38], 2) + d[24] * std::pow(d[39], 2) +
                  2 * d[25] * d[39] * d[40] - d[24] * std::pow(d[40], 2) + 2 * d[26] * d[39] * d[41] -
                  d[24] * std::pow(d[41], 2) + 2 * d[18] * d[36] * d[42] - 2 * d[19] * d[37] * d[42] -
                  2 * d[20] * d[38] * d[42] + 2 * d[21] * d[39] * d[42] - 2 * d[22] * d[40] * d[42] -
                  2 * d[23] * d[41] * d[42] + 3 * d[24] * std::pow(d[42], 2) + 2 * d[19] * d[36] * d[43] +
                  2 * d[18] * d[37] * d[43] + 2 * d[22] * d[39] * d[43] + 2 * d[21] * d[40] * d[43] +
                  2 * d[25] * d[42] * d[43] + d[24] * std::pow(d[43], 2) + 2 * d[20] * d[36] * d[44] +
                  2 * d[18] * d[38] * d[44] + 2 * d[23] * d[39] * d[44] + 2 * d[21] * d[41] * d[44] +
                  2 * d[26] * d[42] * d[44] + d[24] * std::pow(d[44], 2);
    coeffs[313] = d[33] * std::pow(d[36], 2) + 2 * d[34] * d[36] * d[37] - d[33] * std::pow(d[37], 2) +
                  2 * d[35] * d[36] * d[38] - d[33] * std::pow(d[38], 2) + d[33] * std::pow(d[39], 2) +
                  2 * d[34] * d[39] * d[40] - d[33] * std::pow(d[40], 2) + 2 * d[35] * d[39] * d[41] -
                  d[33] * std::pow(d[41], 2) + 2 * d[27] * d[36] * d[42] - 2 * d[28] * d[37] * d[42] -
                  2 * d[29] * d[38] * d[42] + 2 * d[30] * d[39] * d[42] - 2 * d[31] * d[40] * d[42] -
                  2 * d[32] * d[41] * d[42] + 3 * d[33] * std::pow(d[42], 2) + 2 * d[28] * d[36] * d[43] +
                  2 * d[27] * d[37] * d[43] + 2 * d[31] * d[39] * d[43] + 2 * d[30] * d[40] * d[43] +
                  2 * d[34] * d[42] * d[43] + d[33] * std::pow(d[43], 2) + 2 * d[29] * d[36] * d[44] +
                  2 * d[27] * d[38] * d[44] + 2 * d[32] * d[39] * d[44] + 2 * d[30] * d[41] * d[44] +
                  2 * d[35] * d[42] * d[44] + d[33] * std::pow(d[44], 2);
    coeffs[314] = std::pow(d[36], 2) * d[42] - std::pow(d[37], 2) * d[42] - std::pow(d[38], 2) * d[42] +
                  std::pow(d[39], 2) * d[42] - std::pow(d[40], 2) * d[42] - std::pow(d[41], 2) * d[42] +
                  std::pow(d[42], 3) + 2 * d[36] * d[37] * d[43] + 2 * d[39] * d[40] * d[43] +
                  d[42] * std::pow(d[43], 2) + 2 * d[36] * d[38] * d[44] + 2 * d[39] * d[41] * d[44] +
                  d[42] * std::pow(d[44], 2);
    coeffs[315] = 2 * d[0] * d[1] * d[6] + 2 * d[3] * d[4] * d[6] - std::pow(d[0], 2) * d[7] +
                  std::pow(d[1], 2) * d[7] - std::pow(d[2], 2) * d[7] - std::pow(d[3], 2) * d[7] +
                  std::pow(d[4], 2) * d[7] - std::pow(d[5], 2) * d[7] + std::pow(d[6], 2) * d[7] + std::pow(d[7], 3) +
                  2 * d[1] * d[2] * d[8] + 2 * d[4] * d[5] * d[8] + d[7] * std::pow(d[8], 2);
    coeffs[316] =
        2 * d[1] * d[6] * d[9] - 2 * d[0] * d[7] * d[9] + 2 * d[0] * d[6] * d[10] + 2 * d[1] * d[7] * d[10] +
        2 * d[2] * d[8] * d[10] - 2 * d[2] * d[7] * d[11] + 2 * d[1] * d[8] * d[11] + 2 * d[4] * d[6] * d[12] -
        2 * d[3] * d[7] * d[12] + 2 * d[3] * d[6] * d[13] + 2 * d[4] * d[7] * d[13] + 2 * d[5] * d[8] * d[13] -
        2 * d[5] * d[7] * d[14] + 2 * d[4] * d[8] * d[14] + 2 * d[0] * d[1] * d[15] + 2 * d[3] * d[4] * d[15] +
        2 * d[6] * d[7] * d[15] - std::pow(d[0], 2) * d[16] + std::pow(d[1], 2) * d[16] - std::pow(d[2], 2) * d[16] -
        std::pow(d[3], 2) * d[16] + std::pow(d[4], 2) * d[16] - std::pow(d[5], 2) * d[16] + std::pow(d[6], 2) * d[16] +
        3 * std::pow(d[7], 2) * d[16] + std::pow(d[8], 2) * d[16] + 2 * d[1] * d[2] * d[17] + 2 * d[4] * d[5] * d[17] +
        2 * d[7] * d[8] * d[17];
    coeffs[317] =
        -d[7] * std::pow(d[9], 2) + 2 * d[6] * d[9] * d[10] + d[7] * std::pow(d[10], 2) + 2 * d[8] * d[10] * d[11] -
        d[7] * std::pow(d[11], 2) - d[7] * std::pow(d[12], 2) + 2 * d[6] * d[12] * d[13] + d[7] * std::pow(d[13], 2) +
        2 * d[8] * d[13] * d[14] - d[7] * std::pow(d[14], 2) + 2 * d[1] * d[9] * d[15] + 2 * d[0] * d[10] * d[15] +
        2 * d[4] * d[12] * d[15] + 2 * d[3] * d[13] * d[15] + d[7] * std::pow(d[15], 2) - 2 * d[0] * d[9] * d[16] +
        2 * d[1] * d[10] * d[16] - 2 * d[2] * d[11] * d[16] - 2 * d[3] * d[12] * d[16] + 2 * d[4] * d[13] * d[16] -
        2 * d[5] * d[14] * d[16] + 2 * d[6] * d[15] * d[16] + 3 * d[7] * std::pow(d[16], 2) + 2 * d[2] * d[10] * d[17] +
        2 * d[1] * d[11] * d[17] + 2 * d[5] * d[13] * d[17] + 2 * d[4] * d[14] * d[17] + 2 * d[8] * d[16] * d[17] +
        d[7] * std::pow(d[17], 2);
    coeffs[318] = 2 * d[9] * d[10] * d[15] + 2 * d[12] * d[13] * d[15] - std::pow(d[9], 2) * d[16] +
                  std::pow(d[10], 2) * d[16] - std::pow(d[11], 2) * d[16] - std::pow(d[12], 2) * d[16] +
                  std::pow(d[13], 2) * d[16] - std::pow(d[14], 2) * d[16] + std::pow(d[15], 2) * d[16] +
                  std::pow(d[16], 3) + 2 * d[10] * d[11] * d[17] + 2 * d[13] * d[14] * d[17] +
                  d[16] * std::pow(d[17], 2);
    coeffs[319] =
        2 * d[1] * d[6] * d[18] - 2 * d[0] * d[7] * d[18] + 2 * d[0] * d[6] * d[19] + 2 * d[1] * d[7] * d[19] +
        2 * d[2] * d[8] * d[19] - 2 * d[2] * d[7] * d[20] + 2 * d[1] * d[8] * d[20] + 2 * d[4] * d[6] * d[21] -
        2 * d[3] * d[7] * d[21] + 2 * d[3] * d[6] * d[22] + 2 * d[4] * d[7] * d[22] + 2 * d[5] * d[8] * d[22] -
        2 * d[5] * d[7] * d[23] + 2 * d[4] * d[8] * d[23] + 2 * d[0] * d[1] * d[24] + 2 * d[3] * d[4] * d[24] +
        2 * d[6] * d[7] * d[24] - std::pow(d[0], 2) * d[25] + std::pow(d[1], 2) * d[25] - std::pow(d[2], 2) * d[25] -
        std::pow(d[3], 2) * d[25] + std::pow(d[4], 2) * d[25] - std::pow(d[5], 2) * d[25] + std::pow(d[6], 2) * d[25] +
        3 * std::pow(d[7], 2) * d[25] + std::pow(d[8], 2) * d[25] + 2 * d[1] * d[2] * d[26] + 2 * d[4] * d[5] * d[26] +
        2 * d[7] * d[8] * d[26];
    coeffs[320] =
        -2 * d[7] * d[9] * d[18] + 2 * d[6] * d[10] * d[18] + 2 * d[1] * d[15] * d[18] - 2 * d[0] * d[16] * d[18] +
        2 * d[6] * d[9] * d[19] + 2 * d[7] * d[10] * d[19] + 2 * d[8] * d[11] * d[19] + 2 * d[0] * d[15] * d[19] +
        2 * d[1] * d[16] * d[19] + 2 * d[2] * d[17] * d[19] + 2 * d[8] * d[10] * d[20] - 2 * d[7] * d[11] * d[20] -
        2 * d[2] * d[16] * d[20] + 2 * d[1] * d[17] * d[20] - 2 * d[7] * d[12] * d[21] + 2 * d[6] * d[13] * d[21] +
        2 * d[4] * d[15] * d[21] - 2 * d[3] * d[16] * d[21] + 2 * d[6] * d[12] * d[22] + 2 * d[7] * d[13] * d[22] +
        2 * d[8] * d[14] * d[22] + 2 * d[3] * d[15] * d[22] + 2 * d[4] * d[16] * d[22] + 2 * d[5] * d[17] * d[22] +
        2 * d[8] * d[13] * d[23] - 2 * d[7] * d[14] * d[23] - 2 * d[5] * d[16] * d[23] + 2 * d[4] * d[17] * d[23] +
        2 * d[1] * d[9] * d[24] + 2 * d[0] * d[10] * d[24] + 2 * d[4] * d[12] * d[24] + 2 * d[3] * d[13] * d[24] +
        2 * d[7] * d[15] * d[24] + 2 * d[6] * d[16] * d[24] - 2 * d[0] * d[9] * d[25] + 2 * d[1] * d[10] * d[25] -
        2 * d[2] * d[11] * d[25] - 2 * d[3] * d[12] * d[25] + 2 * d[4] * d[13] * d[25] - 2 * d[5] * d[14] * d[25] +
        2 * d[6] * d[15] * d[25] + 6 * d[7] * d[16] * d[25] + 2 * d[8] * d[17] * d[25] + 2 * d[2] * d[10] * d[26] +
        2 * d[1] * d[11] * d[26] + 2 * d[5] * d[13] * d[26] + 2 * d[4] * d[14] * d[26] + 2 * d[8] * d[16] * d[26] +
        2 * d[7] * d[17] * d[26];
    coeffs[321] =
        2 * d[10] * d[15] * d[18] - 2 * d[9] * d[16] * d[18] + 2 * d[9] * d[15] * d[19] + 2 * d[10] * d[16] * d[19] +
        2 * d[11] * d[17] * d[19] - 2 * d[11] * d[16] * d[20] + 2 * d[10] * d[17] * d[20] + 2 * d[13] * d[15] * d[21] -
        2 * d[12] * d[16] * d[21] + 2 * d[12] * d[15] * d[22] + 2 * d[13] * d[16] * d[22] + 2 * d[14] * d[17] * d[22] -
        2 * d[14] * d[16] * d[23] + 2 * d[13] * d[17] * d[23] + 2 * d[9] * d[10] * d[24] + 2 * d[12] * d[13] * d[24] +
        2 * d[15] * d[16] * d[24] - std::pow(d[9], 2) * d[25] + std::pow(d[10], 2) * d[25] -
        std::pow(d[11], 2) * d[25] - std::pow(d[12], 2) * d[25] + std::pow(d[13], 2) * d[25] -
        std::pow(d[14], 2) * d[25] + std::pow(d[15], 2) * d[25] + 3 * std::pow(d[16], 2) * d[25] +
        std::pow(d[17], 2) * d[25] + 2 * d[10] * d[11] * d[26] + 2 * d[13] * d[14] * d[26] + 2 * d[16] * d[17] * d[26];
    coeffs[322] =
        -d[7] * std::pow(d[18], 2) + 2 * d[6] * d[18] * d[19] + d[7] * std::pow(d[19], 2) + 2 * d[8] * d[19] * d[20] -
        d[7] * std::pow(d[20], 2) - d[7] * std::pow(d[21], 2) + 2 * d[6] * d[21] * d[22] + d[7] * std::pow(d[22], 2) +
        2 * d[8] * d[22] * d[23] - d[7] * std::pow(d[23], 2) + 2 * d[1] * d[18] * d[24] + 2 * d[0] * d[19] * d[24] +
        2 * d[4] * d[21] * d[24] + 2 * d[3] * d[22] * d[24] + d[7] * std::pow(d[24], 2) - 2 * d[0] * d[18] * d[25] +
        2 * d[1] * d[19] * d[25] - 2 * d[2] * d[20] * d[25] - 2 * d[3] * d[21] * d[25] + 2 * d[4] * d[22] * d[25] -
        2 * d[5] * d[23] * d[25] + 2 * d[6] * d[24] * d[25] + 3 * d[7] * std::pow(d[25], 2) + 2 * d[2] * d[19] * d[26] +
        2 * d[1] * d[20] * d[26] + 2 * d[5] * d[22] * d[26] + 2 * d[4] * d[23] * d[26] + 2 * d[8] * d[25] * d[26] +
        d[7] * std::pow(d[26], 2);
    coeffs[323] = -d[16] * std::pow(d[18], 2) + 2 * d[15] * d[18] * d[19] + d[16] * std::pow(d[19], 2) +
                  2 * d[17] * d[19] * d[20] - d[16] * std::pow(d[20], 2) - d[16] * std::pow(d[21], 2) +
                  2 * d[15] * d[21] * d[22] + d[16] * std::pow(d[22], 2) + 2 * d[17] * d[22] * d[23] -
                  d[16] * std::pow(d[23], 2) + 2 * d[10] * d[18] * d[24] + 2 * d[9] * d[19] * d[24] +
                  2 * d[13] * d[21] * d[24] + 2 * d[12] * d[22] * d[24] + d[16] * std::pow(d[24], 2) -
                  2 * d[9] * d[18] * d[25] + 2 * d[10] * d[19] * d[25] - 2 * d[11] * d[20] * d[25] -
                  2 * d[12] * d[21] * d[25] + 2 * d[13] * d[22] * d[25] - 2 * d[14] * d[23] * d[25] +
                  2 * d[15] * d[24] * d[25] + 3 * d[16] * std::pow(d[25], 2) + 2 * d[11] * d[19] * d[26] +
                  2 * d[10] * d[20] * d[26] + 2 * d[14] * d[22] * d[26] + 2 * d[13] * d[23] * d[26] +
                  2 * d[17] * d[25] * d[26] + d[16] * std::pow(d[26], 2);
    coeffs[324] = 2 * d[18] * d[19] * d[24] + 2 * d[21] * d[22] * d[24] - std::pow(d[18], 2) * d[25] +
                  std::pow(d[19], 2) * d[25] - std::pow(d[20], 2) * d[25] - std::pow(d[21], 2) * d[25] +
                  std::pow(d[22], 2) * d[25] - std::pow(d[23], 2) * d[25] + std::pow(d[24], 2) * d[25] +
                  std::pow(d[25], 3) + 2 * d[19] * d[20] * d[26] + 2 * d[22] * d[23] * d[26] +
                  d[25] * std::pow(d[26], 2);
    coeffs[325] =
        2 * d[1] * d[6] * d[27] - 2 * d[0] * d[7] * d[27] + 2 * d[0] * d[6] * d[28] + 2 * d[1] * d[7] * d[28] +
        2 * d[2] * d[8] * d[28] - 2 * d[2] * d[7] * d[29] + 2 * d[1] * d[8] * d[29] + 2 * d[4] * d[6] * d[30] -
        2 * d[3] * d[7] * d[30] + 2 * d[3] * d[6] * d[31] + 2 * d[4] * d[7] * d[31] + 2 * d[5] * d[8] * d[31] -
        2 * d[5] * d[7] * d[32] + 2 * d[4] * d[8] * d[32] + 2 * d[0] * d[1] * d[33] + 2 * d[3] * d[4] * d[33] +
        2 * d[6] * d[7] * d[33] - std::pow(d[0], 2) * d[34] + std::pow(d[1], 2) * d[34] - std::pow(d[2], 2) * d[34] -
        std::pow(d[3], 2) * d[34] + std::pow(d[4], 2) * d[34] - std::pow(d[5], 2) * d[34] + std::pow(d[6], 2) * d[34] +
        3 * std::pow(d[7], 2) * d[34] + std::pow(d[8], 2) * d[34] + 2 * d[1] * d[2] * d[35] + 2 * d[4] * d[5] * d[35] +
        2 * d[7] * d[8] * d[35];
    coeffs[326] =
        -2 * d[7] * d[9] * d[27] + 2 * d[6] * d[10] * d[27] + 2 * d[1] * d[15] * d[27] - 2 * d[0] * d[16] * d[27] +
        2 * d[6] * d[9] * d[28] + 2 * d[7] * d[10] * d[28] + 2 * d[8] * d[11] * d[28] + 2 * d[0] * d[15] * d[28] +
        2 * d[1] * d[16] * d[28] + 2 * d[2] * d[17] * d[28] + 2 * d[8] * d[10] * d[29] - 2 * d[7] * d[11] * d[29] -
        2 * d[2] * d[16] * d[29] + 2 * d[1] * d[17] * d[29] - 2 * d[7] * d[12] * d[30] + 2 * d[6] * d[13] * d[30] +
        2 * d[4] * d[15] * d[30] - 2 * d[3] * d[16] * d[30] + 2 * d[6] * d[12] * d[31] + 2 * d[7] * d[13] * d[31] +
        2 * d[8] * d[14] * d[31] + 2 * d[3] * d[15] * d[31] + 2 * d[4] * d[16] * d[31] + 2 * d[5] * d[17] * d[31] +
        2 * d[8] * d[13] * d[32] - 2 * d[7] * d[14] * d[32] - 2 * d[5] * d[16] * d[32] + 2 * d[4] * d[17] * d[32] +
        2 * d[1] * d[9] * d[33] + 2 * d[0] * d[10] * d[33] + 2 * d[4] * d[12] * d[33] + 2 * d[3] * d[13] * d[33] +
        2 * d[7] * d[15] * d[33] + 2 * d[6] * d[16] * d[33] - 2 * d[0] * d[9] * d[34] + 2 * d[1] * d[10] * d[34] -
        2 * d[2] * d[11] * d[34] - 2 * d[3] * d[12] * d[34] + 2 * d[4] * d[13] * d[34] - 2 * d[5] * d[14] * d[34] +
        2 * d[6] * d[15] * d[34] + 6 * d[7] * d[16] * d[34] + 2 * d[8] * d[17] * d[34] + 2 * d[2] * d[10] * d[35] +
        2 * d[1] * d[11] * d[35] + 2 * d[5] * d[13] * d[35] + 2 * d[4] * d[14] * d[35] + 2 * d[8] * d[16] * d[35] +
        2 * d[7] * d[17] * d[35];
    coeffs[327] =
        2 * d[10] * d[15] * d[27] - 2 * d[9] * d[16] * d[27] + 2 * d[9] * d[15] * d[28] + 2 * d[10] * d[16] * d[28] +
        2 * d[11] * d[17] * d[28] - 2 * d[11] * d[16] * d[29] + 2 * d[10] * d[17] * d[29] + 2 * d[13] * d[15] * d[30] -
        2 * d[12] * d[16] * d[30] + 2 * d[12] * d[15] * d[31] + 2 * d[13] * d[16] * d[31] + 2 * d[14] * d[17] * d[31] -
        2 * d[14] * d[16] * d[32] + 2 * d[13] * d[17] * d[32] + 2 * d[9] * d[10] * d[33] + 2 * d[12] * d[13] * d[33] +
        2 * d[15] * d[16] * d[33] - std::pow(d[9], 2) * d[34] + std::pow(d[10], 2) * d[34] -
        std::pow(d[11], 2) * d[34] - std::pow(d[12], 2) * d[34] + std::pow(d[13], 2) * d[34] -
        std::pow(d[14], 2) * d[34] + std::pow(d[15], 2) * d[34] + 3 * std::pow(d[16], 2) * d[34] +
        std::pow(d[17], 2) * d[34] + 2 * d[10] * d[11] * d[35] + 2 * d[13] * d[14] * d[35] + 2 * d[16] * d[17] * d[35];
    coeffs[328] =
        -2 * d[7] * d[18] * d[27] + 2 * d[6] * d[19] * d[27] + 2 * d[1] * d[24] * d[27] - 2 * d[0] * d[25] * d[27] +
        2 * d[6] * d[18] * d[28] + 2 * d[7] * d[19] * d[28] + 2 * d[8] * d[20] * d[28] + 2 * d[0] * d[24] * d[28] +
        2 * d[1] * d[25] * d[28] + 2 * d[2] * d[26] * d[28] + 2 * d[8] * d[19] * d[29] - 2 * d[7] * d[20] * d[29] -
        2 * d[2] * d[25] * d[29] + 2 * d[1] * d[26] * d[29] - 2 * d[7] * d[21] * d[30] + 2 * d[6] * d[22] * d[30] +
        2 * d[4] * d[24] * d[30] - 2 * d[3] * d[25] * d[30] + 2 * d[6] * d[21] * d[31] + 2 * d[7] * d[22] * d[31] +
        2 * d[8] * d[23] * d[31] + 2 * d[3] * d[24] * d[31] + 2 * d[4] * d[25] * d[31] + 2 * d[5] * d[26] * d[31] +
        2 * d[8] * d[22] * d[32] - 2 * d[7] * d[23] * d[32] - 2 * d[5] * d[25] * d[32] + 2 * d[4] * d[26] * d[32] +
        2 * d[1] * d[18] * d[33] + 2 * d[0] * d[19] * d[33] + 2 * d[4] * d[21] * d[33] + 2 * d[3] * d[22] * d[33] +
        2 * d[7] * d[24] * d[33] + 2 * d[6] * d[25] * d[33] - 2 * d[0] * d[18] * d[34] + 2 * d[1] * d[19] * d[34] -
        2 * d[2] * d[20] * d[34] - 2 * d[3] * d[21] * d[34] + 2 * d[4] * d[22] * d[34] - 2 * d[5] * d[23] * d[34] +
        2 * d[6] * d[24] * d[34] + 6 * d[7] * d[25] * d[34] + 2 * d[8] * d[26] * d[34] + 2 * d[2] * d[19] * d[35] +
        2 * d[1] * d[20] * d[35] + 2 * d[5] * d[22] * d[35] + 2 * d[4] * d[23] * d[35] + 2 * d[8] * d[25] * d[35] +
        2 * d[7] * d[26] * d[35];
    coeffs[329] =
        -2 * d[16] * d[18] * d[27] + 2 * d[15] * d[19] * d[27] + 2 * d[10] * d[24] * d[27] - 2 * d[9] * d[25] * d[27] +
        2 * d[15] * d[18] * d[28] + 2 * d[16] * d[19] * d[28] + 2 * d[17] * d[20] * d[28] + 2 * d[9] * d[24] * d[28] +
        2 * d[10] * d[25] * d[28] + 2 * d[11] * d[26] * d[28] + 2 * d[17] * d[19] * d[29] - 2 * d[16] * d[20] * d[29] -
        2 * d[11] * d[25] * d[29] + 2 * d[10] * d[26] * d[29] - 2 * d[16] * d[21] * d[30] + 2 * d[15] * d[22] * d[30] +
        2 * d[13] * d[24] * d[30] - 2 * d[12] * d[25] * d[30] + 2 * d[15] * d[21] * d[31] + 2 * d[16] * d[22] * d[31] +
        2 * d[17] * d[23] * d[31] + 2 * d[12] * d[24] * d[31] + 2 * d[13] * d[25] * d[31] + 2 * d[14] * d[26] * d[31] +
        2 * d[17] * d[22] * d[32] - 2 * d[16] * d[23] * d[32] - 2 * d[14] * d[25] * d[32] + 2 * d[13] * d[26] * d[32] +
        2 * d[10] * d[18] * d[33] + 2 * d[9] * d[19] * d[33] + 2 * d[13] * d[21] * d[33] + 2 * d[12] * d[22] * d[33] +
        2 * d[16] * d[24] * d[33] + 2 * d[15] * d[25] * d[33] - 2 * d[9] * d[18] * d[34] + 2 * d[10] * d[19] * d[34] -
        2 * d[11] * d[20] * d[34] - 2 * d[12] * d[21] * d[34] + 2 * d[13] * d[22] * d[34] - 2 * d[14] * d[23] * d[34] +
        2 * d[15] * d[24] * d[34] + 6 * d[16] * d[25] * d[34] + 2 * d[17] * d[26] * d[34] + 2 * d[11] * d[19] * d[35] +
        2 * d[10] * d[20] * d[35] + 2 * d[14] * d[22] * d[35] + 2 * d[13] * d[23] * d[35] + 2 * d[17] * d[25] * d[35] +
        2 * d[16] * d[26] * d[35];
    coeffs[330] =
        2 * d[19] * d[24] * d[27] - 2 * d[18] * d[25] * d[27] + 2 * d[18] * d[24] * d[28] + 2 * d[19] * d[25] * d[28] +
        2 * d[20] * d[26] * d[28] - 2 * d[20] * d[25] * d[29] + 2 * d[19] * d[26] * d[29] + 2 * d[22] * d[24] * d[30] -
        2 * d[21] * d[25] * d[30] + 2 * d[21] * d[24] * d[31] + 2 * d[22] * d[25] * d[31] + 2 * d[23] * d[26] * d[31] -
        2 * d[23] * d[25] * d[32] + 2 * d[22] * d[26] * d[32] + 2 * d[18] * d[19] * d[33] + 2 * d[21] * d[22] * d[33] +
        2 * d[24] * d[25] * d[33] - std::pow(d[18], 2) * d[34] + std::pow(d[19], 2) * d[34] -
        std::pow(d[20], 2) * d[34] - std::pow(d[21], 2) * d[34] + std::pow(d[22], 2) * d[34] -
        std::pow(d[23], 2) * d[34] + std::pow(d[24], 2) * d[34] + 3 * std::pow(d[25], 2) * d[34] +
        std::pow(d[26], 2) * d[34] + 2 * d[19] * d[20] * d[35] + 2 * d[22] * d[23] * d[35] + 2 * d[25] * d[26] * d[35];
    coeffs[331] =
        -d[7] * std::pow(d[27], 2) + 2 * d[6] * d[27] * d[28] + d[7] * std::pow(d[28], 2) + 2 * d[8] * d[28] * d[29] -
        d[7] * std::pow(d[29], 2) - d[7] * std::pow(d[30], 2) + 2 * d[6] * d[30] * d[31] + d[7] * std::pow(d[31], 2) +
        2 * d[8] * d[31] * d[32] - d[7] * std::pow(d[32], 2) + 2 * d[1] * d[27] * d[33] + 2 * d[0] * d[28] * d[33] +
        2 * d[4] * d[30] * d[33] + 2 * d[3] * d[31] * d[33] + d[7] * std::pow(d[33], 2) - 2 * d[0] * d[27] * d[34] +
        2 * d[1] * d[28] * d[34] - 2 * d[2] * d[29] * d[34] - 2 * d[3] * d[30] * d[34] + 2 * d[4] * d[31] * d[34] -
        2 * d[5] * d[32] * d[34] + 2 * d[6] * d[33] * d[34] + 3 * d[7] * std::pow(d[34], 2) + 2 * d[2] * d[28] * d[35] +
        2 * d[1] * d[29] * d[35] + 2 * d[5] * d[31] * d[35] + 2 * d[4] * d[32] * d[35] + 2 * d[8] * d[34] * d[35] +
        d[7] * std::pow(d[35], 2);
    coeffs[332] = -d[16] * std::pow(d[27], 2) + 2 * d[15] * d[27] * d[28] + d[16] * std::pow(d[28], 2) +
                  2 * d[17] * d[28] * d[29] - d[16] * std::pow(d[29], 2) - d[16] * std::pow(d[30], 2) +
                  2 * d[15] * d[30] * d[31] + d[16] * std::pow(d[31], 2) + 2 * d[17] * d[31] * d[32] -
                  d[16] * std::pow(d[32], 2) + 2 * d[10] * d[27] * d[33] + 2 * d[9] * d[28] * d[33] +
                  2 * d[13] * d[30] * d[33] + 2 * d[12] * d[31] * d[33] + d[16] * std::pow(d[33], 2) -
                  2 * d[9] * d[27] * d[34] + 2 * d[10] * d[28] * d[34] - 2 * d[11] * d[29] * d[34] -
                  2 * d[12] * d[30] * d[34] + 2 * d[13] * d[31] * d[34] - 2 * d[14] * d[32] * d[34] +
                  2 * d[15] * d[33] * d[34] + 3 * d[16] * std::pow(d[34], 2) + 2 * d[11] * d[28] * d[35] +
                  2 * d[10] * d[29] * d[35] + 2 * d[14] * d[31] * d[35] + 2 * d[13] * d[32] * d[35] +
                  2 * d[17] * d[34] * d[35] + d[16] * std::pow(d[35], 2);
    coeffs[333] = -d[25] * std::pow(d[27], 2) + 2 * d[24] * d[27] * d[28] + d[25] * std::pow(d[28], 2) +
                  2 * d[26] * d[28] * d[29] - d[25] * std::pow(d[29], 2) - d[25] * std::pow(d[30], 2) +
                  2 * d[24] * d[30] * d[31] + d[25] * std::pow(d[31], 2) + 2 * d[26] * d[31] * d[32] -
                  d[25] * std::pow(d[32], 2) + 2 * d[19] * d[27] * d[33] + 2 * d[18] * d[28] * d[33] +
                  2 * d[22] * d[30] * d[33] + 2 * d[21] * d[31] * d[33] + d[25] * std::pow(d[33], 2) -
                  2 * d[18] * d[27] * d[34] + 2 * d[19] * d[28] * d[34] - 2 * d[20] * d[29] * d[34] -
                  2 * d[21] * d[30] * d[34] + 2 * d[22] * d[31] * d[34] - 2 * d[23] * d[32] * d[34] +
                  2 * d[24] * d[33] * d[34] + 3 * d[25] * std::pow(d[34], 2) + 2 * d[20] * d[28] * d[35] +
                  2 * d[19] * d[29] * d[35] + 2 * d[23] * d[31] * d[35] + 2 * d[22] * d[32] * d[35] +
                  2 * d[26] * d[34] * d[35] + d[25] * std::pow(d[35], 2);
    coeffs[334] = 2 * d[27] * d[28] * d[33] + 2 * d[30] * d[31] * d[33] - std::pow(d[27], 2) * d[34] +
                  std::pow(d[28], 2) * d[34] - std::pow(d[29], 2) * d[34] - std::pow(d[30], 2) * d[34] +
                  std::pow(d[31], 2) * d[34] - std::pow(d[32], 2) * d[34] + std::pow(d[33], 2) * d[34] +
                  std::pow(d[34], 3) + 2 * d[28] * d[29] * d[35] + 2 * d[31] * d[32] * d[35] +
                  d[34] * std::pow(d[35], 2);
    coeffs[335] =
        2 * d[1] * d[6] * d[36] - 2 * d[0] * d[7] * d[36] + 2 * d[0] * d[6] * d[37] + 2 * d[1] * d[7] * d[37] +
        2 * d[2] * d[8] * d[37] - 2 * d[2] * d[7] * d[38] + 2 * d[1] * d[8] * d[38] + 2 * d[4] * d[6] * d[39] -
        2 * d[3] * d[7] * d[39] + 2 * d[3] * d[6] * d[40] + 2 * d[4] * d[7] * d[40] + 2 * d[5] * d[8] * d[40] -
        2 * d[5] * d[7] * d[41] + 2 * d[4] * d[8] * d[41] + 2 * d[0] * d[1] * d[42] + 2 * d[3] * d[4] * d[42] +
        2 * d[6] * d[7] * d[42] - std::pow(d[0], 2) * d[43] + std::pow(d[1], 2) * d[43] - std::pow(d[2], 2) * d[43] -
        std::pow(d[3], 2) * d[43] + std::pow(d[4], 2) * d[43] - std::pow(d[5], 2) * d[43] + std::pow(d[6], 2) * d[43] +
        3 * std::pow(d[7], 2) * d[43] + std::pow(d[8], 2) * d[43] + 2 * d[1] * d[2] * d[44] + 2 * d[4] * d[5] * d[44] +
        2 * d[7] * d[8] * d[44];
    coeffs[336] =
        -2 * d[7] * d[9] * d[36] + 2 * d[6] * d[10] * d[36] + 2 * d[1] * d[15] * d[36] - 2 * d[0] * d[16] * d[36] +
        2 * d[6] * d[9] * d[37] + 2 * d[7] * d[10] * d[37] + 2 * d[8] * d[11] * d[37] + 2 * d[0] * d[15] * d[37] +
        2 * d[1] * d[16] * d[37] + 2 * d[2] * d[17] * d[37] + 2 * d[8] * d[10] * d[38] - 2 * d[7] * d[11] * d[38] -
        2 * d[2] * d[16] * d[38] + 2 * d[1] * d[17] * d[38] - 2 * d[7] * d[12] * d[39] + 2 * d[6] * d[13] * d[39] +
        2 * d[4] * d[15] * d[39] - 2 * d[3] * d[16] * d[39] + 2 * d[6] * d[12] * d[40] + 2 * d[7] * d[13] * d[40] +
        2 * d[8] * d[14] * d[40] + 2 * d[3] * d[15] * d[40] + 2 * d[4] * d[16] * d[40] + 2 * d[5] * d[17] * d[40] +
        2 * d[8] * d[13] * d[41] - 2 * d[7] * d[14] * d[41] - 2 * d[5] * d[16] * d[41] + 2 * d[4] * d[17] * d[41] +
        2 * d[1] * d[9] * d[42] + 2 * d[0] * d[10] * d[42] + 2 * d[4] * d[12] * d[42] + 2 * d[3] * d[13] * d[42] +
        2 * d[7] * d[15] * d[42] + 2 * d[6] * d[16] * d[42] - 2 * d[0] * d[9] * d[43] + 2 * d[1] * d[10] * d[43] -
        2 * d[2] * d[11] * d[43] - 2 * d[3] * d[12] * d[43] + 2 * d[4] * d[13] * d[43] - 2 * d[5] * d[14] * d[43] +
        2 * d[6] * d[15] * d[43] + 6 * d[7] * d[16] * d[43] + 2 * d[8] * d[17] * d[43] + 2 * d[2] * d[10] * d[44] +
        2 * d[1] * d[11] * d[44] + 2 * d[5] * d[13] * d[44] + 2 * d[4] * d[14] * d[44] + 2 * d[8] * d[16] * d[44] +
        2 * d[7] * d[17] * d[44];
    coeffs[337] =
        2 * d[10] * d[15] * d[36] - 2 * d[9] * d[16] * d[36] + 2 * d[9] * d[15] * d[37] + 2 * d[10] * d[16] * d[37] +
        2 * d[11] * d[17] * d[37] - 2 * d[11] * d[16] * d[38] + 2 * d[10] * d[17] * d[38] + 2 * d[13] * d[15] * d[39] -
        2 * d[12] * d[16] * d[39] + 2 * d[12] * d[15] * d[40] + 2 * d[13] * d[16] * d[40] + 2 * d[14] * d[17] * d[40] -
        2 * d[14] * d[16] * d[41] + 2 * d[13] * d[17] * d[41] + 2 * d[9] * d[10] * d[42] + 2 * d[12] * d[13] * d[42] +
        2 * d[15] * d[16] * d[42] - std::pow(d[9], 2) * d[43] + std::pow(d[10], 2) * d[43] -
        std::pow(d[11], 2) * d[43] - std::pow(d[12], 2) * d[43] + std::pow(d[13], 2) * d[43] -
        std::pow(d[14], 2) * d[43] + std::pow(d[15], 2) * d[43] + 3 * std::pow(d[16], 2) * d[43] +
        std::pow(d[17], 2) * d[43] + 2 * d[10] * d[11] * d[44] + 2 * d[13] * d[14] * d[44] + 2 * d[16] * d[17] * d[44];
    coeffs[338] =
        -2 * d[7] * d[18] * d[36] + 2 * d[6] * d[19] * d[36] + 2 * d[1] * d[24] * d[36] - 2 * d[0] * d[25] * d[36] +
        2 * d[6] * d[18] * d[37] + 2 * d[7] * d[19] * d[37] + 2 * d[8] * d[20] * d[37] + 2 * d[0] * d[24] * d[37] +
        2 * d[1] * d[25] * d[37] + 2 * d[2] * d[26] * d[37] + 2 * d[8] * d[19] * d[38] - 2 * d[7] * d[20] * d[38] -
        2 * d[2] * d[25] * d[38] + 2 * d[1] * d[26] * d[38] - 2 * d[7] * d[21] * d[39] + 2 * d[6] * d[22] * d[39] +
        2 * d[4] * d[24] * d[39] - 2 * d[3] * d[25] * d[39] + 2 * d[6] * d[21] * d[40] + 2 * d[7] * d[22] * d[40] +
        2 * d[8] * d[23] * d[40] + 2 * d[3] * d[24] * d[40] + 2 * d[4] * d[25] * d[40] + 2 * d[5] * d[26] * d[40] +
        2 * d[8] * d[22] * d[41] - 2 * d[7] * d[23] * d[41] - 2 * d[5] * d[25] * d[41] + 2 * d[4] * d[26] * d[41] +
        2 * d[1] * d[18] * d[42] + 2 * d[0] * d[19] * d[42] + 2 * d[4] * d[21] * d[42] + 2 * d[3] * d[22] * d[42] +
        2 * d[7] * d[24] * d[42] + 2 * d[6] * d[25] * d[42] - 2 * d[0] * d[18] * d[43] + 2 * d[1] * d[19] * d[43] -
        2 * d[2] * d[20] * d[43] - 2 * d[3] * d[21] * d[43] + 2 * d[4] * d[22] * d[43] - 2 * d[5] * d[23] * d[43] +
        2 * d[6] * d[24] * d[43] + 6 * d[7] * d[25] * d[43] + 2 * d[8] * d[26] * d[43] + 2 * d[2] * d[19] * d[44] +
        2 * d[1] * d[20] * d[44] + 2 * d[5] * d[22] * d[44] + 2 * d[4] * d[23] * d[44] + 2 * d[8] * d[25] * d[44] +
        2 * d[7] * d[26] * d[44];
    coeffs[339] =
        -2 * d[16] * d[18] * d[36] + 2 * d[15] * d[19] * d[36] + 2 * d[10] * d[24] * d[36] - 2 * d[9] * d[25] * d[36] +
        2 * d[15] * d[18] * d[37] + 2 * d[16] * d[19] * d[37] + 2 * d[17] * d[20] * d[37] + 2 * d[9] * d[24] * d[37] +
        2 * d[10] * d[25] * d[37] + 2 * d[11] * d[26] * d[37] + 2 * d[17] * d[19] * d[38] - 2 * d[16] * d[20] * d[38] -
        2 * d[11] * d[25] * d[38] + 2 * d[10] * d[26] * d[38] - 2 * d[16] * d[21] * d[39] + 2 * d[15] * d[22] * d[39] +
        2 * d[13] * d[24] * d[39] - 2 * d[12] * d[25] * d[39] + 2 * d[15] * d[21] * d[40] + 2 * d[16] * d[22] * d[40] +
        2 * d[17] * d[23] * d[40] + 2 * d[12] * d[24] * d[40] + 2 * d[13] * d[25] * d[40] + 2 * d[14] * d[26] * d[40] +
        2 * d[17] * d[22] * d[41] - 2 * d[16] * d[23] * d[41] - 2 * d[14] * d[25] * d[41] + 2 * d[13] * d[26] * d[41] +
        2 * d[10] * d[18] * d[42] + 2 * d[9] * d[19] * d[42] + 2 * d[13] * d[21] * d[42] + 2 * d[12] * d[22] * d[42] +
        2 * d[16] * d[24] * d[42] + 2 * d[15] * d[25] * d[42] - 2 * d[9] * d[18] * d[43] + 2 * d[10] * d[19] * d[43] -
        2 * d[11] * d[20] * d[43] - 2 * d[12] * d[21] * d[43] + 2 * d[13] * d[22] * d[43] - 2 * d[14] * d[23] * d[43] +
        2 * d[15] * d[24] * d[43] + 6 * d[16] * d[25] * d[43] + 2 * d[17] * d[26] * d[43] + 2 * d[11] * d[19] * d[44] +
        2 * d[10] * d[20] * d[44] + 2 * d[14] * d[22] * d[44] + 2 * d[13] * d[23] * d[44] + 2 * d[17] * d[25] * d[44] +
        2 * d[16] * d[26] * d[44];
    coeffs[340] =
        2 * d[19] * d[24] * d[36] - 2 * d[18] * d[25] * d[36] + 2 * d[18] * d[24] * d[37] + 2 * d[19] * d[25] * d[37] +
        2 * d[20] * d[26] * d[37] - 2 * d[20] * d[25] * d[38] + 2 * d[19] * d[26] * d[38] + 2 * d[22] * d[24] * d[39] -
        2 * d[21] * d[25] * d[39] + 2 * d[21] * d[24] * d[40] + 2 * d[22] * d[25] * d[40] + 2 * d[23] * d[26] * d[40] -
        2 * d[23] * d[25] * d[41] + 2 * d[22] * d[26] * d[41] + 2 * d[18] * d[19] * d[42] + 2 * d[21] * d[22] * d[42] +
        2 * d[24] * d[25] * d[42] - std::pow(d[18], 2) * d[43] + std::pow(d[19], 2) * d[43] -
        std::pow(d[20], 2) * d[43] - std::pow(d[21], 2) * d[43] + std::pow(d[22], 2) * d[43] -
        std::pow(d[23], 2) * d[43] + std::pow(d[24], 2) * d[43] + 3 * std::pow(d[25], 2) * d[43] +
        std::pow(d[26], 2) * d[43] + 2 * d[19] * d[20] * d[44] + 2 * d[22] * d[23] * d[44] + 2 * d[25] * d[26] * d[44];
    coeffs[341] =
        -2 * d[7] * d[27] * d[36] + 2 * d[6] * d[28] * d[36] + 2 * d[1] * d[33] * d[36] - 2 * d[0] * d[34] * d[36] +
        2 * d[6] * d[27] * d[37] + 2 * d[7] * d[28] * d[37] + 2 * d[8] * d[29] * d[37] + 2 * d[0] * d[33] * d[37] +
        2 * d[1] * d[34] * d[37] + 2 * d[2] * d[35] * d[37] + 2 * d[8] * d[28] * d[38] - 2 * d[7] * d[29] * d[38] -
        2 * d[2] * d[34] * d[38] + 2 * d[1] * d[35] * d[38] - 2 * d[7] * d[30] * d[39] + 2 * d[6] * d[31] * d[39] +
        2 * d[4] * d[33] * d[39] - 2 * d[3] * d[34] * d[39] + 2 * d[6] * d[30] * d[40] + 2 * d[7] * d[31] * d[40] +
        2 * d[8] * d[32] * d[40] + 2 * d[3] * d[33] * d[40] + 2 * d[4] * d[34] * d[40] + 2 * d[5] * d[35] * d[40] +
        2 * d[8] * d[31] * d[41] - 2 * d[7] * d[32] * d[41] - 2 * d[5] * d[34] * d[41] + 2 * d[4] * d[35] * d[41] +
        2 * d[1] * d[27] * d[42] + 2 * d[0] * d[28] * d[42] + 2 * d[4] * d[30] * d[42] + 2 * d[3] * d[31] * d[42] +
        2 * d[7] * d[33] * d[42] + 2 * d[6] * d[34] * d[42] - 2 * d[0] * d[27] * d[43] + 2 * d[1] * d[28] * d[43] -
        2 * d[2] * d[29] * d[43] - 2 * d[3] * d[30] * d[43] + 2 * d[4] * d[31] * d[43] - 2 * d[5] * d[32] * d[43] +
        2 * d[6] * d[33] * d[43] + 6 * d[7] * d[34] * d[43] + 2 * d[8] * d[35] * d[43] + 2 * d[2] * d[28] * d[44] +
        2 * d[1] * d[29] * d[44] + 2 * d[5] * d[31] * d[44] + 2 * d[4] * d[32] * d[44] + 2 * d[8] * d[34] * d[44] +
        2 * d[7] * d[35] * d[44];
    coeffs[342] =
        -2 * d[16] * d[27] * d[36] + 2 * d[15] * d[28] * d[36] + 2 * d[10] * d[33] * d[36] - 2 * d[9] * d[34] * d[36] +
        2 * d[15] * d[27] * d[37] + 2 * d[16] * d[28] * d[37] + 2 * d[17] * d[29] * d[37] + 2 * d[9] * d[33] * d[37] +
        2 * d[10] * d[34] * d[37] + 2 * d[11] * d[35] * d[37] + 2 * d[17] * d[28] * d[38] - 2 * d[16] * d[29] * d[38] -
        2 * d[11] * d[34] * d[38] + 2 * d[10] * d[35] * d[38] - 2 * d[16] * d[30] * d[39] + 2 * d[15] * d[31] * d[39] +
        2 * d[13] * d[33] * d[39] - 2 * d[12] * d[34] * d[39] + 2 * d[15] * d[30] * d[40] + 2 * d[16] * d[31] * d[40] +
        2 * d[17] * d[32] * d[40] + 2 * d[12] * d[33] * d[40] + 2 * d[13] * d[34] * d[40] + 2 * d[14] * d[35] * d[40] +
        2 * d[17] * d[31] * d[41] - 2 * d[16] * d[32] * d[41] - 2 * d[14] * d[34] * d[41] + 2 * d[13] * d[35] * d[41] +
        2 * d[10] * d[27] * d[42] + 2 * d[9] * d[28] * d[42] + 2 * d[13] * d[30] * d[42] + 2 * d[12] * d[31] * d[42] +
        2 * d[16] * d[33] * d[42] + 2 * d[15] * d[34] * d[42] - 2 * d[9] * d[27] * d[43] + 2 * d[10] * d[28] * d[43] -
        2 * d[11] * d[29] * d[43] - 2 * d[12] * d[30] * d[43] + 2 * d[13] * d[31] * d[43] - 2 * d[14] * d[32] * d[43] +
        2 * d[15] * d[33] * d[43] + 6 * d[16] * d[34] * d[43] + 2 * d[17] * d[35] * d[43] + 2 * d[11] * d[28] * d[44] +
        2 * d[10] * d[29] * d[44] + 2 * d[14] * d[31] * d[44] + 2 * d[13] * d[32] * d[44] + 2 * d[17] * d[34] * d[44] +
        2 * d[16] * d[35] * d[44];
    coeffs[343] =
        -2 * d[25] * d[27] * d[36] + 2 * d[24] * d[28] * d[36] + 2 * d[19] * d[33] * d[36] - 2 * d[18] * d[34] * d[36] +
        2 * d[24] * d[27] * d[37] + 2 * d[25] * d[28] * d[37] + 2 * d[26] * d[29] * d[37] + 2 * d[18] * d[33] * d[37] +
        2 * d[19] * d[34] * d[37] + 2 * d[20] * d[35] * d[37] + 2 * d[26] * d[28] * d[38] - 2 * d[25] * d[29] * d[38] -
        2 * d[20] * d[34] * d[38] + 2 * d[19] * d[35] * d[38] - 2 * d[25] * d[30] * d[39] + 2 * d[24] * d[31] * d[39] +
        2 * d[22] * d[33] * d[39] - 2 * d[21] * d[34] * d[39] + 2 * d[24] * d[30] * d[40] + 2 * d[25] * d[31] * d[40] +
        2 * d[26] * d[32] * d[40] + 2 * d[21] * d[33] * d[40] + 2 * d[22] * d[34] * d[40] + 2 * d[23] * d[35] * d[40] +
        2 * d[26] * d[31] * d[41] - 2 * d[25] * d[32] * d[41] - 2 * d[23] * d[34] * d[41] + 2 * d[22] * d[35] * d[41] +
        2 * d[19] * d[27] * d[42] + 2 * d[18] * d[28] * d[42] + 2 * d[22] * d[30] * d[42] + 2 * d[21] * d[31] * d[42] +
        2 * d[25] * d[33] * d[42] + 2 * d[24] * d[34] * d[42] - 2 * d[18] * d[27] * d[43] + 2 * d[19] * d[28] * d[43] -
        2 * d[20] * d[29] * d[43] - 2 * d[21] * d[30] * d[43] + 2 * d[22] * d[31] * d[43] - 2 * d[23] * d[32] * d[43] +
        2 * d[24] * d[33] * d[43] + 6 * d[25] * d[34] * d[43] + 2 * d[26] * d[35] * d[43] + 2 * d[20] * d[28] * d[44] +
        2 * d[19] * d[29] * d[44] + 2 * d[23] * d[31] * d[44] + 2 * d[22] * d[32] * d[44] + 2 * d[26] * d[34] * d[44] +
        2 * d[25] * d[35] * d[44];
    coeffs[344] =
        2 * d[28] * d[33] * d[36] - 2 * d[27] * d[34] * d[36] + 2 * d[27] * d[33] * d[37] + 2 * d[28] * d[34] * d[37] +
        2 * d[29] * d[35] * d[37] - 2 * d[29] * d[34] * d[38] + 2 * d[28] * d[35] * d[38] + 2 * d[31] * d[33] * d[39] -
        2 * d[30] * d[34] * d[39] + 2 * d[30] * d[33] * d[40] + 2 * d[31] * d[34] * d[40] + 2 * d[32] * d[35] * d[40] -
        2 * d[32] * d[34] * d[41] + 2 * d[31] * d[35] * d[41] + 2 * d[27] * d[28] * d[42] + 2 * d[30] * d[31] * d[42] +
        2 * d[33] * d[34] * d[42] - std::pow(d[27], 2) * d[43] + std::pow(d[28], 2) * d[43] -
        std::pow(d[29], 2) * d[43] - std::pow(d[30], 2) * d[43] + std::pow(d[31], 2) * d[43] -
        std::pow(d[32], 2) * d[43] + std::pow(d[33], 2) * d[43] + 3 * std::pow(d[34], 2) * d[43] +
        std::pow(d[35], 2) * d[43] + 2 * d[28] * d[29] * d[44] + 2 * d[31] * d[32] * d[44] + 2 * d[34] * d[35] * d[44];
    coeffs[345] =
        -d[7] * std::pow(d[36], 2) + 2 * d[6] * d[36] * d[37] + d[7] * std::pow(d[37], 2) + 2 * d[8] * d[37] * d[38] -
        d[7] * std::pow(d[38], 2) - d[7] * std::pow(d[39], 2) + 2 * d[6] * d[39] * d[40] + d[7] * std::pow(d[40], 2) +
        2 * d[8] * d[40] * d[41] - d[7] * std::pow(d[41], 2) + 2 * d[1] * d[36] * d[42] + 2 * d[0] * d[37] * d[42] +
        2 * d[4] * d[39] * d[42] + 2 * d[3] * d[40] * d[42] + d[7] * std::pow(d[42], 2) - 2 * d[0] * d[36] * d[43] +
        2 * d[1] * d[37] * d[43] - 2 * d[2] * d[38] * d[43] - 2 * d[3] * d[39] * d[43] + 2 * d[4] * d[40] * d[43] -
        2 * d[5] * d[41] * d[43] + 2 * d[6] * d[42] * d[43] + 3 * d[7] * std::pow(d[43], 2) + 2 * d[2] * d[37] * d[44] +
        2 * d[1] * d[38] * d[44] + 2 * d[5] * d[40] * d[44] + 2 * d[4] * d[41] * d[44] + 2 * d[8] * d[43] * d[44] +
        d[7] * std::pow(d[44], 2);
    coeffs[346] = -d[16] * std::pow(d[36], 2) + 2 * d[15] * d[36] * d[37] + d[16] * std::pow(d[37], 2) +
                  2 * d[17] * d[37] * d[38] - d[16] * std::pow(d[38], 2) - d[16] * std::pow(d[39], 2) +
                  2 * d[15] * d[39] * d[40] + d[16] * std::pow(d[40], 2) + 2 * d[17] * d[40] * d[41] -
                  d[16] * std::pow(d[41], 2) + 2 * d[10] * d[36] * d[42] + 2 * d[9] * d[37] * d[42] +
                  2 * d[13] * d[39] * d[42] + 2 * d[12] * d[40] * d[42] + d[16] * std::pow(d[42], 2) -
                  2 * d[9] * d[36] * d[43] + 2 * d[10] * d[37] * d[43] - 2 * d[11] * d[38] * d[43] -
                  2 * d[12] * d[39] * d[43] + 2 * d[13] * d[40] * d[43] - 2 * d[14] * d[41] * d[43] +
                  2 * d[15] * d[42] * d[43] + 3 * d[16] * std::pow(d[43], 2) + 2 * d[11] * d[37] * d[44] +
                  2 * d[10] * d[38] * d[44] + 2 * d[14] * d[40] * d[44] + 2 * d[13] * d[41] * d[44] +
                  2 * d[17] * d[43] * d[44] + d[16] * std::pow(d[44], 2);
    coeffs[347] = -d[25] * std::pow(d[36], 2) + 2 * d[24] * d[36] * d[37] + d[25] * std::pow(d[37], 2) +
                  2 * d[26] * d[37] * d[38] - d[25] * std::pow(d[38], 2) - d[25] * std::pow(d[39], 2) +
                  2 * d[24] * d[39] * d[40] + d[25] * std::pow(d[40], 2) + 2 * d[26] * d[40] * d[41] -
                  d[25] * std::pow(d[41], 2) + 2 * d[19] * d[36] * d[42] + 2 * d[18] * d[37] * d[42] +
                  2 * d[22] * d[39] * d[42] + 2 * d[21] * d[40] * d[42] + d[25] * std::pow(d[42], 2) -
                  2 * d[18] * d[36] * d[43] + 2 * d[19] * d[37] * d[43] - 2 * d[20] * d[38] * d[43] -
                  2 * d[21] * d[39] * d[43] + 2 * d[22] * d[40] * d[43] - 2 * d[23] * d[41] * d[43] +
                  2 * d[24] * d[42] * d[43] + 3 * d[25] * std::pow(d[43], 2) + 2 * d[20] * d[37] * d[44] +
                  2 * d[19] * d[38] * d[44] + 2 * d[23] * d[40] * d[44] + 2 * d[22] * d[41] * d[44] +
                  2 * d[26] * d[43] * d[44] + d[25] * std::pow(d[44], 2);
    coeffs[348] = -d[34] * std::pow(d[36], 2) + 2 * d[33] * d[36] * d[37] + d[34] * std::pow(d[37], 2) +
                  2 * d[35] * d[37] * d[38] - d[34] * std::pow(d[38], 2) - d[34] * std::pow(d[39], 2) +
                  2 * d[33] * d[39] * d[40] + d[34] * std::pow(d[40], 2) + 2 * d[35] * d[40] * d[41] -
                  d[34] * std::pow(d[41], 2) + 2 * d[28] * d[36] * d[42] + 2 * d[27] * d[37] * d[42] +
                  2 * d[31] * d[39] * d[42] + 2 * d[30] * d[40] * d[42] + d[34] * std::pow(d[42], 2) -
                  2 * d[27] * d[36] * d[43] + 2 * d[28] * d[37] * d[43] - 2 * d[29] * d[38] * d[43] -
                  2 * d[30] * d[39] * d[43] + 2 * d[31] * d[40] * d[43] - 2 * d[32] * d[41] * d[43] +
                  2 * d[33] * d[42] * d[43] + 3 * d[34] * std::pow(d[43], 2) + 2 * d[29] * d[37] * d[44] +
                  2 * d[28] * d[38] * d[44] + 2 * d[32] * d[40] * d[44] + 2 * d[31] * d[41] * d[44] +
                  2 * d[35] * d[43] * d[44] + d[34] * std::pow(d[44], 2);
    coeffs[349] = 2 * d[36] * d[37] * d[42] + 2 * d[39] * d[40] * d[42] - std::pow(d[36], 2) * d[43] +
                  std::pow(d[37], 2) * d[43] - std::pow(d[38], 2) * d[43] - std::pow(d[39], 2) * d[43] +
                  std::pow(d[40], 2) * d[43] - std::pow(d[41], 2) * d[43] + std::pow(d[42], 2) * d[43] +
                  std::pow(d[43], 3) + 2 * d[37] * d[38] * d[44] + 2 * d[40] * d[41] * d[44] +
                  d[43] * std::pow(d[44], 2);
    coeffs[350] = 2 * d[0] * d[2] * d[6] + 2 * d[3] * d[5] * d[6] + 2 * d[1] * d[2] * d[7] + 2 * d[4] * d[5] * d[7] -
                  std::pow(d[0], 2) * d[8] - std::pow(d[1], 2) * d[8] + std::pow(d[2], 2) * d[8] -
                  std::pow(d[3], 2) * d[8] - std::pow(d[4], 2) * d[8] + std::pow(d[5], 2) * d[8] +
                  std::pow(d[6], 2) * d[8] + std::pow(d[7], 2) * d[8] + std::pow(d[8], 3);
    coeffs[351] =
        2 * d[2] * d[6] * d[9] - 2 * d[0] * d[8] * d[9] + 2 * d[2] * d[7] * d[10] - 2 * d[1] * d[8] * d[10] +
        2 * d[0] * d[6] * d[11] + 2 * d[1] * d[7] * d[11] + 2 * d[2] * d[8] * d[11] + 2 * d[5] * d[6] * d[12] -
        2 * d[3] * d[8] * d[12] + 2 * d[5] * d[7] * d[13] - 2 * d[4] * d[8] * d[13] + 2 * d[3] * d[6] * d[14] +
        2 * d[4] * d[7] * d[14] + 2 * d[5] * d[8] * d[14] + 2 * d[0] * d[2] * d[15] + 2 * d[3] * d[5] * d[15] +
        2 * d[6] * d[8] * d[15] + 2 * d[1] * d[2] * d[16] + 2 * d[4] * d[5] * d[16] + 2 * d[7] * d[8] * d[16] -
        std::pow(d[0], 2) * d[17] - std::pow(d[1], 2) * d[17] + std::pow(d[2], 2) * d[17] - std::pow(d[3], 2) * d[17] -
        std::pow(d[4], 2) * d[17] + std::pow(d[5], 2) * d[17] + std::pow(d[6], 2) * d[17] + std::pow(d[7], 2) * d[17] +
        3 * std::pow(d[8], 2) * d[17];
    coeffs[352] =
        -d[8] * std::pow(d[9], 2) - d[8] * std::pow(d[10], 2) + 2 * d[6] * d[9] * d[11] + 2 * d[7] * d[10] * d[11] +
        d[8] * std::pow(d[11], 2) - d[8] * std::pow(d[12], 2) - d[8] * std::pow(d[13], 2) + 2 * d[6] * d[12] * d[14] +
        2 * d[7] * d[13] * d[14] + d[8] * std::pow(d[14], 2) + 2 * d[2] * d[9] * d[15] + 2 * d[0] * d[11] * d[15] +
        2 * d[5] * d[12] * d[15] + 2 * d[3] * d[14] * d[15] + d[8] * std::pow(d[15], 2) + 2 * d[2] * d[10] * d[16] +
        2 * d[1] * d[11] * d[16] + 2 * d[5] * d[13] * d[16] + 2 * d[4] * d[14] * d[16] + d[8] * std::pow(d[16], 2) -
        2 * d[0] * d[9] * d[17] - 2 * d[1] * d[10] * d[17] + 2 * d[2] * d[11] * d[17] - 2 * d[3] * d[12] * d[17] -
        2 * d[4] * d[13] * d[17] + 2 * d[5] * d[14] * d[17] + 2 * d[6] * d[15] * d[17] + 2 * d[7] * d[16] * d[17] +
        3 * d[8] * std::pow(d[17], 2);
    coeffs[353] = 2 * d[9] * d[11] * d[15] + 2 * d[12] * d[14] * d[15] + 2 * d[10] * d[11] * d[16] +
                  2 * d[13] * d[14] * d[16] - std::pow(d[9], 2) * d[17] - std::pow(d[10], 2) * d[17] +
                  std::pow(d[11], 2) * d[17] - std::pow(d[12], 2) * d[17] - std::pow(d[13], 2) * d[17] +
                  std::pow(d[14], 2) * d[17] + std::pow(d[15], 2) * d[17] + std::pow(d[16], 2) * d[17] +
                  std::pow(d[17], 3);
    coeffs[354] =
        2 * d[2] * d[6] * d[18] - 2 * d[0] * d[8] * d[18] + 2 * d[2] * d[7] * d[19] - 2 * d[1] * d[8] * d[19] +
        2 * d[0] * d[6] * d[20] + 2 * d[1] * d[7] * d[20] + 2 * d[2] * d[8] * d[20] + 2 * d[5] * d[6] * d[21] -
        2 * d[3] * d[8] * d[21] + 2 * d[5] * d[7] * d[22] - 2 * d[4] * d[8] * d[22] + 2 * d[3] * d[6] * d[23] +
        2 * d[4] * d[7] * d[23] + 2 * d[5] * d[8] * d[23] + 2 * d[0] * d[2] * d[24] + 2 * d[3] * d[5] * d[24] +
        2 * d[6] * d[8] * d[24] + 2 * d[1] * d[2] * d[25] + 2 * d[4] * d[5] * d[25] + 2 * d[7] * d[8] * d[25] -
        std::pow(d[0], 2) * d[26] - std::pow(d[1], 2) * d[26] + std::pow(d[2], 2) * d[26] - std::pow(d[3], 2) * d[26] -
        std::pow(d[4], 2) * d[26] + std::pow(d[5], 2) * d[26] + std::pow(d[6], 2) * d[26] + std::pow(d[7], 2) * d[26] +
        3 * std::pow(d[8], 2) * d[26];
    coeffs[355] =
        -2 * d[8] * d[9] * d[18] + 2 * d[6] * d[11] * d[18] + 2 * d[2] * d[15] * d[18] - 2 * d[0] * d[17] * d[18] -
        2 * d[8] * d[10] * d[19] + 2 * d[7] * d[11] * d[19] + 2 * d[2] * d[16] * d[19] - 2 * d[1] * d[17] * d[19] +
        2 * d[6] * d[9] * d[20] + 2 * d[7] * d[10] * d[20] + 2 * d[8] * d[11] * d[20] + 2 * d[0] * d[15] * d[20] +
        2 * d[1] * d[16] * d[20] + 2 * d[2] * d[17] * d[20] - 2 * d[8] * d[12] * d[21] + 2 * d[6] * d[14] * d[21] +
        2 * d[5] * d[15] * d[21] - 2 * d[3] * d[17] * d[21] - 2 * d[8] * d[13] * d[22] + 2 * d[7] * d[14] * d[22] +
        2 * d[5] * d[16] * d[22] - 2 * d[4] * d[17] * d[22] + 2 * d[6] * d[12] * d[23] + 2 * d[7] * d[13] * d[23] +
        2 * d[8] * d[14] * d[23] + 2 * d[3] * d[15] * d[23] + 2 * d[4] * d[16] * d[23] + 2 * d[5] * d[17] * d[23] +
        2 * d[2] * d[9] * d[24] + 2 * d[0] * d[11] * d[24] + 2 * d[5] * d[12] * d[24] + 2 * d[3] * d[14] * d[24] +
        2 * d[8] * d[15] * d[24] + 2 * d[6] * d[17] * d[24] + 2 * d[2] * d[10] * d[25] + 2 * d[1] * d[11] * d[25] +
        2 * d[5] * d[13] * d[25] + 2 * d[4] * d[14] * d[25] + 2 * d[8] * d[16] * d[25] + 2 * d[7] * d[17] * d[25] -
        2 * d[0] * d[9] * d[26] - 2 * d[1] * d[10] * d[26] + 2 * d[2] * d[11] * d[26] - 2 * d[3] * d[12] * d[26] -
        2 * d[4] * d[13] * d[26] + 2 * d[5] * d[14] * d[26] + 2 * d[6] * d[15] * d[26] + 2 * d[7] * d[16] * d[26] +
        6 * d[8] * d[17] * d[26];
    coeffs[356] =
        2 * d[11] * d[15] * d[18] - 2 * d[9] * d[17] * d[18] + 2 * d[11] * d[16] * d[19] - 2 * d[10] * d[17] * d[19] +
        2 * d[9] * d[15] * d[20] + 2 * d[10] * d[16] * d[20] + 2 * d[11] * d[17] * d[20] + 2 * d[14] * d[15] * d[21] -
        2 * d[12] * d[17] * d[21] + 2 * d[14] * d[16] * d[22] - 2 * d[13] * d[17] * d[22] + 2 * d[12] * d[15] * d[23] +
        2 * d[13] * d[16] * d[23] + 2 * d[14] * d[17] * d[23] + 2 * d[9] * d[11] * d[24] + 2 * d[12] * d[14] * d[24] +
        2 * d[15] * d[17] * d[24] + 2 * d[10] * d[11] * d[25] + 2 * d[13] * d[14] * d[25] + 2 * d[16] * d[17] * d[25] -
        std::pow(d[9], 2) * d[26] - std::pow(d[10], 2) * d[26] + std::pow(d[11], 2) * d[26] -
        std::pow(d[12], 2) * d[26] - std::pow(d[13], 2) * d[26] + std::pow(d[14], 2) * d[26] +
        std::pow(d[15], 2) * d[26] + std::pow(d[16], 2) * d[26] + 3 * std::pow(d[17], 2) * d[26];
    coeffs[357] =
        -d[8] * std::pow(d[18], 2) - d[8] * std::pow(d[19], 2) + 2 * d[6] * d[18] * d[20] + 2 * d[7] * d[19] * d[20] +
        d[8] * std::pow(d[20], 2) - d[8] * std::pow(d[21], 2) - d[8] * std::pow(d[22], 2) + 2 * d[6] * d[21] * d[23] +
        2 * d[7] * d[22] * d[23] + d[8] * std::pow(d[23], 2) + 2 * d[2] * d[18] * d[24] + 2 * d[0] * d[20] * d[24] +
        2 * d[5] * d[21] * d[24] + 2 * d[3] * d[23] * d[24] + d[8] * std::pow(d[24], 2) + 2 * d[2] * d[19] * d[25] +
        2 * d[1] * d[20] * d[25] + 2 * d[5] * d[22] * d[25] + 2 * d[4] * d[23] * d[25] + d[8] * std::pow(d[25], 2) -
        2 * d[0] * d[18] * d[26] - 2 * d[1] * d[19] * d[26] + 2 * d[2] * d[20] * d[26] - 2 * d[3] * d[21] * d[26] -
        2 * d[4] * d[22] * d[26] + 2 * d[5] * d[23] * d[26] + 2 * d[6] * d[24] * d[26] + 2 * d[7] * d[25] * d[26] +
        3 * d[8] * std::pow(d[26], 2);
    coeffs[358] = -d[17] * std::pow(d[18], 2) - d[17] * std::pow(d[19], 2) + 2 * d[15] * d[18] * d[20] +
                  2 * d[16] * d[19] * d[20] + d[17] * std::pow(d[20], 2) - d[17] * std::pow(d[21], 2) -
                  d[17] * std::pow(d[22], 2) + 2 * d[15] * d[21] * d[23] + 2 * d[16] * d[22] * d[23] +
                  d[17] * std::pow(d[23], 2) + 2 * d[11] * d[18] * d[24] + 2 * d[9] * d[20] * d[24] +
                  2 * d[14] * d[21] * d[24] + 2 * d[12] * d[23] * d[24] + d[17] * std::pow(d[24], 2) +
                  2 * d[11] * d[19] * d[25] + 2 * d[10] * d[20] * d[25] + 2 * d[14] * d[22] * d[25] +
                  2 * d[13] * d[23] * d[25] + d[17] * std::pow(d[25], 2) - 2 * d[9] * d[18] * d[26] -
                  2 * d[10] * d[19] * d[26] + 2 * d[11] * d[20] * d[26] - 2 * d[12] * d[21] * d[26] -
                  2 * d[13] * d[22] * d[26] + 2 * d[14] * d[23] * d[26] + 2 * d[15] * d[24] * d[26] +
                  2 * d[16] * d[25] * d[26] + 3 * d[17] * std::pow(d[26], 2);
    coeffs[359] = 2 * d[18] * d[20] * d[24] + 2 * d[21] * d[23] * d[24] + 2 * d[19] * d[20] * d[25] +
                  2 * d[22] * d[23] * d[25] - std::pow(d[18], 2) * d[26] - std::pow(d[19], 2) * d[26] +
                  std::pow(d[20], 2) * d[26] - std::pow(d[21], 2) * d[26] - std::pow(d[22], 2) * d[26] +
                  std::pow(d[23], 2) * d[26] + std::pow(d[24], 2) * d[26] + std::pow(d[25], 2) * d[26] +
                  std::pow(d[26], 3);
    coeffs[360] =
        2 * d[2] * d[6] * d[27] - 2 * d[0] * d[8] * d[27] + 2 * d[2] * d[7] * d[28] - 2 * d[1] * d[8] * d[28] +
        2 * d[0] * d[6] * d[29] + 2 * d[1] * d[7] * d[29] + 2 * d[2] * d[8] * d[29] + 2 * d[5] * d[6] * d[30] -
        2 * d[3] * d[8] * d[30] + 2 * d[5] * d[7] * d[31] - 2 * d[4] * d[8] * d[31] + 2 * d[3] * d[6] * d[32] +
        2 * d[4] * d[7] * d[32] + 2 * d[5] * d[8] * d[32] + 2 * d[0] * d[2] * d[33] + 2 * d[3] * d[5] * d[33] +
        2 * d[6] * d[8] * d[33] + 2 * d[1] * d[2] * d[34] + 2 * d[4] * d[5] * d[34] + 2 * d[7] * d[8] * d[34] -
        std::pow(d[0], 2) * d[35] - std::pow(d[1], 2) * d[35] + std::pow(d[2], 2) * d[35] - std::pow(d[3], 2) * d[35] -
        std::pow(d[4], 2) * d[35] + std::pow(d[5], 2) * d[35] + std::pow(d[6], 2) * d[35] + std::pow(d[7], 2) * d[35] +
        3 * std::pow(d[8], 2) * d[35];
    coeffs[361] =
        -2 * d[8] * d[9] * d[27] + 2 * d[6] * d[11] * d[27] + 2 * d[2] * d[15] * d[27] - 2 * d[0] * d[17] * d[27] -
        2 * d[8] * d[10] * d[28] + 2 * d[7] * d[11] * d[28] + 2 * d[2] * d[16] * d[28] - 2 * d[1] * d[17] * d[28] +
        2 * d[6] * d[9] * d[29] + 2 * d[7] * d[10] * d[29] + 2 * d[8] * d[11] * d[29] + 2 * d[0] * d[15] * d[29] +
        2 * d[1] * d[16] * d[29] + 2 * d[2] * d[17] * d[29] - 2 * d[8] * d[12] * d[30] + 2 * d[6] * d[14] * d[30] +
        2 * d[5] * d[15] * d[30] - 2 * d[3] * d[17] * d[30] - 2 * d[8] * d[13] * d[31] + 2 * d[7] * d[14] * d[31] +
        2 * d[5] * d[16] * d[31] - 2 * d[4] * d[17] * d[31] + 2 * d[6] * d[12] * d[32] + 2 * d[7] * d[13] * d[32] +
        2 * d[8] * d[14] * d[32] + 2 * d[3] * d[15] * d[32] + 2 * d[4] * d[16] * d[32] + 2 * d[5] * d[17] * d[32] +
        2 * d[2] * d[9] * d[33] + 2 * d[0] * d[11] * d[33] + 2 * d[5] * d[12] * d[33] + 2 * d[3] * d[14] * d[33] +
        2 * d[8] * d[15] * d[33] + 2 * d[6] * d[17] * d[33] + 2 * d[2] * d[10] * d[34] + 2 * d[1] * d[11] * d[34] +
        2 * d[5] * d[13] * d[34] + 2 * d[4] * d[14] * d[34] + 2 * d[8] * d[16] * d[34] + 2 * d[7] * d[17] * d[34] -
        2 * d[0] * d[9] * d[35] - 2 * d[1] * d[10] * d[35] + 2 * d[2] * d[11] * d[35] - 2 * d[3] * d[12] * d[35] -
        2 * d[4] * d[13] * d[35] + 2 * d[5] * d[14] * d[35] + 2 * d[6] * d[15] * d[35] + 2 * d[7] * d[16] * d[35] +
        6 * d[8] * d[17] * d[35];
    coeffs[362] =
        2 * d[11] * d[15] * d[27] - 2 * d[9] * d[17] * d[27] + 2 * d[11] * d[16] * d[28] - 2 * d[10] * d[17] * d[28] +
        2 * d[9] * d[15] * d[29] + 2 * d[10] * d[16] * d[29] + 2 * d[11] * d[17] * d[29] + 2 * d[14] * d[15] * d[30] -
        2 * d[12] * d[17] * d[30] + 2 * d[14] * d[16] * d[31] - 2 * d[13] * d[17] * d[31] + 2 * d[12] * d[15] * d[32] +
        2 * d[13] * d[16] * d[32] + 2 * d[14] * d[17] * d[32] + 2 * d[9] * d[11] * d[33] + 2 * d[12] * d[14] * d[33] +
        2 * d[15] * d[17] * d[33] + 2 * d[10] * d[11] * d[34] + 2 * d[13] * d[14] * d[34] + 2 * d[16] * d[17] * d[34] -
        std::pow(d[9], 2) * d[35] - std::pow(d[10], 2) * d[35] + std::pow(d[11], 2) * d[35] -
        std::pow(d[12], 2) * d[35] - std::pow(d[13], 2) * d[35] + std::pow(d[14], 2) * d[35] +
        std::pow(d[15], 2) * d[35] + std::pow(d[16], 2) * d[35] + 3 * std::pow(d[17], 2) * d[35];
    coeffs[363] =
        -2 * d[8] * d[18] * d[27] + 2 * d[6] * d[20] * d[27] + 2 * d[2] * d[24] * d[27] - 2 * d[0] * d[26] * d[27] -
        2 * d[8] * d[19] * d[28] + 2 * d[7] * d[20] * d[28] + 2 * d[2] * d[25] * d[28] - 2 * d[1] * d[26] * d[28] +
        2 * d[6] * d[18] * d[29] + 2 * d[7] * d[19] * d[29] + 2 * d[8] * d[20] * d[29] + 2 * d[0] * d[24] * d[29] +
        2 * d[1] * d[25] * d[29] + 2 * d[2] * d[26] * d[29] - 2 * d[8] * d[21] * d[30] + 2 * d[6] * d[23] * d[30] +
        2 * d[5] * d[24] * d[30] - 2 * d[3] * d[26] * d[30] - 2 * d[8] * d[22] * d[31] + 2 * d[7] * d[23] * d[31] +
        2 * d[5] * d[25] * d[31] - 2 * d[4] * d[26] * d[31] + 2 * d[6] * d[21] * d[32] + 2 * d[7] * d[22] * d[32] +
        2 * d[8] * d[23] * d[32] + 2 * d[3] * d[24] * d[32] + 2 * d[4] * d[25] * d[32] + 2 * d[5] * d[26] * d[32] +
        2 * d[2] * d[18] * d[33] + 2 * d[0] * d[20] * d[33] + 2 * d[5] * d[21] * d[33] + 2 * d[3] * d[23] * d[33] +
        2 * d[8] * d[24] * d[33] + 2 * d[6] * d[26] * d[33] + 2 * d[2] * d[19] * d[34] + 2 * d[1] * d[20] * d[34] +
        2 * d[5] * d[22] * d[34] + 2 * d[4] * d[23] * d[34] + 2 * d[8] * d[25] * d[34] + 2 * d[7] * d[26] * d[34] -
        2 * d[0] * d[18] * d[35] - 2 * d[1] * d[19] * d[35] + 2 * d[2] * d[20] * d[35] - 2 * d[3] * d[21] * d[35] -
        2 * d[4] * d[22] * d[35] + 2 * d[5] * d[23] * d[35] + 2 * d[6] * d[24] * d[35] + 2 * d[7] * d[25] * d[35] +
        6 * d[8] * d[26] * d[35];
    coeffs[364] =
        -2 * d[17] * d[18] * d[27] + 2 * d[15] * d[20] * d[27] + 2 * d[11] * d[24] * d[27] - 2 * d[9] * d[26] * d[27] -
        2 * d[17] * d[19] * d[28] + 2 * d[16] * d[20] * d[28] + 2 * d[11] * d[25] * d[28] - 2 * d[10] * d[26] * d[28] +
        2 * d[15] * d[18] * d[29] + 2 * d[16] * d[19] * d[29] + 2 * d[17] * d[20] * d[29] + 2 * d[9] * d[24] * d[29] +
        2 * d[10] * d[25] * d[29] + 2 * d[11] * d[26] * d[29] - 2 * d[17] * d[21] * d[30] + 2 * d[15] * d[23] * d[30] +
        2 * d[14] * d[24] * d[30] - 2 * d[12] * d[26] * d[30] - 2 * d[17] * d[22] * d[31] + 2 * d[16] * d[23] * d[31] +
        2 * d[14] * d[25] * d[31] - 2 * d[13] * d[26] * d[31] + 2 * d[15] * d[21] * d[32] + 2 * d[16] * d[22] * d[32] +
        2 * d[17] * d[23] * d[32] + 2 * d[12] * d[24] * d[32] + 2 * d[13] * d[25] * d[32] + 2 * d[14] * d[26] * d[32] +
        2 * d[11] * d[18] * d[33] + 2 * d[9] * d[20] * d[33] + 2 * d[14] * d[21] * d[33] + 2 * d[12] * d[23] * d[33] +
        2 * d[17] * d[24] * d[33] + 2 * d[15] * d[26] * d[33] + 2 * d[11] * d[19] * d[34] + 2 * d[10] * d[20] * d[34] +
        2 * d[14] * d[22] * d[34] + 2 * d[13] * d[23] * d[34] + 2 * d[17] * d[25] * d[34] + 2 * d[16] * d[26] * d[34] -
        2 * d[9] * d[18] * d[35] - 2 * d[10] * d[19] * d[35] + 2 * d[11] * d[20] * d[35] - 2 * d[12] * d[21] * d[35] -
        2 * d[13] * d[22] * d[35] + 2 * d[14] * d[23] * d[35] + 2 * d[15] * d[24] * d[35] + 2 * d[16] * d[25] * d[35] +
        6 * d[17] * d[26] * d[35];
    coeffs[365] =
        2 * d[20] * d[24] * d[27] - 2 * d[18] * d[26] * d[27] + 2 * d[20] * d[25] * d[28] - 2 * d[19] * d[26] * d[28] +
        2 * d[18] * d[24] * d[29] + 2 * d[19] * d[25] * d[29] + 2 * d[20] * d[26] * d[29] + 2 * d[23] * d[24] * d[30] -
        2 * d[21] * d[26] * d[30] + 2 * d[23] * d[25] * d[31] - 2 * d[22] * d[26] * d[31] + 2 * d[21] * d[24] * d[32] +
        2 * d[22] * d[25] * d[32] + 2 * d[23] * d[26] * d[32] + 2 * d[18] * d[20] * d[33] + 2 * d[21] * d[23] * d[33] +
        2 * d[24] * d[26] * d[33] + 2 * d[19] * d[20] * d[34] + 2 * d[22] * d[23] * d[34] + 2 * d[25] * d[26] * d[34] -
        std::pow(d[18], 2) * d[35] - std::pow(d[19], 2) * d[35] + std::pow(d[20], 2) * d[35] -
        std::pow(d[21], 2) * d[35] - std::pow(d[22], 2) * d[35] + std::pow(d[23], 2) * d[35] +
        std::pow(d[24], 2) * d[35] + std::pow(d[25], 2) * d[35] + 3 * std::pow(d[26], 2) * d[35];
    coeffs[366] =
        -d[8] * std::pow(d[27], 2) - d[8] * std::pow(d[28], 2) + 2 * d[6] * d[27] * d[29] + 2 * d[7] * d[28] * d[29] +
        d[8] * std::pow(d[29], 2) - d[8] * std::pow(d[30], 2) - d[8] * std::pow(d[31], 2) + 2 * d[6] * d[30] * d[32] +
        2 * d[7] * d[31] * d[32] + d[8] * std::pow(d[32], 2) + 2 * d[2] * d[27] * d[33] + 2 * d[0] * d[29] * d[33] +
        2 * d[5] * d[30] * d[33] + 2 * d[3] * d[32] * d[33] + d[8] * std::pow(d[33], 2) + 2 * d[2] * d[28] * d[34] +
        2 * d[1] * d[29] * d[34] + 2 * d[5] * d[31] * d[34] + 2 * d[4] * d[32] * d[34] + d[8] * std::pow(d[34], 2) -
        2 * d[0] * d[27] * d[35] - 2 * d[1] * d[28] * d[35] + 2 * d[2] * d[29] * d[35] - 2 * d[3] * d[30] * d[35] -
        2 * d[4] * d[31] * d[35] + 2 * d[5] * d[32] * d[35] + 2 * d[6] * d[33] * d[35] + 2 * d[7] * d[34] * d[35] +
        3 * d[8] * std::pow(d[35], 2);
    coeffs[367] = -d[17] * std::pow(d[27], 2) - d[17] * std::pow(d[28], 2) + 2 * d[15] * d[27] * d[29] +
                  2 * d[16] * d[28] * d[29] + d[17] * std::pow(d[29], 2) - d[17] * std::pow(d[30], 2) -
                  d[17] * std::pow(d[31], 2) + 2 * d[15] * d[30] * d[32] + 2 * d[16] * d[31] * d[32] +
                  d[17] * std::pow(d[32], 2) + 2 * d[11] * d[27] * d[33] + 2 * d[9] * d[29] * d[33] +
                  2 * d[14] * d[30] * d[33] + 2 * d[12] * d[32] * d[33] + d[17] * std::pow(d[33], 2) +
                  2 * d[11] * d[28] * d[34] + 2 * d[10] * d[29] * d[34] + 2 * d[14] * d[31] * d[34] +
                  2 * d[13] * d[32] * d[34] + d[17] * std::pow(d[34], 2) - 2 * d[9] * d[27] * d[35] -
                  2 * d[10] * d[28] * d[35] + 2 * d[11] * d[29] * d[35] - 2 * d[12] * d[30] * d[35] -
                  2 * d[13] * d[31] * d[35] + 2 * d[14] * d[32] * d[35] + 2 * d[15] * d[33] * d[35] +
                  2 * d[16] * d[34] * d[35] + 3 * d[17] * std::pow(d[35], 2);
    coeffs[368] = -d[26] * std::pow(d[27], 2) - d[26] * std::pow(d[28], 2) + 2 * d[24] * d[27] * d[29] +
                  2 * d[25] * d[28] * d[29] + d[26] * std::pow(d[29], 2) - d[26] * std::pow(d[30], 2) -
                  d[26] * std::pow(d[31], 2) + 2 * d[24] * d[30] * d[32] + 2 * d[25] * d[31] * d[32] +
                  d[26] * std::pow(d[32], 2) + 2 * d[20] * d[27] * d[33] + 2 * d[18] * d[29] * d[33] +
                  2 * d[23] * d[30] * d[33] + 2 * d[21] * d[32] * d[33] + d[26] * std::pow(d[33], 2) +
                  2 * d[20] * d[28] * d[34] + 2 * d[19] * d[29] * d[34] + 2 * d[23] * d[31] * d[34] +
                  2 * d[22] * d[32] * d[34] + d[26] * std::pow(d[34], 2) - 2 * d[18] * d[27] * d[35] -
                  2 * d[19] * d[28] * d[35] + 2 * d[20] * d[29] * d[35] - 2 * d[21] * d[30] * d[35] -
                  2 * d[22] * d[31] * d[35] + 2 * d[23] * d[32] * d[35] + 2 * d[24] * d[33] * d[35] +
                  2 * d[25] * d[34] * d[35] + 3 * d[26] * std::pow(d[35], 2);
    coeffs[369] = 2 * d[27] * d[29] * d[33] + 2 * d[30] * d[32] * d[33] + 2 * d[28] * d[29] * d[34] +
                  2 * d[31] * d[32] * d[34] - std::pow(d[27], 2) * d[35] - std::pow(d[28], 2) * d[35] +
                  std::pow(d[29], 2) * d[35] - std::pow(d[30], 2) * d[35] - std::pow(d[31], 2) * d[35] +
                  std::pow(d[32], 2) * d[35] + std::pow(d[33], 2) * d[35] + std::pow(d[34], 2) * d[35] +
                  std::pow(d[35], 3);
    coeffs[370] =
        2 * d[2] * d[6] * d[36] - 2 * d[0] * d[8] * d[36] + 2 * d[2] * d[7] * d[37] - 2 * d[1] * d[8] * d[37] +
        2 * d[0] * d[6] * d[38] + 2 * d[1] * d[7] * d[38] + 2 * d[2] * d[8] * d[38] + 2 * d[5] * d[6] * d[39] -
        2 * d[3] * d[8] * d[39] + 2 * d[5] * d[7] * d[40] - 2 * d[4] * d[8] * d[40] + 2 * d[3] * d[6] * d[41] +
        2 * d[4] * d[7] * d[41] + 2 * d[5] * d[8] * d[41] + 2 * d[0] * d[2] * d[42] + 2 * d[3] * d[5] * d[42] +
        2 * d[6] * d[8] * d[42] + 2 * d[1] * d[2] * d[43] + 2 * d[4] * d[5] * d[43] + 2 * d[7] * d[8] * d[43] -
        std::pow(d[0], 2) * d[44] - std::pow(d[1], 2) * d[44] + std::pow(d[2], 2) * d[44] - std::pow(d[3], 2) * d[44] -
        std::pow(d[4], 2) * d[44] + std::pow(d[5], 2) * d[44] + std::pow(d[6], 2) * d[44] + std::pow(d[7], 2) * d[44] +
        3 * std::pow(d[8], 2) * d[44];
    coeffs[371] =
        -2 * d[8] * d[9] * d[36] + 2 * d[6] * d[11] * d[36] + 2 * d[2] * d[15] * d[36] - 2 * d[0] * d[17] * d[36] -
        2 * d[8] * d[10] * d[37] + 2 * d[7] * d[11] * d[37] + 2 * d[2] * d[16] * d[37] - 2 * d[1] * d[17] * d[37] +
        2 * d[6] * d[9] * d[38] + 2 * d[7] * d[10] * d[38] + 2 * d[8] * d[11] * d[38] + 2 * d[0] * d[15] * d[38] +
        2 * d[1] * d[16] * d[38] + 2 * d[2] * d[17] * d[38] - 2 * d[8] * d[12] * d[39] + 2 * d[6] * d[14] * d[39] +
        2 * d[5] * d[15] * d[39] - 2 * d[3] * d[17] * d[39] - 2 * d[8] * d[13] * d[40] + 2 * d[7] * d[14] * d[40] +
        2 * d[5] * d[16] * d[40] - 2 * d[4] * d[17] * d[40] + 2 * d[6] * d[12] * d[41] + 2 * d[7] * d[13] * d[41] +
        2 * d[8] * d[14] * d[41] + 2 * d[3] * d[15] * d[41] + 2 * d[4] * d[16] * d[41] + 2 * d[5] * d[17] * d[41] +
        2 * d[2] * d[9] * d[42] + 2 * d[0] * d[11] * d[42] + 2 * d[5] * d[12] * d[42] + 2 * d[3] * d[14] * d[42] +
        2 * d[8] * d[15] * d[42] + 2 * d[6] * d[17] * d[42] + 2 * d[2] * d[10] * d[43] + 2 * d[1] * d[11] * d[43] +
        2 * d[5] * d[13] * d[43] + 2 * d[4] * d[14] * d[43] + 2 * d[8] * d[16] * d[43] + 2 * d[7] * d[17] * d[43] -
        2 * d[0] * d[9] * d[44] - 2 * d[1] * d[10] * d[44] + 2 * d[2] * d[11] * d[44] - 2 * d[3] * d[12] * d[44] -
        2 * d[4] * d[13] * d[44] + 2 * d[5] * d[14] * d[44] + 2 * d[6] * d[15] * d[44] + 2 * d[7] * d[16] * d[44] +
        6 * d[8] * d[17] * d[44];
    coeffs[372] =
        2 * d[11] * d[15] * d[36] - 2 * d[9] * d[17] * d[36] + 2 * d[11] * d[16] * d[37] - 2 * d[10] * d[17] * d[37] +
        2 * d[9] * d[15] * d[38] + 2 * d[10] * d[16] * d[38] + 2 * d[11] * d[17] * d[38] + 2 * d[14] * d[15] * d[39] -
        2 * d[12] * d[17] * d[39] + 2 * d[14] * d[16] * d[40] - 2 * d[13] * d[17] * d[40] + 2 * d[12] * d[15] * d[41] +
        2 * d[13] * d[16] * d[41] + 2 * d[14] * d[17] * d[41] + 2 * d[9] * d[11] * d[42] + 2 * d[12] * d[14] * d[42] +
        2 * d[15] * d[17] * d[42] + 2 * d[10] * d[11] * d[43] + 2 * d[13] * d[14] * d[43] + 2 * d[16] * d[17] * d[43] -
        std::pow(d[9], 2) * d[44] - std::pow(d[10], 2) * d[44] + std::pow(d[11], 2) * d[44] -
        std::pow(d[12], 2) * d[44] - std::pow(d[13], 2) * d[44] + std::pow(d[14], 2) * d[44] +
        std::pow(d[15], 2) * d[44] + std::pow(d[16], 2) * d[44] + 3 * std::pow(d[17], 2) * d[44];
    coeffs[373] =
        -2 * d[8] * d[18] * d[36] + 2 * d[6] * d[20] * d[36] + 2 * d[2] * d[24] * d[36] - 2 * d[0] * d[26] * d[36] -
        2 * d[8] * d[19] * d[37] + 2 * d[7] * d[20] * d[37] + 2 * d[2] * d[25] * d[37] - 2 * d[1] * d[26] * d[37] +
        2 * d[6] * d[18] * d[38] + 2 * d[7] * d[19] * d[38] + 2 * d[8] * d[20] * d[38] + 2 * d[0] * d[24] * d[38] +
        2 * d[1] * d[25] * d[38] + 2 * d[2] * d[26] * d[38] - 2 * d[8] * d[21] * d[39] + 2 * d[6] * d[23] * d[39] +
        2 * d[5] * d[24] * d[39] - 2 * d[3] * d[26] * d[39] - 2 * d[8] * d[22] * d[40] + 2 * d[7] * d[23] * d[40] +
        2 * d[5] * d[25] * d[40] - 2 * d[4] * d[26] * d[40] + 2 * d[6] * d[21] * d[41] + 2 * d[7] * d[22] * d[41] +
        2 * d[8] * d[23] * d[41] + 2 * d[3] * d[24] * d[41] + 2 * d[4] * d[25] * d[41] + 2 * d[5] * d[26] * d[41] +
        2 * d[2] * d[18] * d[42] + 2 * d[0] * d[20] * d[42] + 2 * d[5] * d[21] * d[42] + 2 * d[3] * d[23] * d[42] +
        2 * d[8] * d[24] * d[42] + 2 * d[6] * d[26] * d[42] + 2 * d[2] * d[19] * d[43] + 2 * d[1] * d[20] * d[43] +
        2 * d[5] * d[22] * d[43] + 2 * d[4] * d[23] * d[43] + 2 * d[8] * d[25] * d[43] + 2 * d[7] * d[26] * d[43] -
        2 * d[0] * d[18] * d[44] - 2 * d[1] * d[19] * d[44] + 2 * d[2] * d[20] * d[44] - 2 * d[3] * d[21] * d[44] -
        2 * d[4] * d[22] * d[44] + 2 * d[5] * d[23] * d[44] + 2 * d[6] * d[24] * d[44] + 2 * d[7] * d[25] * d[44] +
        6 * d[8] * d[26] * d[44];
    coeffs[374] =
        -2 * d[17] * d[18] * d[36] + 2 * d[15] * d[20] * d[36] + 2 * d[11] * d[24] * d[36] - 2 * d[9] * d[26] * d[36] -
        2 * d[17] * d[19] * d[37] + 2 * d[16] * d[20] * d[37] + 2 * d[11] * d[25] * d[37] - 2 * d[10] * d[26] * d[37] +
        2 * d[15] * d[18] * d[38] + 2 * d[16] * d[19] * d[38] + 2 * d[17] * d[20] * d[38] + 2 * d[9] * d[24] * d[38] +
        2 * d[10] * d[25] * d[38] + 2 * d[11] * d[26] * d[38] - 2 * d[17] * d[21] * d[39] + 2 * d[15] * d[23] * d[39] +
        2 * d[14] * d[24] * d[39] - 2 * d[12] * d[26] * d[39] - 2 * d[17] * d[22] * d[40] + 2 * d[16] * d[23] * d[40] +
        2 * d[14] * d[25] * d[40] - 2 * d[13] * d[26] * d[40] + 2 * d[15] * d[21] * d[41] + 2 * d[16] * d[22] * d[41] +
        2 * d[17] * d[23] * d[41] + 2 * d[12] * d[24] * d[41] + 2 * d[13] * d[25] * d[41] + 2 * d[14] * d[26] * d[41] +
        2 * d[11] * d[18] * d[42] + 2 * d[9] * d[20] * d[42] + 2 * d[14] * d[21] * d[42] + 2 * d[12] * d[23] * d[42] +
        2 * d[17] * d[24] * d[42] + 2 * d[15] * d[26] * d[42] + 2 * d[11] * d[19] * d[43] + 2 * d[10] * d[20] * d[43] +
        2 * d[14] * d[22] * d[43] + 2 * d[13] * d[23] * d[43] + 2 * d[17] * d[25] * d[43] + 2 * d[16] * d[26] * d[43] -
        2 * d[9] * d[18] * d[44] - 2 * d[10] * d[19] * d[44] + 2 * d[11] * d[20] * d[44] - 2 * d[12] * d[21] * d[44] -
        2 * d[13] * d[22] * d[44] + 2 * d[14] * d[23] * d[44] + 2 * d[15] * d[24] * d[44] + 2 * d[16] * d[25] * d[44] +
        6 * d[17] * d[26] * d[44];
    coeffs[375] =
        2 * d[20] * d[24] * d[36] - 2 * d[18] * d[26] * d[36] + 2 * d[20] * d[25] * d[37] - 2 * d[19] * d[26] * d[37] +
        2 * d[18] * d[24] * d[38] + 2 * d[19] * d[25] * d[38] + 2 * d[20] * d[26] * d[38] + 2 * d[23] * d[24] * d[39] -
        2 * d[21] * d[26] * d[39] + 2 * d[23] * d[25] * d[40] - 2 * d[22] * d[26] * d[40] + 2 * d[21] * d[24] * d[41] +
        2 * d[22] * d[25] * d[41] + 2 * d[23] * d[26] * d[41] + 2 * d[18] * d[20] * d[42] + 2 * d[21] * d[23] * d[42] +
        2 * d[24] * d[26] * d[42] + 2 * d[19] * d[20] * d[43] + 2 * d[22] * d[23] * d[43] + 2 * d[25] * d[26] * d[43] -
        std::pow(d[18], 2) * d[44] - std::pow(d[19], 2) * d[44] + std::pow(d[20], 2) * d[44] -
        std::pow(d[21], 2) * d[44] - std::pow(d[22], 2) * d[44] + std::pow(d[23], 2) * d[44] +
        std::pow(d[24], 2) * d[44] + std::pow(d[25], 2) * d[44] + 3 * std::pow(d[26], 2) * d[44];
    coeffs[376] =
        -2 * d[8] * d[27] * d[36] + 2 * d[6] * d[29] * d[36] + 2 * d[2] * d[33] * d[36] - 2 * d[0] * d[35] * d[36] -
        2 * d[8] * d[28] * d[37] + 2 * d[7] * d[29] * d[37] + 2 * d[2] * d[34] * d[37] - 2 * d[1] * d[35] * d[37] +
        2 * d[6] * d[27] * d[38] + 2 * d[7] * d[28] * d[38] + 2 * d[8] * d[29] * d[38] + 2 * d[0] * d[33] * d[38] +
        2 * d[1] * d[34] * d[38] + 2 * d[2] * d[35] * d[38] - 2 * d[8] * d[30] * d[39] + 2 * d[6] * d[32] * d[39] +
        2 * d[5] * d[33] * d[39] - 2 * d[3] * d[35] * d[39] - 2 * d[8] * d[31] * d[40] + 2 * d[7] * d[32] * d[40] +
        2 * d[5] * d[34] * d[40] - 2 * d[4] * d[35] * d[40] + 2 * d[6] * d[30] * d[41] + 2 * d[7] * d[31] * d[41] +
        2 * d[8] * d[32] * d[41] + 2 * d[3] * d[33] * d[41] + 2 * d[4] * d[34] * d[41] + 2 * d[5] * d[35] * d[41] +
        2 * d[2] * d[27] * d[42] + 2 * d[0] * d[29] * d[42] + 2 * d[5] * d[30] * d[42] + 2 * d[3] * d[32] * d[42] +
        2 * d[8] * d[33] * d[42] + 2 * d[6] * d[35] * d[42] + 2 * d[2] * d[28] * d[43] + 2 * d[1] * d[29] * d[43] +
        2 * d[5] * d[31] * d[43] + 2 * d[4] * d[32] * d[43] + 2 * d[8] * d[34] * d[43] + 2 * d[7] * d[35] * d[43] -
        2 * d[0] * d[27] * d[44] - 2 * d[1] * d[28] * d[44] + 2 * d[2] * d[29] * d[44] - 2 * d[3] * d[30] * d[44] -
        2 * d[4] * d[31] * d[44] + 2 * d[5] * d[32] * d[44] + 2 * d[6] * d[33] * d[44] + 2 * d[7] * d[34] * d[44] +
        6 * d[8] * d[35] * d[44];
    coeffs[377] =
        -2 * d[17] * d[27] * d[36] + 2 * d[15] * d[29] * d[36] + 2 * d[11] * d[33] * d[36] - 2 * d[9] * d[35] * d[36] -
        2 * d[17] * d[28] * d[37] + 2 * d[16] * d[29] * d[37] + 2 * d[11] * d[34] * d[37] - 2 * d[10] * d[35] * d[37] +
        2 * d[15] * d[27] * d[38] + 2 * d[16] * d[28] * d[38] + 2 * d[17] * d[29] * d[38] + 2 * d[9] * d[33] * d[38] +
        2 * d[10] * d[34] * d[38] + 2 * d[11] * d[35] * d[38] - 2 * d[17] * d[30] * d[39] + 2 * d[15] * d[32] * d[39] +
        2 * d[14] * d[33] * d[39] - 2 * d[12] * d[35] * d[39] - 2 * d[17] * d[31] * d[40] + 2 * d[16] * d[32] * d[40] +
        2 * d[14] * d[34] * d[40] - 2 * d[13] * d[35] * d[40] + 2 * d[15] * d[30] * d[41] + 2 * d[16] * d[31] * d[41] +
        2 * d[17] * d[32] * d[41] + 2 * d[12] * d[33] * d[41] + 2 * d[13] * d[34] * d[41] + 2 * d[14] * d[35] * d[41] +
        2 * d[11] * d[27] * d[42] + 2 * d[9] * d[29] * d[42] + 2 * d[14] * d[30] * d[42] + 2 * d[12] * d[32] * d[42] +
        2 * d[17] * d[33] * d[42] + 2 * d[15] * d[35] * d[42] + 2 * d[11] * d[28] * d[43] + 2 * d[10] * d[29] * d[43] +
        2 * d[14] * d[31] * d[43] + 2 * d[13] * d[32] * d[43] + 2 * d[17] * d[34] * d[43] + 2 * d[16] * d[35] * d[43] -
        2 * d[9] * d[27] * d[44] - 2 * d[10] * d[28] * d[44] + 2 * d[11] * d[29] * d[44] - 2 * d[12] * d[30] * d[44] -
        2 * d[13] * d[31] * d[44] + 2 * d[14] * d[32] * d[44] + 2 * d[15] * d[33] * d[44] + 2 * d[16] * d[34] * d[44] +
        6 * d[17] * d[35] * d[44];
    coeffs[378] =
        -2 * d[26] * d[27] * d[36] + 2 * d[24] * d[29] * d[36] + 2 * d[20] * d[33] * d[36] - 2 * d[18] * d[35] * d[36] -
        2 * d[26] * d[28] * d[37] + 2 * d[25] * d[29] * d[37] + 2 * d[20] * d[34] * d[37] - 2 * d[19] * d[35] * d[37] +
        2 * d[24] * d[27] * d[38] + 2 * d[25] * d[28] * d[38] + 2 * d[26] * d[29] * d[38] + 2 * d[18] * d[33] * d[38] +
        2 * d[19] * d[34] * d[38] + 2 * d[20] * d[35] * d[38] - 2 * d[26] * d[30] * d[39] + 2 * d[24] * d[32] * d[39] +
        2 * d[23] * d[33] * d[39] - 2 * d[21] * d[35] * d[39] - 2 * d[26] * d[31] * d[40] + 2 * d[25] * d[32] * d[40] +
        2 * d[23] * d[34] * d[40] - 2 * d[22] * d[35] * d[40] + 2 * d[24] * d[30] * d[41] + 2 * d[25] * d[31] * d[41] +
        2 * d[26] * d[32] * d[41] + 2 * d[21] * d[33] * d[41] + 2 * d[22] * d[34] * d[41] + 2 * d[23] * d[35] * d[41] +
        2 * d[20] * d[27] * d[42] + 2 * d[18] * d[29] * d[42] + 2 * d[23] * d[30] * d[42] + 2 * d[21] * d[32] * d[42] +
        2 * d[26] * d[33] * d[42] + 2 * d[24] * d[35] * d[42] + 2 * d[20] * d[28] * d[43] + 2 * d[19] * d[29] * d[43] +
        2 * d[23] * d[31] * d[43] + 2 * d[22] * d[32] * d[43] + 2 * d[26] * d[34] * d[43] + 2 * d[25] * d[35] * d[43] -
        2 * d[18] * d[27] * d[44] - 2 * d[19] * d[28] * d[44] + 2 * d[20] * d[29] * d[44] - 2 * d[21] * d[30] * d[44] -
        2 * d[22] * d[31] * d[44] + 2 * d[23] * d[32] * d[44] + 2 * d[24] * d[33] * d[44] + 2 * d[25] * d[34] * d[44] +
        6 * d[26] * d[35] * d[44];
    coeffs[379] =
        2 * d[29] * d[33] * d[36] - 2 * d[27] * d[35] * d[36] + 2 * d[29] * d[34] * d[37] - 2 * d[28] * d[35] * d[37] +
        2 * d[27] * d[33] * d[38] + 2 * d[28] * d[34] * d[38] + 2 * d[29] * d[35] * d[38] + 2 * d[32] * d[33] * d[39] -
        2 * d[30] * d[35] * d[39] + 2 * d[32] * d[34] * d[40] - 2 * d[31] * d[35] * d[40] + 2 * d[30] * d[33] * d[41] +
        2 * d[31] * d[34] * d[41] + 2 * d[32] * d[35] * d[41] + 2 * d[27] * d[29] * d[42] + 2 * d[30] * d[32] * d[42] +
        2 * d[33] * d[35] * d[42] + 2 * d[28] * d[29] * d[43] + 2 * d[31] * d[32] * d[43] + 2 * d[34] * d[35] * d[43] -
        std::pow(d[27], 2) * d[44] - std::pow(d[28], 2) * d[44] + std::pow(d[29], 2) * d[44] -
        std::pow(d[30], 2) * d[44] - std::pow(d[31], 2) * d[44] + std::pow(d[32], 2) * d[44] +
        std::pow(d[33], 2) * d[44] + std::pow(d[34], 2) * d[44] + 3 * std::pow(d[35], 2) * d[44];
    coeffs[380] =
        -d[8] * std::pow(d[36], 2) - d[8] * std::pow(d[37], 2) + 2 * d[6] * d[36] * d[38] + 2 * d[7] * d[37] * d[38] +
        d[8] * std::pow(d[38], 2) - d[8] * std::pow(d[39], 2) - d[8] * std::pow(d[40], 2) + 2 * d[6] * d[39] * d[41] +
        2 * d[7] * d[40] * d[41] + d[8] * std::pow(d[41], 2) + 2 * d[2] * d[36] * d[42] + 2 * d[0] * d[38] * d[42] +
        2 * d[5] * d[39] * d[42] + 2 * d[3] * d[41] * d[42] + d[8] * std::pow(d[42], 2) + 2 * d[2] * d[37] * d[43] +
        2 * d[1] * d[38] * d[43] + 2 * d[5] * d[40] * d[43] + 2 * d[4] * d[41] * d[43] + d[8] * std::pow(d[43], 2) -
        2 * d[0] * d[36] * d[44] - 2 * d[1] * d[37] * d[44] + 2 * d[2] * d[38] * d[44] - 2 * d[3] * d[39] * d[44] -
        2 * d[4] * d[40] * d[44] + 2 * d[5] * d[41] * d[44] + 2 * d[6] * d[42] * d[44] + 2 * d[7] * d[43] * d[44] +
        3 * d[8] * std::pow(d[44], 2);
    coeffs[381] = -d[17] * std::pow(d[36], 2) - d[17] * std::pow(d[37], 2) + 2 * d[15] * d[36] * d[38] +
                  2 * d[16] * d[37] * d[38] + d[17] * std::pow(d[38], 2) - d[17] * std::pow(d[39], 2) -
                  d[17] * std::pow(d[40], 2) + 2 * d[15] * d[39] * d[41] + 2 * d[16] * d[40] * d[41] +
                  d[17] * std::pow(d[41], 2) + 2 * d[11] * d[36] * d[42] + 2 * d[9] * d[38] * d[42] +
                  2 * d[14] * d[39] * d[42] + 2 * d[12] * d[41] * d[42] + d[17] * std::pow(d[42], 2) +
                  2 * d[11] * d[37] * d[43] + 2 * d[10] * d[38] * d[43] + 2 * d[14] * d[40] * d[43] +
                  2 * d[13] * d[41] * d[43] + d[17] * std::pow(d[43], 2) - 2 * d[9] * d[36] * d[44] -
                  2 * d[10] * d[37] * d[44] + 2 * d[11] * d[38] * d[44] - 2 * d[12] * d[39] * d[44] -
                  2 * d[13] * d[40] * d[44] + 2 * d[14] * d[41] * d[44] + 2 * d[15] * d[42] * d[44] +
                  2 * d[16] * d[43] * d[44] + 3 * d[17] * std::pow(d[44], 2);
    coeffs[382] = -d[26] * std::pow(d[36], 2) - d[26] * std::pow(d[37], 2) + 2 * d[24] * d[36] * d[38] +
                  2 * d[25] * d[37] * d[38] + d[26] * std::pow(d[38], 2) - d[26] * std::pow(d[39], 2) -
                  d[26] * std::pow(d[40], 2) + 2 * d[24] * d[39] * d[41] + 2 * d[25] * d[40] * d[41] +
                  d[26] * std::pow(d[41], 2) + 2 * d[20] * d[36] * d[42] + 2 * d[18] * d[38] * d[42] +
                  2 * d[23] * d[39] * d[42] + 2 * d[21] * d[41] * d[42] + d[26] * std::pow(d[42], 2) +
                  2 * d[20] * d[37] * d[43] + 2 * d[19] * d[38] * d[43] + 2 * d[23] * d[40] * d[43] +
                  2 * d[22] * d[41] * d[43] + d[26] * std::pow(d[43], 2) - 2 * d[18] * d[36] * d[44] -
                  2 * d[19] * d[37] * d[44] + 2 * d[20] * d[38] * d[44] - 2 * d[21] * d[39] * d[44] -
                  2 * d[22] * d[40] * d[44] + 2 * d[23] * d[41] * d[44] + 2 * d[24] * d[42] * d[44] +
                  2 * d[25] * d[43] * d[44] + 3 * d[26] * std::pow(d[44], 2);
    coeffs[383] = -d[35] * std::pow(d[36], 2) - d[35] * std::pow(d[37], 2) + 2 * d[33] * d[36] * d[38] +
                  2 * d[34] * d[37] * d[38] + d[35] * std::pow(d[38], 2) - d[35] * std::pow(d[39], 2) -
                  d[35] * std::pow(d[40], 2) + 2 * d[33] * d[39] * d[41] + 2 * d[34] * d[40] * d[41] +
                  d[35] * std::pow(d[41], 2) + 2 * d[29] * d[36] * d[42] + 2 * d[27] * d[38] * d[42] +
                  2 * d[32] * d[39] * d[42] + 2 * d[30] * d[41] * d[42] + d[35] * std::pow(d[42], 2) +
                  2 * d[29] * d[37] * d[43] + 2 * d[28] * d[38] * d[43] + 2 * d[32] * d[40] * d[43] +
                  2 * d[31] * d[41] * d[43] + d[35] * std::pow(d[43], 2) - 2 * d[27] * d[36] * d[44] -
                  2 * d[28] * d[37] * d[44] + 2 * d[29] * d[38] * d[44] - 2 * d[30] * d[39] * d[44] -
                  2 * d[31] * d[40] * d[44] + 2 * d[32] * d[41] * d[44] + 2 * d[33] * d[42] * d[44] +
                  2 * d[34] * d[43] * d[44] + 3 * d[35] * std::pow(d[44], 2);
    coeffs[384] = 2 * d[36] * d[38] * d[42] + 2 * d[39] * d[41] * d[42] + 2 * d[37] * d[38] * d[43] +
                  2 * d[40] * d[41] * d[43] - std::pow(d[36], 2) * d[44] - std::pow(d[37], 2) * d[44] +
                  std::pow(d[38], 2) * d[44] - std::pow(d[39], 2) * d[44] - std::pow(d[40], 2) * d[44] +
                  std::pow(d[41], 2) * d[44] + std::pow(d[42], 2) * d[44] + std::pow(d[43], 2) * d[44] +
                  std::pow(d[44], 3);

    // Setup elimination template
    static const int coeffs0_ind[] = {
        0,   35,  1,   0,   35,  36,  2,   1,   0,   35,  36,  37,  70,  105, 140, 3,   2,   1,   36,  37,  38,  71,
        106, 141, 3,   2,   37,  38,  72,  107, 142, 3,   38,  73,  108, 143, 4,   39,  0,   35,  5,   4,   39,  40,
        1,   0,   35,  36,  70,  105, 140, 6,   5,   4,   39,  40,  41,  2,   1,   36,  37,  71,  74,  106, 109, 141,
        144, 6,   5,   40,  41,  3,   2,   37,  38,  72,  75,  107, 110, 142, 145, 6,   41,  3,   38,  73,  76,  108,
        111, 143, 146, 7,   42,  4,   39,  35,  0,   70,  105, 140, 175, 210, 245, 280, 315, 350, 8,   7,   42,  43,
        5,   4,   39,  40,  74,  36,  1,   71,  106, 109, 141, 144, 176, 211, 246, 281, 316, 351, 8,   7,   42,  43,
        6,   5,   40,  41,  75,  37,  2,   72,  77,  107, 110, 112, 142, 145, 147, 177, 212, 247, 282, 317, 352, 8,
        43,  6,   41,  76,  38,  3,   73,  78,  108, 111, 113, 143, 146, 148, 178, 213, 248, 283, 318, 353, 9,   44,
        7,   42,  39,  4,   74,  109, 144, 179, 214, 249, 284, 319, 354, 9,   44,  8,   7,   42,  43,  77,  40,  5,
        75,  110, 112, 145, 147, 180, 215, 250, 285, 320, 355, 9,   44,  8,   43,  78,  41,  6,   76,  79,  111, 113,
        114, 146, 148, 149, 181, 216, 251, 286, 321, 356, 9,   44,  42,  7,   77,  112, 147, 182, 217, 252, 287, 322,
        357, 9,   44,  79,  43,  8,   78,  113, 114, 148, 149, 183, 218, 253, 288, 323, 358, 44,  9,   79,  114, 149,
        184, 219, 254, 289, 324, 359, 10,  45,  35,  0,   11,  10,  45,  46,  36,  35,  0,   70,  105, 140, 1,   12,
        11,  10,  45,  46,  47,  80,  115, 150, 37,  36,  1,   71,  106, 141, 2,   12,  11,  46,  47,  81,  116, 151,
        38,  37,  2,   72,  107, 142, 3,   12,  47,  82,  117, 152, 38,  3,   73,  108, 143, 13,  48,  10,  45,  39,
        0,   70,  105, 140, 35,  4,   175, 210, 245, 280, 315, 350, 14,  13,  48,  49,  11,  10,  45,  46,  80,  115,
        150, 40,  39,  4,   74,  109, 1,   71,  106, 141, 144, 36,  5,   176, 211, 246, 281, 316, 351, 14,  13,  48,
        49,  12,  11,  46,  47,  81,  83,  116, 118, 151, 153, 41,  40,  5,   75,  110, 2,   72,  107, 142, 145, 37,
        6,   177, 212, 247, 282, 317, 352, 14,  49,  12,  47,  82,  84,  117, 119, 152, 154, 41,  6,   76,  111, 3,
        73,  108, 143, 146, 38,  178, 213, 248, 283, 318, 353, 15,  50,  13,  48,  45,  10,  80,  115, 150, 42,  4,
        74,  109, 144, 39,  7,   179, 185, 214, 220, 249, 255, 284, 290, 319, 325, 354, 360, 15,  50,  14,  13,  48,
        49,  83,  46,  11,  81,  116, 118, 151, 153, 43,  42,  7,   77,  112, 5,   75,  110, 145, 147, 40,  8,   180,
        186, 215, 221, 250, 256, 285, 291, 320, 326, 355, 361, 15,  50,  14,  49,  84,  47,  12,  82,  85,  117, 119,
        120, 152, 154, 155, 43,  8,   78,  113, 6,   76,  111, 146, 148, 41,  181, 187, 216, 222, 251, 257, 286, 292,
        321, 327, 356, 362, 15,  50,  48,  13,  83,  118, 153, 44,  7,   77,  112, 147, 42,  9,   182, 188, 217, 223,
        252, 258, 287, 293, 322, 328, 357, 363, 15,  50,  85,  49,  14,  84,  119, 120, 154, 155, 44,  9,   79,  114,
        8,   78,  113, 148, 149, 43,  183, 189, 218, 224, 253, 259, 288, 294, 323, 329, 358, 364, 50,  15,  85,  120,
        155, 9,   79,  114, 149, 44,  184, 190, 219, 225, 254, 260, 289, 295, 324, 330, 359, 365, 16,  51,  45,  10,
        105, 35,  140, 70,  175, 0,   210, 245, 280, 315, 350, 17,  16,  51,  52,  46,  45,  10,  80,  115, 150, 11,
        106, 36,  141, 71,  176, 1,   211, 246, 281, 316, 351, 17,  16,  51,  52,  86,  121, 156, 47,  46,  11,  81,
        116, 151, 12,  107, 37,  142, 72,  177, 2,   212, 247, 282, 317, 352, 17,  52,  87,  122, 157, 47,  12,  82,
        117, 152, 108, 38,  143, 73,  178, 3,   213, 248, 283, 318, 353, 18,  53,  16,  51,  48,  10,  80,  115, 150,
        45,  13,  185, 220, 109, 39,  144, 74,  179, 4,   214, 249, 255, 284, 290, 319, 325, 354, 360, 18,  53,  17,
        16,  51,  52,  86,  121, 156, 49,  48,  13,  83,  118, 11,  81,  116, 151, 153, 46,  14,  186, 221, 110, 40,
        145, 75,  180, 5,   215, 250, 256, 285, 291, 320, 326, 355, 361, 18,  53,  17,  52,  87,  88,  122, 123, 157,
        158, 49,  14,  84,  119, 12,  82,  117, 152, 154, 47,  187, 222, 111, 41,  146, 76,  181, 6,   216, 251, 257,
        286, 292, 321, 327, 356, 362, 18,  53,  51,  16,  86,  121, 156, 50,  13,  83,  118, 153, 48,  15,  188, 191,
        223, 112, 42,  147, 77,  182, 7,   217, 226, 252, 258, 261, 287, 293, 296, 322, 328, 331, 357, 363, 366, 18,
        53,  88,  52,  17,  87,  122, 123, 157, 158, 50,  15,  85,  120, 14,  84,  119, 154, 155, 49,  189, 192, 224,
        113, 43,  148, 78,  183, 8,   218, 227, 253, 259, 262, 288, 294, 297, 323, 329, 332, 358, 364, 367, 53,  18,
        88,  123, 158, 15,  85,  120, 155, 50,  190, 193, 225, 114, 44,  149, 79,  184, 9,   219, 228, 254, 260, 263,
        289, 295, 298, 324, 330, 333, 359, 365, 368, 19,  54,  51,  16,  115, 45,  150, 80,  185, 10,  220, 255, 290,
        325, 360, 19,  54,  52,  51,  16,  86,  121, 156, 17,  116, 46,  151, 81,  186, 11,  221, 256, 291, 326, 361,
        19,  54,  89,  124, 159, 52,  17,  87,  122, 157, 117, 47,  152, 82,  187, 12,  222, 257, 292, 327, 362, 19,
        54,  53,  16,  86,  121, 156, 51,  18,  191, 226, 118, 48,  153, 83,  188, 13,  223, 258, 261, 293, 296, 328,
        331, 363, 366, 20,  55,  0,   35,  21,  20,  55,  56,  1,   0,   140, 70,  35,  105, 36,  22,  21,  20,  55,
        56,  57,  90,  125, 160, 2,   1,   141, 71,  36,  106, 37,  22,  21,  56,  57,  91,  126, 161, 3,   2,   142,
        72,  37,  107, 38,  22,  57,  92,  127, 162, 3,   143, 73,  38,  108, 23,  58,  20,  55,  4,   35,  175, 70,
        140, 210, 105, 245, 39,  280, 0,   315, 350, 24,  23,  58,  59,  21,  20,  55,  56,  90,  125, 160, 5,   4,
        144, 74,  39,  36,  176, 71,  141, 211, 106, 246, 109, 40,  281, 1,   316, 351, 24,  23,  58,  59,  22,  21,
        56,  57,  91,  93,  126, 128, 161, 163, 6,   5,   145, 75,  40,  37,  177, 72,  142, 212, 107, 247, 110, 41,
        282, 2,   317, 352, 24,  59,  22,  57,  92,  94,  127, 129, 162, 164, 6,   146, 76,  41,  38,  178, 73,  143,
        213, 108, 248, 111, 283, 3,   318, 353, 25,  60,  23,  58,  55,  20,  90,  125, 160, 195, 230, 265, 300, 7,
        39,  179, 74,  144, 214, 109, 249, 42,  284, 4,   319, 335, 354, 370, 25,  60,  24,  23,  58,  59,  93,  56,
        21,  91,  126, 128, 161, 163, 196, 231, 266, 301, 8,   7,   147, 77,  42,  40,  180, 75,  145, 215, 110, 250,
        112, 43,  285, 5,   320, 336, 355, 371, 25,  60,  24,  59,  94,  57,  22,  92,  95,  127, 129, 130, 162, 164,
        165, 197, 232, 267, 302, 8,   148, 78,  43,  41,  181, 76,  146, 216, 111, 251, 113, 286, 6,   321, 337, 356,
        372, 25,  60,  58,  23,  93,  128, 163, 198, 233, 268, 303, 9,   42,  182, 77,  147, 217, 112, 252, 44,  287,
        7,   322, 338, 357, 373, 25,  60,  95,  59,  24,  94,  129, 130, 164, 165, 199, 234, 269, 304, 9,   149, 79,
        44,  43,  183, 78,  148, 218, 113, 253, 114, 288, 8,   323, 339, 358, 374, 60,  25,  95,  130, 165, 200, 235,
        270, 305, 44,  184, 79,  149, 219, 114, 254, 289, 9,   324, 340, 359, 375, 26,  61,  55,  20,  10,  45,  105,
        210, 70,  0,   245, 140, 175, 280, 35,  315, 350, 27,  26,  61,  62,  56,  55,  20,  90,  125, 160, 21,  11,
        10,  150, 80,  45,  115, 46,  106, 211, 71,  1,   246, 141, 176, 281, 36,  316, 351, 27,  26,  61,  62,  96,
        131, 166, 57,  56,  21,  91,  126, 161, 22,  12,  11,  151, 81,  46,  116, 47,  107, 212, 72,  2,   247, 142,
        177, 282, 37,  317, 352, 27,  62,  97,  132, 167, 57,  22,  92,  127, 162, 12,  152, 82,  47,  117, 108, 213,
        73,  3,   248, 143, 178, 283, 38,  318, 353, 28,  63,  26,  61,  58,  20,  90,  125, 160, 55,  23,  195, 230,
        265, 300, 13,  45,  185, 80,  150, 220, 115, 255, 48,  290, 109, 214, 74,  4,   249, 144, 179, 284, 39,  10,
        319, 325, 335, 354, 360, 370, 28,  63,  27,  26,  61,  62,  96,  131, 166, 59,  58,  23,  93,  128, 21,  91,
        126, 161, 163, 56,  24,  196, 231, 266, 301, 14,  13,  153, 83,  48,  46,  186, 81,  151, 221, 116, 256, 118,
        49,  291, 110, 215, 75,  5,   250, 145, 180, 285, 40,  11,  320, 326, 336, 355, 361, 371, 28,  63,  27,  62,
        97,  98,  132, 133, 167, 168, 59,  24,  94,  129, 22,  92,  127, 162, 164, 57,  197, 232, 267, 302, 14,  154,
        84,  49,  47,  187, 82,  152, 222, 117, 257, 119, 292, 111, 216, 76,  6,   251, 146, 181, 286, 41,  12,  321,
        327, 337, 356, 362, 372, 28,  63,  61,  26,  96,  131, 166, 60,  23,  93,  128, 163, 58,  25,  198, 201, 233,
        236, 268, 271, 303, 306, 15,  48,  188, 83,  153, 223, 118, 258, 50,  293, 112, 217, 77,  7,   252, 147, 182,
        287, 42,  13,  322, 328, 338, 341, 357, 363, 373, 376, 28,  63,  98,  62,  27,  97,  132, 133, 167, 168, 60,
        25,  95,  130, 24,  94,  129, 164, 165, 59,  199, 202, 234, 237, 269, 272, 304, 307, 15,  155, 85,  50,  49,
        189, 84,  154, 224, 119, 259, 120, 294, 113, 218, 78,  8,   253, 148, 183, 288, 43,  14,  323, 329, 339, 342,
        358, 364, 374, 377, 63,  28,  98,  133, 168, 25,  95,  130, 165, 60,  200, 203, 235, 238, 270, 273, 305, 308,
        50,  190, 85,  155, 225, 120, 260, 295, 114, 219, 79,  9,   254, 149, 184, 289, 44,  15,  324, 330, 340, 343,
        359, 365, 375, 378, 29,  64,  61,  26,  125, 55,  160, 90,  195, 20,  230, 265, 300, 16,  51,  115, 220, 80,
        10,  255, 150, 185, 290, 45,  325, 335, 360, 370, 30,  65,  20,  55,  35,  105, 280, 245, 140, 70,  0,   210,
        315, 175, 350, 31,  30,  65,  66,  21,  20,  160, 90,  55,  125, 56,  36,  106, 281, 246, 141, 71,  1,   211,
        316, 176, 351, 31,  30,  65,  66,  100, 135, 170, 22,  21,  161, 91,  56,  126, 57,  37,  107, 282, 247, 142,
        72,  2,   212, 317, 177, 352, 31,  66,  101, 136, 171, 22,  162, 92,  57,  127, 38,  108, 283, 248, 143, 73,
        3,   213, 318, 178, 353, 32,  67,  30,  65,  23,  55,  195, 90,  160, 230, 125, 265, 58,  300, 20,  39,  109,
        284, 249, 144, 74,  4,   214, 319, 179, 335, 354, 370, 32,  67,  31,  30,  65,  66,  100, 135, 170, 24,  23,
        163, 93,  58,  56,  196, 91,  161, 231, 126, 266, 128, 59,  301, 21,  40,  110, 285, 250, 145, 75,  5,   215,
        320, 180, 336, 355, 371, 32,  67,  31,  66,  101, 102, 136, 137, 171, 172, 24,  164, 94,  59,  57,  197, 92,
        162, 232, 127, 267, 129, 302, 22,  41,  111, 286, 251, 146, 76,  6,   216, 321, 181, 337, 356, 372, 32,  67,
        65,  30,  100, 135, 170, 205, 240, 275, 310, 25,  58,  198, 93,  163, 233, 128, 268, 60,  303, 23,  42,  112,
        287, 252, 147, 77,  7,   217, 322, 182, 338, 345, 357, 373, 380, 32,  67,  102, 66,  31,  101, 136, 137, 171,
        172, 206, 241, 276, 311, 25,  165, 95,  60,  59,  199, 94,  164, 234, 129, 269, 130, 304, 24,  43,  113, 288,
        253, 148, 78,  8,   218, 323, 183, 339, 346, 358, 374, 381, 67,  32,  102, 137, 172, 207, 242, 277, 312, 60,
        200, 95,  165, 235, 130, 270, 305, 25,  44,  114, 289, 254, 149, 79,  9,   219, 324, 184, 340, 347, 359, 375,
        382, 33,  68,  65,  30,  26,  61,  125, 230, 90,  20,  265, 160, 195, 300, 55,  335, 45,  115, 290, 255, 150,
        80,  10,  220, 325, 185, 360, 370, 29,  64,  62,  61,  26,  96,  131, 166, 27,  126, 56,  161, 91,  196, 21,
        231, 266, 301, 17,  16,  156, 86,  51,  121, 52,  116, 221, 81,  11,  256, 151, 186, 291, 46,  326, 336, 361,
        371, 29,  64,  63,  26,  96,  131, 166, 61,  28,  201, 236, 128, 58,  163, 93,  198, 23,  233, 268, 271, 303,
        306, 18,  51,  191, 86,  156, 226, 121, 261, 53,  296, 118, 223, 83,  13,  258, 153, 188, 293, 48,  16,  328,
        338, 331, 341, 363, 373, 366, 376, 54,  19,  121, 51,  156, 86,  191, 16,  226, 261, 296, 331, 366, 29,  64,
        99,  134, 169, 62,  27,  97,  132, 167, 127, 57,  162, 92,  197, 22,  232, 267, 302, 17,  157, 87,  52,  122,
        117, 222, 82,  12,  257, 152, 187, 292, 47,  327, 337, 362, 372, 19,  54,  89,  124, 159, 53,  18,  88,  123,
        17,  87,  122, 157, 158, 52,  192, 227, 119, 49,  154, 84,  189, 14,  224, 259, 262, 294, 297, 329, 332, 364,
        367, 54,  19,  89,  124, 159, 122, 52,  157, 87,  192, 17,  227, 262, 297, 332, 367, 54,  19,  89,  124, 159,
        18,  88,  123, 158, 53,  193, 194, 228, 120, 50,  155, 85,  190, 15,  225, 229, 260, 263, 264, 295, 298, 299,
        330, 333, 334, 365, 368, 369, 19,  89,  124, 159, 54,  194, 229, 123, 53,  158, 88,  193, 18,  228, 263, 264,
        298, 299, 333, 334, 368, 369, 124, 54,  159, 89,  194, 19,  229, 264, 299, 334, 369};
    static const int coeffs1_ind[] = {
        69,  139, 314, 279, 174, 104, 34,  244, 349, 209, 384, 34,  69,  65,  135, 310, 275, 170, 100, 30,  240, 345,
        205, 380, 34,  69,  30,  65,  55,  125, 300, 265, 160, 90,  20,  230, 335, 195, 370, 34,  69,  31,  30,  170,
        100, 65,  135, 66,  56,  126, 301, 266, 161, 91,  21,  231, 336, 196, 371, 33,  68,  66,  65,  30,  100, 135,
        170, 31,  27,  26,  166, 96,  61,  131, 62,  126, 231, 91,  21,  266, 161, 196, 301, 56,  336, 46,  116, 291,
        256, 151, 81,  11,  221, 326, 186, 361, 371, 34,  69,  32,  65,  205, 100, 170, 240, 135, 275, 67,  310, 30,
        58,  128, 303, 268, 163, 93,  23,  233, 338, 198, 345, 373, 380, 33,  68,  67,  30,  100, 135, 170, 65,  32,
        205, 240, 275, 310, 28,  61,  201, 96,  166, 236, 131, 271, 63,  306, 128, 233, 93,  23,  268, 163, 198, 303,
        58,  26,  338, 48,  118, 293, 258, 153, 83,  13,  223, 328, 188, 341, 345, 363, 373, 376, 380, 69,  34,  33,
        68,  135, 240, 100, 30,  275, 170, 205, 310, 65,  345, 61,  131, 306, 271, 166, 96,  26,  236, 341, 201, 376,
        380, 68,  33,  135, 65,  170, 100, 205, 30,  240, 275, 310, 29,  64,  131, 236, 96,  26,  271, 166, 201, 306,
        61,  341, 51,  121, 296, 261, 156, 86,  16,  226, 331, 191, 345, 366, 376, 380, 64,  29,  131, 61,  166, 96,
        201, 26,  236, 271, 306, 19,  54,  121, 226, 86,  16,  261, 156, 191, 296, 51,  331, 341, 366, 376, 34,  174,
        104, 69,  139, 66,  136, 311, 276, 171, 101, 31,  241, 346, 206, 381, 34,  69,  104, 139, 174, 31,  171, 101,
        66,  136, 57,  127, 302, 267, 162, 92,  22,  232, 337, 197, 372, 33,  68,  103, 138, 173, 66,  31,  101, 136,
        171, 27,  167, 97,  62,  132, 127, 232, 92,  22,  267, 162, 197, 302, 57,  337, 47,  117, 292, 257, 152, 82,
        12,  222, 327, 187, 362, 372, 34,  69,  104, 139, 174, 32,  172, 102, 67,  66,  206, 101, 171, 241, 136, 276,
        137, 311, 31,  59,  129, 304, 269, 164, 94,  24,  234, 339, 199, 346, 374, 381, 33,  68,  103, 138, 173, 67,
        32,  102, 137, 31,  101, 136, 171, 172, 66,  206, 241, 276, 311, 28,  168, 98,  63,  62,  202, 97,  167, 237,
        132, 272, 133, 307, 129, 234, 94,  24,  269, 164, 199, 304, 59,  27,  339, 49,  119, 294, 259, 154, 84,  14,
        224, 329, 189, 342, 346, 364, 374, 377, 381, 29,  64,  99,  134, 169, 63,  28,  98,  133, 27,  97,  132, 167,
        168, 62,  202, 237, 129, 59,  164, 94,  199, 24,  234, 269, 272, 304, 307, 18,  158, 88,  53,  52,  192, 87,
        157, 227, 122, 262, 123, 297, 119, 224, 84,  14,  259, 154, 189, 294, 49,  17,  329, 339, 332, 342, 364, 374,
        367, 377, 69,  34,  104, 139, 174, 33,  173, 103, 68,  138, 136, 241, 101, 31,  276, 171, 206, 311, 66,  346,
        62,  132, 307, 272, 167, 97,  27,  237, 342, 202, 377, 381, 68,  33,  103, 138, 173, 136, 66,  171, 101, 206,
        31,  241, 276, 311, 29,  169, 99,  64,  134, 132, 237, 97,  27,  272, 167, 202, 307, 62,  342, 52,  122, 297,
        262, 157, 87,  17,  227, 332, 192, 346, 367, 377, 381, 64,  29,  99,  134, 169, 132, 62,  167, 97,  202, 27,
        237, 272, 307, 19,  159, 89,  54,  124, 122, 227, 87,  17,  262, 157, 192, 297, 52,  332, 342, 367, 377, 69,
        209, 104, 174, 244, 139, 279, 314, 34,  67,  137, 312, 277, 172, 102, 32,  242, 347, 207, 349, 382, 384, 69,
        34,  104, 139, 174, 209, 244, 279, 314, 67,  207, 102, 172, 242, 137, 277, 312, 32,  60,  130, 305, 270, 165,
        95,  25,  235, 340, 200, 347, 349, 375, 382, 384, 68,  33,  103, 138, 173, 32,  102, 137, 172, 67,  207, 208,
        242, 243, 277, 278, 312, 313, 63,  203, 98,  168, 238, 133, 273, 308, 130, 235, 95,  25,  270, 165, 200, 305,
        60,  28,  340, 50,  120, 295, 260, 155, 85,  15,  225, 330, 190, 343, 347, 348, 365, 375, 378, 382, 383, 64,
        29,  99,  134, 169, 28,  98,  133, 168, 63,  203, 204, 238, 130, 60,  165, 95,  200, 25,  235, 239, 270, 273,
        274, 305, 308, 309, 53,  193, 88,  158, 228, 123, 263, 298, 120, 225, 85,  15,  260, 155, 190, 295, 50,  18,
        330, 340, 333, 343, 344, 365, 375, 368, 378, 379, 34,  104, 139, 174, 69,  209, 244, 279, 314, 68,  208, 103,
        173, 243, 138, 278, 313, 137, 242, 102, 32,  277, 172, 207, 312, 67,  33,  347, 63,  133, 308, 273, 168, 98,
        28,  238, 343, 203, 348, 349, 378, 382, 383, 384, 33,  103, 138, 173, 68,  208, 243, 137, 67,  172, 102, 207,
        32,  242, 277, 278, 312, 313, 64,  204, 99,  169, 239, 134, 274, 309, 133, 238, 98,  28,  273, 168, 203, 308,
        63,  29,  343, 53,  123, 298, 263, 158, 88,  18,  228, 333, 193, 347, 344, 348, 368, 378, 382, 379, 383, 29,
        99,  134, 169, 64,  204, 239, 133, 63,  168, 98,  203, 28,  238, 273, 274, 308, 309, 54,  194, 89,  159, 229,
        124, 264, 299, 123, 228, 88,  18,  263, 158, 193, 298, 53,  19,  333, 343, 334, 344, 368, 378, 369, 379, 139,
        244, 104, 34,  279, 174, 209, 314, 69,  349, 68,  138, 313, 278, 173, 103, 33,  243, 348, 208, 383, 384, 139,
        69,  174, 104, 209, 34,  244, 279, 314, 138, 243, 103, 33,  278, 173, 208, 313, 68,  348, 64,  134, 309, 274,
        169, 99,  29,  239, 344, 204, 349, 379, 383, 384, 138, 68,  173, 103, 208, 33,  243, 278, 313, 134, 239, 99,
        29,  274, 169, 204, 309, 64,  344, 54,  124, 299, 264, 159, 89,  19,  229, 334, 194, 348, 369, 379, 383, 134,
        64,  169, 99,  204, 29,  239, 274, 309, 124, 229, 89,  19,  264, 159, 194, 299, 54,  334, 344, 369, 379};

    static const int C0_ind[] = {
        0,    5,    96,   97,   100,  101,  192,  193,  194,  195,  196,  197,  206,  209,  212,  288,  289,  290,
        291,  292,  293,  302,  305,  308,  385,  386,  387,  388,  398,  401,  404,  482,  483,  494,  497,  500,
        576,  581,  582,  585,  672,  673,  676,  677,  678,  679,  680,  681,  682,  688,  691,  768,  769,  770,
        771,  772,  773,  774,  775,  776,  777,  778,  782,  784,  785,  787,  788,  865,  866,  867,  868,  870,
        871,  872,  873,  874,  878,  880,  881,  883,  884,  962,  963,  967,  968,  970,  974,  976,  977,  979,
        980,  1056, 1061, 1062, 1065, 1067, 1068, 1069, 1071, 1074, 1090, 1099, 1102, 1105, 1145, 1151, 1152, 1153,
        1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1167, 1168, 1170, 1171, 1186, 1195, 1198, 1201,
        1241, 1247, 1249, 1250, 1251, 1252, 1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1262, 1263, 1264, 1265,
        1266, 1267, 1268, 1282, 1291, 1294, 1297, 1337, 1343, 1346, 1347, 1351, 1352, 1354, 1355, 1356, 1357, 1358,
        1359, 1360, 1361, 1362, 1363, 1364, 1378, 1387, 1390, 1393, 1433, 1439, 1440, 1445, 1446, 1449, 1451, 1452,
        1453, 1455, 1458, 1474, 1483, 1486, 1489, 1529, 1535, 1537, 1540, 1542, 1543, 1544, 1545, 1546, 1547, 1548,
        1549, 1551, 1552, 1554, 1555, 1570, 1579, 1582, 1585, 1625, 1631, 1634, 1635, 1639, 1640, 1642, 1643, 1644,
        1645, 1646, 1647, 1648, 1649, 1650, 1651, 1652, 1666, 1675, 1678, 1681, 1721, 1727, 1734, 1737, 1739, 1740,
        1741, 1743, 1746, 1762, 1771, 1774, 1777, 1817, 1823, 1831, 1832, 1834, 1835, 1836, 1837, 1839, 1840, 1842,
        1843, 1858, 1867, 1870, 1873, 1913, 1919, 1931, 1932, 1933, 1935, 1938, 1954, 1963, 1966, 1969, 2009, 2015,
        2016, 2021, 2037, 2048, 2112, 2113, 2116, 2117, 2133, 2134, 2135, 2136, 2137, 2142, 2144, 2208, 2209, 2210,
        2211, 2212, 2213, 2222, 2225, 2228, 2229, 2230, 2231, 2232, 2233, 2238, 2240, 2305, 2306, 2307, 2308, 2318,
        2321, 2324, 2325, 2326, 2327, 2328, 2329, 2334, 2336, 2402, 2403, 2414, 2417, 2420, 2422, 2423, 2424, 2425,
        2430, 2496, 2501, 2502, 2505, 2517, 2522, 2523, 2524, 2525, 2527, 2528, 2529, 2531, 2541, 2544, 2584, 2590,
        2592, 2593, 2596, 2597, 2598, 2599, 2600, 2601, 2602, 2608, 2611, 2613, 2614, 2615, 2616, 2617, 2618, 2619,
        2620, 2621, 2622, 2623, 2624, 2625, 2627, 2637, 2640, 2680, 2686, 2689, 2690, 2691, 2692, 2694, 2695, 2696,
        2697, 2698, 2702, 2704, 2705, 2707, 2708, 2709, 2710, 2711, 2712, 2713, 2714, 2715, 2716, 2717, 2718, 2719,
        2720, 2721, 2723, 2733, 2736, 2776, 2782, 2786, 2787, 2791, 2792, 2794, 2798, 2800, 2801, 2803, 2804, 2806,
        2807, 2808, 2809, 2810, 2811, 2812, 2813, 2814, 2815, 2817, 2819, 2829, 2832, 2872, 2878, 2880, 2885, 2886,
        2889, 2891, 2892, 2893, 2895, 2898, 2901, 2906, 2907, 2908, 2909, 2911, 2912, 2913, 2914, 2915, 2923, 2925,
        2926, 2928, 2929, 2968, 2969, 2974, 2975, 2977, 2980, 2982, 2983, 2984, 2985, 2986, 2987, 2988, 2989, 2991,
        2992, 2994, 2995, 2997, 2998, 2999, 3000, 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3008, 3009, 3010, 3011,
        3019, 3021, 3022, 3024, 3025, 3064, 3065, 3070, 3071, 3074, 3075, 3079, 3080, 3082, 3083, 3084, 3085, 3086,
        3087, 3088, 3089, 3090, 3091, 3092, 3094, 3095, 3096, 3097, 3098, 3099, 3100, 3101, 3102, 3103, 3105, 3106,
        3107, 3115, 3117, 3118, 3120, 3121, 3160, 3161, 3166, 3167, 3174, 3177, 3179, 3180, 3181, 3183, 3186, 3189,
        3194, 3195, 3196, 3197, 3199, 3200, 3201, 3202, 3203, 3211, 3213, 3214, 3216, 3217, 3256, 3257, 3262, 3263,
        3271, 3272, 3274, 3275, 3276, 3277, 3279, 3280, 3282, 3283, 3286, 3287, 3288, 3289, 3290, 3291, 3292, 3293,
        3294, 3295, 3297, 3298, 3299, 3307, 3309, 3310, 3312, 3313, 3352, 3353, 3358, 3359, 3371, 3372, 3373, 3375,
        3378, 3386, 3387, 3388, 3389, 3391, 3393, 3394, 3395, 3403, 3405, 3406, 3408, 3409, 3448, 3449, 3454, 3455,
        3456, 3461, 3477, 3488, 3492, 3493, 3494, 3495, 3496, 3497, 3498, 3500, 3503, 3542, 3548, 3552, 3553, 3556,
        3557, 3573, 3574, 3575, 3576, 3577, 3582, 3584, 3588, 3589, 3590, 3591, 3592, 3593, 3594, 3596, 3599, 3638,
        3644, 3649, 3650, 3651, 3652, 3662, 3665, 3668, 3669, 3670, 3671, 3672, 3673, 3678, 3680, 3684, 3685, 3686,
        3687, 3688, 3689, 3690, 3692, 3695, 3734, 3740, 3746, 3747, 3758, 3761, 3764, 3766, 3767, 3768, 3769, 3774,
        3780, 3781, 3782, 3783, 3784, 3785, 3786, 3788, 3791, 3830, 3836, 3840, 3845, 3846, 3849, 3861, 3866, 3867,
        3868, 3869, 3871, 3872, 3873, 3875, 3876, 3877, 3878, 3879, 3880, 3881, 3882, 3884, 3885, 3887, 3888, 3926,
        3928, 3932, 3934, 3937, 3940, 3942, 3943, 3944, 3945, 3946, 3952, 3955, 3957, 3958, 3959, 3960, 3961, 3962,
        3963, 3964, 3965, 3966, 3967, 3968, 3969, 3971, 3972, 3973, 3974, 3975, 3976, 3977, 3978, 3980, 3981, 3983,
        3984, 4022, 4024, 4028, 4030, 4034, 4035, 4039, 4040, 4042, 4046, 4048, 4049, 4051, 4052, 4054, 4055, 4056,
        4057, 4058, 4059, 4060, 4061, 4062, 4063, 4065, 4067, 4068, 4069, 4070, 4071, 4072, 4073, 4074, 4076, 4077,
        4079, 4080, 4118, 4120, 4124, 4126, 4134, 4137, 4139, 4140, 4141, 4143, 4146, 4149, 4154, 4155, 4156, 4157,
        4159, 4160, 4161, 4162, 4163, 4164, 4165, 4166, 4167, 4168, 4169, 4170, 4171, 4172, 4173, 4174, 4175, 4176,
        4177, 4214, 4216, 4217, 4220, 4222, 4223, 4231, 4232, 4234, 4235, 4236, 4237, 4239, 4240, 4242, 4243, 4246,
        4247, 4248, 4249, 4250, 4251, 4252, 4253, 4254, 4255, 4257, 4258, 4259, 4260, 4261, 4262, 4263, 4264, 4265,
        4266, 4267, 4268, 4269, 4270, 4271, 4272, 4273, 4310, 4312, 4313, 4316, 4318, 4319, 4331, 4332, 4333, 4335,
        4338, 4346, 4347, 4348, 4349, 4351, 4353, 4354, 4355, 4356, 4357, 4358, 4359, 4360, 4361, 4362, 4363, 4364,
        4365, 4366, 4367, 4368, 4369, 4406, 4408, 4409, 4412, 4414, 4415, 4416, 4421, 4437, 4448, 4452, 4453, 4454,
        4455, 4456, 4457, 4458, 4460, 4463, 4502, 4508, 4513, 4516, 4533, 4534, 4535, 4536, 4537, 4542, 4544, 4548,
        4549, 4550, 4551, 4552, 4553, 4554, 4556, 4559, 4598, 4604, 4610, 4611, 4622, 4625, 4628, 4630, 4631, 4632,
        4633, 4638, 4644, 4645, 4646, 4647, 4648, 4649, 4650, 4652, 4655, 4694, 4700, 4710, 4713, 4725, 4730, 4731,
        4732, 4733, 4735, 4736, 4737, 4739, 4740, 4741, 4742, 4743, 4744, 4745, 4746, 4748, 4749, 4751, 4752, 4790,
        4792, 4796, 4798, 4800, 4805, 4850, 4863, 4896, 4897, 4900, 4901, 4946, 4947, 4948, 4949, 4950, 4958, 4959,
        4992, 4993, 4994, 4995, 4996, 4997, 5006, 5009, 5012, 5042, 5043, 5044, 5045, 5046, 5054, 5055, 5089, 5090,
        5091, 5092, 5102, 5105, 5108, 5138, 5139, 5140, 5141, 5142, 5150, 5151, 5186, 5187, 5198, 5201, 5204, 5235,
        5236, 5237, 5238, 5246, 5280, 5285, 5286, 5289, 5330, 5335, 5336, 5337, 5338, 5339, 5340, 5341, 5343, 5344,
        5354, 5367, 5373, 5376, 5377, 5380, 5381, 5382, 5383, 5384, 5385, 5386, 5392, 5395, 5426, 5427, 5428, 5429,
        5430, 5431, 5432, 5433, 5434, 5435, 5436, 5437, 5438, 5439, 5440, 5450, 5463, 5469, 5473, 5474, 5475, 5476,
        5478, 5479, 5480, 5481, 5482, 5486, 5488, 5489, 5491, 5492, 5522, 5523, 5524, 5525, 5526, 5527, 5528, 5529,
        5530, 5531, 5532, 5533, 5534, 5535, 5536, 5546, 5559, 5565, 5570, 5571, 5575, 5576, 5578, 5582, 5584, 5585,
        5587, 5588, 5619, 5620, 5621, 5622, 5623, 5624, 5625, 5626, 5627, 5628, 5629, 5630, 5632, 5642, 5655, 5661,
        5664, 5669, 5670, 5673, 5675, 5676, 5677, 5679, 5682, 5698, 5707, 5710, 5713, 5714, 5719, 5720, 5721, 5722,
        5723, 5724, 5725, 5727, 5728, 5738, 5751, 5753, 5757, 5759, 5761, 5764, 5766, 5767, 5768, 5769, 5770, 5771,
        5772, 5773, 5775, 5776, 5778, 5779, 5794, 5803, 5806, 5809, 5810, 5811, 5812, 5813, 5814, 5815, 5816, 5817,
        5818, 5819, 5820, 5821, 5822, 5823, 5824, 5834, 5847, 5849, 5853, 5855, 5858, 5859, 5863, 5864, 5866, 5867,
        5868, 5869, 5870, 5871, 5872, 5873, 5874, 5875, 5876, 5890, 5899, 5902, 5905, 5907, 5908, 5909, 5910, 5911,
        5912, 5913, 5914, 5915, 5916, 5917, 5918, 5920, 5930, 5943, 5945, 5949, 5951, 5958, 5961, 5963, 5964, 5965,
        5967, 5970, 5986, 5995, 5998, 6001, 6002, 6007, 6008, 6009, 6010, 6011, 6012, 6013, 6015, 6016, 6026, 6039,
        6041, 6045, 6047, 6055, 6056, 6058, 6059, 6060, 6061, 6063, 6064, 6066, 6067, 6082, 6091, 6094, 6097, 6099,
        6100, 6101, 6102, 6103, 6104, 6105, 6106, 6107, 6108, 6109, 6110, 6112, 6122, 6135, 6137, 6141, 6143, 6155,
        6156, 6157, 6159, 6162, 6178, 6187, 6190, 6193, 6199, 6200, 6201, 6202, 6203, 6204, 6205, 6208, 6218, 6231,
        6233, 6237, 6239, 6240, 6245, 6261, 6272, 6290, 6303, 6305, 6306, 6307, 6308, 6309, 6310, 6311, 6312, 6313,
        6315, 6331, 6336, 6337, 6340, 6341, 6357, 6358, 6359, 6360, 6361, 6366, 6368, 6386, 6387, 6388, 6389, 6390,
        6398, 6399, 6401, 6402, 6403, 6404, 6405, 6406, 6407, 6408, 6409, 6411, 6427, 6433, 6434, 6435, 6436, 6446,
        6449, 6452, 6453, 6454, 6455, 6456, 6457, 6462, 6464, 6482, 6483, 6484, 6485, 6486, 6494, 6495, 6497, 6498,
        6499, 6500, 6501, 6502, 6503, 6504, 6505, 6507, 6523, 6530, 6531, 6542, 6545, 6548, 6550, 6551, 6552, 6553,
        6558, 6579, 6580, 6581, 6582, 6590, 6593, 6594, 6595, 6596, 6597, 6598, 6599, 6600, 6601, 6603, 6619, 6624,
        6629, 6630, 6633, 6645, 6650, 6651, 6652, 6653, 6655, 6656, 6657, 6659, 6669, 6672, 6674, 6679, 6680, 6681,
        6682, 6683, 6684, 6685, 6687, 6688, 6689, 6690, 6691, 6692, 6693, 6694, 6695, 6696, 6697, 6698, 6699, 6711,
        6712, 6715, 6717, 6718, 6721, 6724, 6726, 6727, 6728, 6729, 6730, 6736, 6739, 6741, 6742, 6743, 6744, 6745,
        6746, 6747, 6748, 6749, 6750, 6751, 6752, 6753, 6755, 6765, 6768, 6770, 6771, 6772, 6773, 6774, 6775, 6776,
        6777, 6778, 6779, 6780, 6781, 6782, 6783, 6784, 6785, 6786, 6787, 6788, 6789, 6790, 6791, 6792, 6793, 6794,
        6795, 6807, 6808, 6811, 6813, 6814, 6818, 6819, 6823, 6824, 6826, 6830, 6832, 6833, 6835, 6836, 6838, 6839,
        6840, 6841, 6842, 6843, 6844, 6845, 6846, 6847, 6849, 6851, 6861, 6864, 6867, 6868, 6869, 6870, 6871, 6872,
        6873, 6874, 6875, 6876, 6877, 6878, 6880, 6881, 6882, 6883, 6884, 6885, 6886, 6887, 6888, 6889, 6890, 6891,
        6903, 6904, 6907, 6909, 6910, 6918, 6921, 6923, 6924, 6925, 6927, 6930, 6933, 6938, 6939, 6940, 6941, 6943,
        6944, 6945, 6946, 6947, 6955, 6957, 6958, 6960, 6961, 6962, 6967, 6968, 6969, 6970, 6971, 6972, 6973, 6975,
        6976, 6977, 6978, 6979, 6980, 6981, 6982, 6983, 6984, 6985, 6986, 6987, 6999, 7000, 7001, 7003, 7005, 7006,
        7007, 7015, 7016, 7018, 7019, 7020, 7021, 7023, 7024, 7026, 7027, 7030, 7031, 7032, 7033, 7034, 7035, 7036,
        7037, 7038, 7039, 7041, 7042, 7043, 7051, 7053, 7054, 7056, 7057, 7059, 7060, 7061, 7062, 7063, 7064, 7065,
        7066, 7067, 7068, 7069, 7070, 7072, 7073, 7074, 7075, 7076, 7077, 7078, 7079, 7080, 7081, 7082, 7083, 7095,
        7096, 7097, 7099, 7101, 7102, 7103, 7115, 7116, 7117, 7119, 7122, 7130, 7131, 7132, 7133, 7135, 7137, 7138,
        7139, 7147, 7149, 7150, 7152, 7153, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7168, 7169, 7170, 7171, 7172,
        7173, 7174, 7175, 7176, 7177, 7178, 7179, 7191, 7192, 7193, 7195, 7197, 7198, 7199, 7200, 7205, 7221, 7232,
        7236, 7237, 7238, 7239, 7240, 7241, 7242, 7244, 7247, 7250, 7263, 7265, 7266, 7267, 7268, 7269, 7270, 7271,
        7272, 7273, 7275, 7286, 7291, 7292, 7296, 7301, 7346, 7359, 7372, 7373, 7374, 7375, 7376, 7377, 7378, 7379,
        7380, 7381, 7386, 7392, 7393, 7396, 7397, 7442, 7443, 7444, 7445, 7446, 7454, 7455, 7468, 7469, 7470, 7471,
        7472, 7473, 7474, 7475, 7476, 7477, 7482, 7489, 7490, 7491, 7492, 7502, 7505, 7508, 7538, 7539, 7540, 7541,
        7542, 7550, 7551, 7564, 7565, 7566, 7567, 7568, 7569, 7570, 7571, 7572, 7573, 7578, 7586, 7587, 7598, 7601,
        7604, 7635, 7636, 7637, 7638, 7646, 7660, 7661, 7662, 7663, 7664, 7665, 7666, 7667, 7668, 7669, 7674, 7680,
        7685, 7686, 7689, 7730, 7735, 7736, 7737, 7738, 7739, 7740, 7741, 7743, 7744, 7754, 7756, 7757, 7758, 7759,
        7760, 7761, 7762, 7763, 7764, 7765, 7767, 7770, 7773, 7777, 7780, 7782, 7783, 7784, 7785, 7786, 7792, 7795,
        7826, 7827, 7828, 7829, 7830, 7831, 7832, 7833, 7834, 7835, 7836, 7837, 7838, 7839, 7840, 7850, 7852, 7853,
        7854, 7855, 7856, 7857, 7858, 7859, 7860, 7861, 7863, 7866, 7869, 7874, 7875, 7879, 7880, 7882, 7886, 7888,
        7889, 7891, 7892, 7923, 7924, 7925, 7926, 7927, 7928, 7929, 7930, 7931, 7932, 7933, 7934, 7936, 7946, 7948,
        7949, 7950, 7951, 7952, 7953, 7954, 7955, 7956, 7957, 7959, 7962, 7965, 7974, 7977, 7979, 7980, 7981, 7983,
        7986, 8002, 8011, 8014, 8017, 8018, 8023, 8024, 8025, 8026, 8027, 8028, 8029, 8031, 8032, 8042, 8044, 8045,
        8046, 8047, 8048, 8049, 8050, 8051, 8052, 8053, 8055, 8057, 8058, 8061, 8063, 8071, 8072, 8074, 8075, 8076,
        8077, 8079, 8080, 8082, 8083, 8098, 8107, 8110, 8113, 8115, 8116, 8117, 8118, 8119, 8120, 8121, 8122, 8123,
        8124, 8125, 8126, 8128, 8138, 8140, 8141, 8142, 8143, 8144, 8145, 8146, 8147, 8148, 8149, 8151, 8153, 8154,
        8157, 8159, 8171, 8172, 8173, 8175, 8178, 8194, 8203, 8206, 8209, 8215, 8216, 8217, 8218, 8219, 8220, 8221,
        8224, 8234, 8236, 8237, 8238, 8239, 8240, 8241, 8242, 8243, 8244, 8245, 8247, 8249, 8250, 8253, 8255, 8256,
        8261, 8277, 8288, 8306, 8319, 8321, 8322, 8323, 8324, 8325, 8326, 8327, 8328, 8329, 8331, 8332, 8333, 8334,
        8335, 8336, 8337, 8338, 8339, 8340, 8341, 8346, 8347, 8353, 8356, 8373, 8374, 8375, 8376, 8377, 8382, 8384,
        8388, 8389, 8390, 8391, 8392, 8393, 8394, 8396, 8399, 8402, 8403, 8404, 8405, 8406, 8414, 8415, 8417, 8418,
        8419, 8420, 8421, 8422, 8423, 8424, 8425, 8427, 8438, 8443, 8444, 8454, 8457, 8469, 8474, 8475, 8476, 8477,
        8479, 8480, 8481, 8483, 8484, 8485, 8486, 8487, 8488, 8489, 8490, 8492, 8493, 8495, 8496, 8498, 8503, 8504,
        8505, 8506, 8507, 8508, 8509, 8511, 8512, 8513, 8514, 8515, 8516, 8517, 8518, 8519, 8520, 8521, 8522, 8523,
        8534, 8535, 8536, 8539, 8540, 8541, 8542, 8565, 8576, 8580, 8581, 8582, 8583, 8584, 8585, 8586, 8588, 8591,
        8630, 8636, 8642, 8643, 8654, 8657, 8660, 8662, 8663, 8664, 8665, 8670, 8676, 8677, 8678, 8679, 8680, 8681,
        8682, 8684, 8687, 8691, 8692, 8693, 8694, 8702, 8705, 8706, 8707, 8708, 8709, 8710, 8711, 8712, 8713, 8715,
        8726, 8731, 8732, 8743, 8744, 8746, 8752, 8755, 8758, 8759, 8760, 8761, 8762, 8763, 8764, 8765, 8766, 8767,
        8769, 8771, 8772, 8773, 8774, 8775, 8776, 8777, 8778, 8780, 8781, 8783, 8784, 8822, 8824, 8828, 8830, 8854,
        8855, 8856, 8857, 8862, 8868, 8869, 8870, 8871, 8872, 8873, 8874, 8876, 8879, 8918, 8924, 8939, 8940, 8941,
        8943, 8946, 8954, 8955, 8956, 8957, 8959, 8961, 8962, 8963, 8964, 8965, 8966, 8967, 8968, 8969, 8970, 8971,
        8972, 8973, 8974, 8975, 8976, 8977, 9014, 9016, 9017, 9020, 9022, 9023, 9050, 9051, 9052, 9053, 9055, 9057,
        9059, 9060, 9061, 9062, 9063, 9064, 9065, 9066, 9068, 9069, 9071, 9072, 9110, 9112, 9116, 9118, 9156, 9157,
        9158, 9159, 9160, 9161, 9162, 9164, 9167, 9206, 9212};
    static const int C1_ind[] = {
        76,   77,   78,   79,   80,   81,   82,   83,   84,   85,   90,   146,  159,  172,  173,  174,  175,  176,
        177,  178,  179,  180,  181,  186,  192,  197,  242,  255,  268,  269,  270,  271,  272,  273,  274,  275,
        276,  277,  282,  289,  292,  338,  339,  340,  341,  342,  350,  351,  364,  365,  366,  367,  368,  369,
        370,  371,  372,  373,  378,  385,  388,  405,  406,  407,  408,  409,  414,  416,  434,  435,  436,  437,
        438,  446,  447,  449,  450,  451,  452,  453,  454,  455,  456,  457,  459,  460,  461,  462,  463,  464,
        465,  466,  467,  468,  469,  474,  475,  486,  489,  530,  535,  536,  537,  538,  539,  540,  541,  543,
        544,  554,  556,  557,  558,  559,  560,  561,  562,  563,  564,  565,  567,  570,  573,  582,  585,  597,
        602,  603,  604,  605,  607,  608,  609,  611,  621,  624,  626,  631,  632,  633,  634,  635,  636,  637,
        639,  640,  641,  642,  643,  644,  645,  646,  647,  648,  649,  650,  651,  652,  653,  654,  655,  656,
        657,  658,  659,  660,  661,  663,  664,  666,  667,  669,  670,  693,  704,  722,  735,  737,  738,  739,
        740,  741,  742,  743,  744,  745,  747,  748,  749,  750,  751,  752,  753,  754,  755,  756,  757,  762,
        763,  789,  800,  804,  805,  806,  807,  808,  809,  810,  812,  815,  818,  831,  833,  834,  835,  836,
        837,  838,  839,  840,  841,  843,  844,  845,  846,  847,  848,  849,  850,  851,  852,  853,  854,  858,
        859,  860,  885,  896,  900,  901,  902,  903,  904,  905,  906,  908,  911,  914,  927,  929,  930,  931,
        932,  933,  934,  935,  936,  937,  939,  950,  955,  956,  1011, 1012, 1013, 1014, 1022, 1036, 1037, 1038,
        1039, 1040, 1041, 1042, 1043, 1044, 1045, 1050, 1058, 1059, 1070, 1073, 1076, 1107, 1108, 1109, 1110, 1118,
        1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1146, 1154, 1155, 1166, 1169, 1172, 1174, 1175,
        1176, 1177, 1182, 1203, 1204, 1205, 1206, 1214, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 1227,
        1228, 1229, 1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, 1242, 1243, 1255, 1256, 1258, 1264, 1267, 1299,
        1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308, 1309, 1310, 1312, 1322, 1324, 1325, 1326, 1327, 1328,
        1329, 1330, 1331, 1332, 1333, 1335, 1338, 1341, 1351, 1352, 1354, 1360, 1363, 1366, 1367, 1368, 1369, 1370,
        1371, 1372, 1373, 1374, 1375, 1377, 1379, 1389, 1392, 1395, 1396, 1397, 1398, 1399, 1400, 1401, 1402, 1403,
        1404, 1405, 1406, 1408, 1409, 1410, 1411, 1412, 1413, 1414, 1415, 1416, 1417, 1418, 1419, 1420, 1421, 1422,
        1423, 1424, 1425, 1426, 1427, 1428, 1429, 1431, 1432, 1434, 1435, 1437, 1438, 1447, 1448, 1450, 1456, 1459,
        1462, 1463, 1464, 1465, 1466, 1467, 1468, 1469, 1470, 1471, 1473, 1475, 1476, 1477, 1478, 1479, 1480, 1481,
        1482, 1484, 1485, 1487, 1488, 1491, 1492, 1493, 1494, 1495, 1496, 1497, 1498, 1499, 1500, 1501, 1502, 1504,
        1505, 1506, 1507, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 1526, 1527, 1528, 1531, 1532, 1533, 1534,
        1558, 1559, 1560, 1561, 1566, 1587, 1588, 1589, 1590, 1598, 1601, 1602, 1603, 1604, 1605, 1606, 1607, 1608,
        1609, 1611, 1612, 1613, 1614, 1615, 1616, 1617, 1618, 1619, 1620, 1621, 1626, 1627, 1654, 1655, 1656, 1657,
        1662, 1668, 1669, 1670, 1671, 1672, 1673, 1674, 1676, 1679, 1683, 1684, 1685, 1686, 1694, 1697, 1698, 1699,
        1700, 1701, 1702, 1703, 1704, 1705, 1707, 1708, 1709, 1710, 1711, 1712, 1713, 1714, 1715, 1716, 1717, 1718,
        1722, 1723, 1724, 1750, 1751, 1752, 1753, 1758, 1764, 1765, 1766, 1767, 1768, 1769, 1770, 1772, 1775, 1779,
        1780, 1781, 1782, 1790, 1793, 1794, 1795, 1796, 1797, 1798, 1799, 1800, 1801, 1803, 1814, 1819, 1820, 1879,
        1880, 1881, 1882, 1883, 1884, 1885, 1888, 1898, 1900, 1901, 1902, 1903, 1904, 1905, 1906, 1907, 1908, 1909,
        1911, 1914, 1917, 1931, 1932, 1933, 1935, 1938, 1954, 1963, 1966, 1969, 1975, 1976, 1977, 1978, 1979, 1980,
        1981, 1984, 1994, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2007, 2009, 2010, 2013, 2015,
        2027, 2028, 2029, 2031, 2034, 2042, 2043, 2044, 2045, 2047, 2049, 2050, 2051, 2059, 2061, 2062, 2064, 2065,
        2071, 2072, 2073, 2074, 2075, 2076, 2077, 2080, 2081, 2082, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090,
        2091, 2092, 2093, 2094, 2095, 2096, 2097, 2098, 2099, 2100, 2101, 2103, 2104, 2105, 2106, 2107, 2109, 2110,
        2111, 2123, 2124, 2125, 2127, 2130, 2138, 2139, 2140, 2141, 2143, 2145, 2146, 2147, 2148, 2149, 2150, 2151,
        2152, 2153, 2154, 2155, 2156, 2157, 2158, 2159, 2160, 2161, 2167, 2168, 2169, 2170, 2171, 2172, 2173, 2176,
        2177, 2178, 2179, 2180, 2181, 2182, 2183, 2184, 2185, 2186, 2187, 2198, 2199, 2200, 2201, 2203, 2204, 2205,
        2206, 2207, 2234, 2235, 2236, 2237, 2239, 2241, 2243, 2253, 2256, 2263, 2264, 2265, 2266, 2267, 2268, 2269,
        2272, 2273, 2274, 2275, 2276, 2277, 2278, 2279, 2280, 2281, 2282, 2283, 2284, 2285, 2286, 2287, 2288, 2289,
        2290, 2291, 2292, 2293, 2295, 2296, 2298, 2299, 2301, 2302, 2330, 2331, 2332, 2333, 2335, 2337, 2339, 2340,
        2341, 2342, 2343, 2344, 2345, 2346, 2348, 2349, 2351, 2352, 2359, 2360, 2361, 2362, 2363, 2364, 2365, 2368,
        2369, 2370, 2371, 2372, 2373, 2374, 2375, 2376, 2377, 2378, 2379, 2380, 2381, 2382, 2383, 2384, 2385, 2386,
        2387, 2388, 2389, 2390, 2391, 2392, 2394, 2395, 2396, 2397, 2398, 2426, 2427, 2428, 2429, 2431, 2433, 2435,
        2436, 2437, 2438, 2439, 2440, 2441, 2442, 2444, 2445, 2447, 2448, 2455, 2456, 2457, 2458, 2459, 2460, 2461,
        2464, 2465, 2466, 2467, 2468, 2469, 2470, 2471, 2472, 2473, 2474, 2475, 2486, 2487, 2488, 2491, 2492, 2493,
        2494, 2561, 2562, 2563, 2564, 2565, 2566, 2567, 2568, 2569, 2571, 2572, 2573, 2574, 2575, 2576, 2577, 2578,
        2579, 2580, 2581, 2586, 2587, 2628, 2629, 2630, 2631, 2632, 2633, 2634, 2636, 2639, 2657, 2658, 2659, 2660,
        2661, 2662, 2663, 2664, 2665, 2667, 2668, 2669, 2670, 2671, 2672, 2673, 2674, 2675, 2676, 2677, 2678, 2682,
        2683, 2684, 2724, 2725, 2726, 2727, 2728, 2729, 2730, 2732, 2735, 2753, 2754, 2755, 2756, 2757, 2758, 2759,
        2760, 2761, 2763, 2764, 2765, 2766, 2767, 2768, 2769, 2770, 2771, 2772, 2773, 2774, 2778, 2779, 2780, 2820,
        2821, 2822, 2823, 2824, 2825, 2826, 2828, 2831, 2849, 2850, 2851, 2852, 2853, 2854, 2855, 2856, 2857, 2859,
        2870, 2875, 2876};

    Eigen::MatrixXd C0 = Eigen::MatrixXd::Zero(96, 96);
    Eigen::MatrixXd C1 = Eigen::MatrixXd::Zero(96, 30);
    for (int i = 0; i < 2349; i++) {
        C0(C0_ind[i]) = coeffs(coeffs0_ind[i]);
    }
    for (int i = 0; i < 1011; i++) {
        C1(C1_ind[i]) = coeffs(coeffs1_ind[i]);
    }

    Eigen::MatrixXd C12 = C0.partialPivLu().solve(C1);

    // Setup action matrix
    Eigen::Matrix<double, 40, 30> RR;
    RR << -C12.bottomRows(10), Eigen::Matrix<double, 30, 30>::Identity(30, 30);

    static const int AM_ind[] = {36, 17, 0,  14, 1,  16, 2,  18, 19, 3,  26, 22, 4,  24, 25,
                                 5,  27, 28, 6,  33, 31, 32, 7,  34, 35, 8,  37, 38, 39, 9};
    Eigen::Matrix<double, 30, 30> AM;
    for (int i = 0; i < 30; i++) {
        AM.row(i) = RR.row(AM_ind[i]);
    }

    // Eigen::Matrix<std::complex<double>, 4, 30> sols;
    sols.setZero();

    // Solve eigenvalue problem

    double p[1 + 30];
    Eigen::Matrix<double, 30, 30> AMp = AM;
    charpoly_danilevsky_piv(AMp, p);
    double roots[30];
    int nroots;
    find_real_roots_sturm(p, 30, roots, &nroots, 8, 0);
    fast_eigenvector_solver(roots, nroots, AM, sols);

    return nroots;
}
// Action =  w
// Quotient ring basis (V) =
// 1,x,x^2,x*y,x*y*w,x*z,x*z*w,x*w,x*w^2,x*w^3,y,y^2,y^2*w,y*z,y*z*w,y*z*w^2,y*w,y*w^2,y*w^3,z,z^2,z^2*w,z^2*w^2,z*w,z*w^2,z*w^3,w,w^2,w^3,w^4,
// Available monomials (RR*V) =
// x^2*w,x*y*w^2,x*z*w^2,x*w^4,y^2*w^2,y*z*w^3,y*w^4,z^2*w^3,z*w^4,w^5,1,x,x^2,x*y,x*y*w,x*z,x*z*w,x*w,x*w^2,x*w^3,y,y^2,y^2*w,y*z,y*z*w,y*z*w^2,y*w,y*w^2,y*w^3,z,z^2,z^2*w,z^2*w^2,z*w,z*w^2,z*w^3,w,w^2,w^3,w^4,

namespace poselib {

int relpose_4pt_planar(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                       std::vector<Eigen::Matrix3d> *essential_matrices) {

    // Compute nullspace to epipolar constraints
    Eigen::Matrix<double, 9, 4> epipolar_constraints;
    for (size_t i = 0; i < 4; ++i) {
        epipolar_constraints.col(i) << x1[i](0) * x2[i], x1[i](1) * x2[i], x1[i](2) * x2[i];
    }
    Eigen::Matrix<double, 9, 9> Q = epipolar_constraints.fullPivHouseholderQr().matrixQ();
    Eigen::Matrix<double, 9, 5> N = Q.rightCols(5);

    Eigen::VectorXd B(Eigen::Map<Eigen::VectorXd>(N.data(), N.cols() * N.rows()));

    Eigen::Matrix<std::complex<double>, 4, 30> sols;
    int n_sols = get_sols_relpose_4pt_planar(B, sols);

    essential_matrices->clear();
    essential_matrices->reserve(n_sols);

    for (int i = 0; i < n_sols; i++) {
        Eigen::Vector<double, 9> essential_matrix_vector = sols(0, i).real() * N.col(0) + sols(1, i).real() * N.col(1) + sols(2, i).real() * N.col(2) + sols(3, i).real() * N.col(3) + N.col(4);
        essential_matrix_vector.normalize();
        Eigen::Matrix3d essential_matrix = Eigen::Map<Eigen::Matrix3d>(essential_matrix_vector.data());
        essential_matrices->push_back(essential_matrix);
    }

    return n_sols;
}

int relpose_4pt_planar(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                       std::vector<CameraPose> *output) {
    std::vector<Eigen::Matrix3d> essential_matrices;
    int n_sols = relpose_4pt_planar(x1, x2, &essential_matrices);

    output->clear();
    output->reserve(n_sols);
    for (int i = 0; i < n_sols; ++i) {
        motion_from_essential(essential_matrices[i], x1[0], x2[0], output);
    }

    return output->size();
}

} // namespace poselib