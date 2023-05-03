#include "relpose_6pt_planar.h"
#include <Eigen/Dense>
#include <math.h>
#include <stdio.h>
#include <PoseLib/misc/essential.h>

#define RELERROR      1.0e-12   /* smallest relative error we want */
//#define MAXPOW        0        /* max power of 10 we wish to search to */
#define MAXIT         800       /* max number of iterations */
#define SMALL_ENOUGH  1.0e-12   /* a coefficient smaller than SMALL_ENOUGH 
/* is considered to be zero (0.0). */
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
namespace
{
	
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
} //namespace

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

void fast_eigenvector_solver(double *eigv, int neig, Eigen::Matrix<double, 9, 9> &AM,
                             Eigen::Matrix<std::complex<double>, 2, 9> &sols) {
    static const int ind[] = {2, 4, 6, 9, 12, 15, 18, 22, 25, 29};
    // Truncated action matrix containing non-trivial rows
    Eigen::Matrix<double, 10, 9> AMs;
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

int solver_6pt_planar(const Eigen::VectorXd &data, Eigen::Matrix<std::complex<double>, 2, 9> &sols)
{
	// Compute coefficients
    const double* d = data.data();
    Eigen::VectorXd coeffs(20);
    coeffs[0] = -d[2]*d[4]*d[6] + d[1]*d[5]*d[6] + d[2]*d[3]*d[7] - d[0]*d[5]*d[7] - d[1]*d[3]*d[8] + d[0]*d[4]*d[8];
    coeffs[1] = -d[5]*d[7]*d[9] + d[4]*d[8]*d[9] + d[5]*d[6]*d[10] - d[3]*d[8]*d[10] - d[4]*d[6]*d[11] + d[3]*d[7]*d[11] + d[2]*d[7]*d[12] - d[1]*d[8]*d[12] - d[2]*d[6]*d[13] + d[0]*d[8]*d[13] + d[1]*d[6]*d[14] - d[0]*d[7]*d[14] - d[2]*d[4]*d[15] + d[1]*d[5]*d[15] + d[2]*d[3]*d[16] - d[0]*d[5]*d[16] - d[1]*d[3]*d[17] + d[0]*d[4]*d[17];
    coeffs[2] = -d[8]*d[10]*d[12] + d[7]*d[11]*d[12] + d[8]*d[9]*d[13] - d[6]*d[11]*d[13] - d[7]*d[9]*d[14] + d[6]*d[10]*d[14] + d[5]*d[10]*d[15] - d[4]*d[11]*d[15] - d[2]*d[13]*d[15] + d[1]*d[14]*d[15] - d[5]*d[9]*d[16] + d[3]*d[11]*d[16] + d[2]*d[12]*d[16] - d[0]*d[14]*d[16] + d[4]*d[9]*d[17] - d[3]*d[10]*d[17] - d[1]*d[12]*d[17] + d[0]*d[13]*d[17];
    coeffs[3] = -d[11]*d[13]*d[15] + d[10]*d[14]*d[15] + d[11]*d[12]*d[16] - d[9]*d[14]*d[16] - d[10]*d[12]*d[17] + d[9]*d[13]*d[17];
    coeffs[4] = -d[5]*d[7]*d[18] + d[4]*d[8]*d[18] + d[5]*d[6]*d[19] - d[3]*d[8]*d[19] - d[4]*d[6]*d[20] + d[3]*d[7]*d[20] + d[2]*d[7]*d[21] - d[1]*d[8]*d[21] - d[2]*d[6]*d[22] + d[0]*d[8]*d[22] + d[1]*d[6]*d[23] - d[0]*d[7]*d[23] - d[2]*d[4]*d[24] + d[1]*d[5]*d[24] + d[2]*d[3]*d[25] - d[0]*d[5]*d[25] - d[1]*d[3]*d[26] + d[0]*d[4]*d[26];
    coeffs[5] = d[8]*d[13]*d[18] - d[7]*d[14]*d[18] - d[5]*d[16]*d[18] + d[4]*d[17]*d[18] - d[8]*d[12]*d[19] + d[6]*d[14]*d[19] + d[5]*d[15]*d[19] - d[3]*d[17]*d[19] + d[7]*d[12]*d[20] - d[6]*d[13]*d[20] - d[4]*d[15]*d[20] + d[3]*d[16]*d[20] - d[8]*d[10]*d[21] + d[7]*d[11]*d[21] + d[2]*d[16]*d[21] - d[1]*d[17]*d[21] + d[8]*d[9]*d[22] - d[6]*d[11]*d[22] - d[2]*d[15]*d[22] + d[0]*d[17]*d[22] - d[7]*d[9]*d[23] + d[6]*d[10]*d[23] + d[1]*d[15]*d[23] - d[0]*d[16]*d[23] + d[5]*d[10]*d[24] - d[4]*d[11]*d[24] - d[2]*d[13]*d[24] + d[1]*d[14]*d[24] - d[5]*d[9]*d[25] + d[3]*d[11]*d[25] + d[2]*d[12]*d[25] - d[0]*d[14]*d[25] + d[4]*d[9]*d[26] - d[3]*d[10]*d[26] - d[1]*d[12]*d[26] + d[0]*d[13]*d[26];
    coeffs[6] = -d[14]*d[16]*d[18] + d[13]*d[17]*d[18] + d[14]*d[15]*d[19] - d[12]*d[17]*d[19] - d[13]*d[15]*d[20] + d[12]*d[16]*d[20] + d[11]*d[16]*d[21] - d[10]*d[17]*d[21] - d[11]*d[15]*d[22] + d[9]*d[17]*d[22] + d[10]*d[15]*d[23] - d[9]*d[16]*d[23] - d[11]*d[13]*d[24] + d[10]*d[14]*d[24] + d[11]*d[12]*d[25] - d[9]*d[14]*d[25] - d[10]*d[12]*d[26] + d[9]*d[13]*d[26];
    coeffs[7] = -d[8]*d[19]*d[21] + d[7]*d[20]*d[21] + d[8]*d[18]*d[22] - d[6]*d[20]*d[22] - d[7]*d[18]*d[23] + d[6]*d[19]*d[23] + d[5]*d[19]*d[24] - d[4]*d[20]*d[24] - d[2]*d[22]*d[24] + d[1]*d[23]*d[24] - d[5]*d[18]*d[25] + d[3]*d[20]*d[25] + d[2]*d[21]*d[25] - d[0]*d[23]*d[25] + d[4]*d[18]*d[26] - d[3]*d[19]*d[26] - d[1]*d[21]*d[26] + d[0]*d[22]*d[26];
    coeffs[8] = -d[17]*d[19]*d[21] + d[16]*d[20]*d[21] + d[17]*d[18]*d[22] - d[15]*d[20]*d[22] - d[16]*d[18]*d[23] + d[15]*d[19]*d[23] + d[14]*d[19]*d[24] - d[13]*d[20]*d[24] - d[11]*d[22]*d[24] + d[10]*d[23]*d[24] - d[14]*d[18]*d[25] + d[12]*d[20]*d[25] + d[11]*d[21]*d[25] - d[9]*d[23]*d[25] + d[13]*d[18]*d[26] - d[12]*d[19]*d[26] - d[10]*d[21]*d[26] + d[9]*d[22]*d[26];
    coeffs[9] = -d[20]*d[22]*d[24] + d[19]*d[23]*d[24] + d[20]*d[21]*d[25] - d[18]*d[23]*d[25] - d[19]*d[21]*d[26] + d[18]*d[22]*d[26];
    coeffs[10] = -2*std::pow(d[2],2)*d[4] + 2*d[1]*d[2]*d[5] + 2*d[2]*d[3]*d[5] - 2*d[0]*std::pow(d[5],2) - 4*d[2]*d[4]*d[6] + 2*d[1]*d[5]*d[6] + 2*d[3]*d[5]*d[6] - 2*d[4]*std::pow(d[6],2) + 2*d[1]*d[2]*d[7] + 2*d[2]*d[3]*d[7] - 4*d[0]*d[5]*d[7] + 2*d[1]*d[6]*d[7] + 2*d[3]*d[6]*d[7] - 2*d[0]*std::pow(d[7],2) - 2*std::pow(d[1],2)*d[8] - 4*d[1]*d[3]*d[8] - 2*std::pow(d[3],2)*d[8] + 8*d[0]*d[4]*d[8];
    coeffs[11] = -2*std::pow(d[5],2)*d[9] - 4*d[5]*d[7]*d[9] - 2*std::pow(d[7],2)*d[9] + 8*d[4]*d[8]*d[9] + 2*d[2]*d[5]*d[10] + 2*d[5]*d[6]*d[10] + 2*d[2]*d[7]*d[10] + 2*d[6]*d[7]*d[10] - 4*d[1]*d[8]*d[10] - 4*d[3]*d[8]*d[10] - 4*d[2]*d[4]*d[11] + 2*d[1]*d[5]*d[11] + 2*d[3]*d[5]*d[11] - 4*d[4]*d[6]*d[11] + 2*d[1]*d[7]*d[11] + 2*d[3]*d[7]*d[11] + 2*d[2]*d[5]*d[12] + 2*d[5]*d[6]*d[12] + 2*d[2]*d[7]*d[12] + 2*d[6]*d[7]*d[12] - 4*d[1]*d[8]*d[12] - 4*d[3]*d[8]*d[12] - 2*std::pow(d[2],2)*d[13] - 4*d[2]*d[6]*d[13] - 2*std::pow(d[6],2)*d[13] + 8*d[0]*d[8]*d[13] + 2*d[1]*d[2]*d[14] + 2*d[2]*d[3]*d[14] - 4*d[0]*d[5]*d[14] + 2*d[1]*d[6]*d[14] + 2*d[3]*d[6]*d[14] - 4*d[0]*d[7]*d[14] - 4*d[2]*d[4]*d[15] + 2*d[1]*d[5]*d[15] + 2*d[3]*d[5]*d[15] - 4*d[4]*d[6]*d[15] + 2*d[1]*d[7]*d[15] + 2*d[3]*d[7]*d[15] + 2*d[1]*d[2]*d[16] + 2*d[2]*d[3]*d[16] - 4*d[0]*d[5]*d[16] + 2*d[1]*d[6]*d[16] + 2*d[3]*d[6]*d[16] - 4*d[0]*d[7]*d[16] - 2*std::pow(d[1],2)*d[17] - 4*d[1]*d[3]*d[17] - 2*std::pow(d[3],2)*d[17] + 8*d[0]*d[4]*d[17];
    coeffs[12] = -2*d[8]*std::pow(d[10],2) + 2*d[5]*d[10]*d[11] + 2*d[7]*d[10]*d[11] - 2*d[4]*std::pow(d[11],2) - 4*d[8]*d[10]*d[12] + 2*d[5]*d[11]*d[12] + 2*d[7]*d[11]*d[12] - 2*d[8]*std::pow(d[12],2) + 8*d[8]*d[9]*d[13] - 4*d[2]*d[11]*d[13] - 4*d[6]*d[11]*d[13] - 4*d[5]*d[9]*d[14] - 4*d[7]*d[9]*d[14] + 2*d[2]*d[10]*d[14] + 2*d[6]*d[10]*d[14] + 2*d[1]*d[11]*d[14] + 2*d[3]*d[11]*d[14] + 2*d[2]*d[12]*d[14] + 2*d[6]*d[12]*d[14] - 2*d[0]*std::pow(d[14],2) + 2*d[5]*d[10]*d[15] + 2*d[7]*d[10]*d[15] - 4*d[4]*d[11]*d[15] + 2*d[5]*d[12]*d[15] + 2*d[7]*d[12]*d[15] - 4*d[2]*d[13]*d[15] - 4*d[6]*d[13]*d[15] + 2*d[1]*d[14]*d[15] + 2*d[3]*d[14]*d[15] - 2*d[4]*std::pow(d[15],2) - 4*d[5]*d[9]*d[16] - 4*d[7]*d[9]*d[16] + 2*d[2]*d[10]*d[16] + 2*d[6]*d[10]*d[16] + 2*d[1]*d[11]*d[16] + 2*d[3]*d[11]*d[16] + 2*d[2]*d[12]*d[16] + 2*d[6]*d[12]*d[16] - 4*d[0]*d[14]*d[16] + 2*d[1]*d[15]*d[16] + 2*d[3]*d[15]*d[16] - 2*d[0]*std::pow(d[16],2) + 8*d[4]*d[9]*d[17] - 4*d[1]*d[10]*d[17] - 4*d[3]*d[10]*d[17] - 4*d[1]*d[12]*d[17] - 4*d[3]*d[12]*d[17] + 8*d[0]*d[13]*d[17];
    coeffs[13] = -2*std::pow(d[11],2)*d[13] + 2*d[10]*d[11]*d[14] + 2*d[11]*d[12]*d[14] - 2*d[9]*std::pow(d[14],2) - 4*d[11]*d[13]*d[15] + 2*d[10]*d[14]*d[15] + 2*d[12]*d[14]*d[15] - 2*d[13]*std::pow(d[15],2) + 2*d[10]*d[11]*d[16] + 2*d[11]*d[12]*d[16] - 4*d[9]*d[14]*d[16] + 2*d[10]*d[15]*d[16] + 2*d[12]*d[15]*d[16] - 2*d[9]*std::pow(d[16],2) - 2*std::pow(d[10],2)*d[17] - 4*d[10]*d[12]*d[17] - 2*std::pow(d[12],2)*d[17] + 8*d[9]*d[13]*d[17];
    coeffs[14] = -2*std::pow(d[5],2)*d[18] - 4*d[5]*d[7]*d[18] - 2*std::pow(d[7],2)*d[18] + 8*d[4]*d[8]*d[18] + 2*d[2]*d[5]*d[19] + 2*d[5]*d[6]*d[19] + 2*d[2]*d[7]*d[19] + 2*d[6]*d[7]*d[19] - 4*d[1]*d[8]*d[19] - 4*d[3]*d[8]*d[19] - 4*d[2]*d[4]*d[20] + 2*d[1]*d[5]*d[20] + 2*d[3]*d[5]*d[20] - 4*d[4]*d[6]*d[20] + 2*d[1]*d[7]*d[20] + 2*d[3]*d[7]*d[20] + 2*d[2]*d[5]*d[21] + 2*d[5]*d[6]*d[21] + 2*d[2]*d[7]*d[21] + 2*d[6]*d[7]*d[21] - 4*d[1]*d[8]*d[21] - 4*d[3]*d[8]*d[21] - 2*std::pow(d[2],2)*d[22] - 4*d[2]*d[6]*d[22] - 2*std::pow(d[6],2)*d[22] + 8*d[0]*d[8]*d[22] + 2*d[1]*d[2]*d[23] + 2*d[2]*d[3]*d[23] - 4*d[0]*d[5]*d[23] + 2*d[1]*d[6]*d[23] + 2*d[3]*d[6]*d[23] - 4*d[0]*d[7]*d[23] - 4*d[2]*d[4]*d[24] + 2*d[1]*d[5]*d[24] + 2*d[3]*d[5]*d[24] - 4*d[4]*d[6]*d[24] + 2*d[1]*d[7]*d[24] + 2*d[3]*d[7]*d[24] + 2*d[1]*d[2]*d[25] + 2*d[2]*d[3]*d[25] - 4*d[0]*d[5]*d[25] + 2*d[1]*d[6]*d[25] + 2*d[3]*d[6]*d[25] - 4*d[0]*d[7]*d[25] - 2*std::pow(d[1],2)*d[26] - 4*d[1]*d[3]*d[26] - 2*std::pow(d[3],2)*d[26] + 8*d[0]*d[4]*d[26];
    coeffs[15] = 8*d[8]*d[13]*d[18] - 4*d[5]*d[14]*d[18] - 4*d[7]*d[14]*d[18] - 4*d[5]*d[16]*d[18] - 4*d[7]*d[16]*d[18] + 8*d[4]*d[17]*d[18] - 4*d[8]*d[10]*d[19] + 2*d[5]*d[11]*d[19] + 2*d[7]*d[11]*d[19] - 4*d[8]*d[12]*d[19] + 2*d[2]*d[14]*d[19] + 2*d[6]*d[14]*d[19] + 2*d[5]*d[15]*d[19] + 2*d[7]*d[15]*d[19] + 2*d[2]*d[16]*d[19] + 2*d[6]*d[16]*d[19] - 4*d[1]*d[17]*d[19] - 4*d[3]*d[17]*d[19] + 2*d[5]*d[10]*d[20] + 2*d[7]*d[10]*d[20] - 4*d[4]*d[11]*d[20] + 2*d[5]*d[12]*d[20] + 2*d[7]*d[12]*d[20] - 4*d[2]*d[13]*d[20] - 4*d[6]*d[13]*d[20] + 2*d[1]*d[14]*d[20] + 2*d[3]*d[14]*d[20] - 4*d[4]*d[15]*d[20] + 2*d[1]*d[16]*d[20] + 2*d[3]*d[16]*d[20] - 4*d[8]*d[10]*d[21] + 2*d[5]*d[11]*d[21] + 2*d[7]*d[11]*d[21] - 4*d[8]*d[12]*d[21] + 2*d[2]*d[14]*d[21] + 2*d[6]*d[14]*d[21] + 2*d[5]*d[15]*d[21] + 2*d[7]*d[15]*d[21] + 2*d[2]*d[16]*d[21] + 2*d[6]*d[16]*d[21] - 4*d[1]*d[17]*d[21] - 4*d[3]*d[17]*d[21] + 8*d[8]*d[9]*d[22] - 4*d[2]*d[11]*d[22] - 4*d[6]*d[11]*d[22] - 4*d[2]*d[15]*d[22] - 4*d[6]*d[15]*d[22] + 8*d[0]*d[17]*d[22] - 4*d[5]*d[9]*d[23] - 4*d[7]*d[9]*d[23] + 2*d[2]*d[10]*d[23] + 2*d[6]*d[10]*d[23] + 2*d[1]*d[11]*d[23] + 2*d[3]*d[11]*d[23] + 2*d[2]*d[12]*d[23] + 2*d[6]*d[12]*d[23] - 4*d[0]*d[14]*d[23] + 2*d[1]*d[15]*d[23] + 2*d[3]*d[15]*d[23] - 4*d[0]*d[16]*d[23] + 2*d[5]*d[10]*d[24] + 2*d[7]*d[10]*d[24] - 4*d[4]*d[11]*d[24] + 2*d[5]*d[12]*d[24] + 2*d[7]*d[12]*d[24] - 4*d[2]*d[13]*d[24] - 4*d[6]*d[13]*d[24] + 2*d[1]*d[14]*d[24] + 2*d[3]*d[14]*d[24] - 4*d[4]*d[15]*d[24] + 2*d[1]*d[16]*d[24] + 2*d[3]*d[16]*d[24] - 4*d[5]*d[9]*d[25] - 4*d[7]*d[9]*d[25] + 2*d[2]*d[10]*d[25] + 2*d[6]*d[10]*d[25] + 2*d[1]*d[11]*d[25] + 2*d[3]*d[11]*d[25] + 2*d[2]*d[12]*d[25] + 2*d[6]*d[12]*d[25] - 4*d[0]*d[14]*d[25] + 2*d[1]*d[15]*d[25] + 2*d[3]*d[15]*d[25] - 4*d[0]*d[16]*d[25] + 8*d[4]*d[9]*d[26] - 4*d[1]*d[10]*d[26] - 4*d[3]*d[10]*d[26] - 4*d[1]*d[12]*d[26] - 4*d[3]*d[12]*d[26] + 8*d[0]*d[13]*d[26];
    coeffs[16] = -2*std::pow(d[14],2)*d[18] - 4*d[14]*d[16]*d[18] - 2*std::pow(d[16],2)*d[18] + 8*d[13]*d[17]*d[18] + 2*d[11]*d[14]*d[19] + 2*d[14]*d[15]*d[19] + 2*d[11]*d[16]*d[19] + 2*d[15]*d[16]*d[19] - 4*d[10]*d[17]*d[19] - 4*d[12]*d[17]*d[19] - 4*d[11]*d[13]*d[20] + 2*d[10]*d[14]*d[20] + 2*d[12]*d[14]*d[20] - 4*d[13]*d[15]*d[20] + 2*d[10]*d[16]*d[20] + 2*d[12]*d[16]*d[20] + 2*d[11]*d[14]*d[21] + 2*d[14]*d[15]*d[21] + 2*d[11]*d[16]*d[21] + 2*d[15]*d[16]*d[21] - 4*d[10]*d[17]*d[21] - 4*d[12]*d[17]*d[21] - 2*std::pow(d[11],2)*d[22] - 4*d[11]*d[15]*d[22] - 2*std::pow(d[15],2)*d[22] + 8*d[9]*d[17]*d[22] + 2*d[10]*d[11]*d[23] + 2*d[11]*d[12]*d[23] - 4*d[9]*d[14]*d[23] + 2*d[10]*d[15]*d[23] + 2*d[12]*d[15]*d[23] - 4*d[9]*d[16]*d[23] - 4*d[11]*d[13]*d[24] + 2*d[10]*d[14]*d[24] + 2*d[12]*d[14]*d[24] - 4*d[13]*d[15]*d[24] + 2*d[10]*d[16]*d[24] + 2*d[12]*d[16]*d[24] + 2*d[10]*d[11]*d[25] + 2*d[11]*d[12]*d[25] - 4*d[9]*d[14]*d[25] + 2*d[10]*d[15]*d[25] + 2*d[12]*d[15]*d[25] - 4*d[9]*d[16]*d[25] - 2*std::pow(d[10],2)*d[26] - 4*d[10]*d[12]*d[26] - 2*std::pow(d[12],2)*d[26] + 8*d[9]*d[13]*d[26];
    coeffs[17] = -2*d[8]*std::pow(d[19],2) + 2*d[5]*d[19]*d[20] + 2*d[7]*d[19]*d[20] - 2*d[4]*std::pow(d[20],2) - 4*d[8]*d[19]*d[21] + 2*d[5]*d[20]*d[21] + 2*d[7]*d[20]*d[21] - 2*d[8]*std::pow(d[21],2) + 8*d[8]*d[18]*d[22] - 4*d[2]*d[20]*d[22] - 4*d[6]*d[20]*d[22] - 4*d[5]*d[18]*d[23] - 4*d[7]*d[18]*d[23] + 2*d[2]*d[19]*d[23] + 2*d[6]*d[19]*d[23] + 2*d[1]*d[20]*d[23] + 2*d[3]*d[20]*d[23] + 2*d[2]*d[21]*d[23] + 2*d[6]*d[21]*d[23] - 2*d[0]*std::pow(d[23],2) + 2*d[5]*d[19]*d[24] + 2*d[7]*d[19]*d[24] - 4*d[4]*d[20]*d[24] + 2*d[5]*d[21]*d[24] + 2*d[7]*d[21]*d[24] - 4*d[2]*d[22]*d[24] - 4*d[6]*d[22]*d[24] + 2*d[1]*d[23]*d[24] + 2*d[3]*d[23]*d[24] - 2*d[4]*std::pow(d[24],2) - 4*d[5]*d[18]*d[25] - 4*d[7]*d[18]*d[25] + 2*d[2]*d[19]*d[25] + 2*d[6]*d[19]*d[25] + 2*d[1]*d[20]*d[25] + 2*d[3]*d[20]*d[25] + 2*d[2]*d[21]*d[25] + 2*d[6]*d[21]*d[25] - 4*d[0]*d[23]*d[25] + 2*d[1]*d[24]*d[25] + 2*d[3]*d[24]*d[25] - 2*d[0]*std::pow(d[25],2) + 8*d[4]*d[18]*d[26] - 4*d[1]*d[19]*d[26] - 4*d[3]*d[19]*d[26] - 4*d[1]*d[21]*d[26] - 4*d[3]*d[21]*d[26] + 8*d[0]*d[22]*d[26];
    coeffs[18] = -2*d[17]*std::pow(d[19],2) + 2*d[14]*d[19]*d[20] + 2*d[16]*d[19]*d[20] - 2*d[13]*std::pow(d[20],2) - 4*d[17]*d[19]*d[21] + 2*d[14]*d[20]*d[21] + 2*d[16]*d[20]*d[21] - 2*d[17]*std::pow(d[21],2) + 8*d[17]*d[18]*d[22] - 4*d[11]*d[20]*d[22] - 4*d[15]*d[20]*d[22] - 4*d[14]*d[18]*d[23] - 4*d[16]*d[18]*d[23] + 2*d[11]*d[19]*d[23] + 2*d[15]*d[19]*d[23] + 2*d[10]*d[20]*d[23] + 2*d[12]*d[20]*d[23] + 2*d[11]*d[21]*d[23] + 2*d[15]*d[21]*d[23] - 2*d[9]*std::pow(d[23],2) + 2*d[14]*d[19]*d[24] + 2*d[16]*d[19]*d[24] - 4*d[13]*d[20]*d[24] + 2*d[14]*d[21]*d[24] + 2*d[16]*d[21]*d[24] - 4*d[11]*d[22]*d[24] - 4*d[15]*d[22]*d[24] + 2*d[10]*d[23]*d[24] + 2*d[12]*d[23]*d[24] - 2*d[13]*std::pow(d[24],2) - 4*d[14]*d[18]*d[25] - 4*d[16]*d[18]*d[25] + 2*d[11]*d[19]*d[25] + 2*d[15]*d[19]*d[25] + 2*d[10]*d[20]*d[25] + 2*d[12]*d[20]*d[25] + 2*d[11]*d[21]*d[25] + 2*d[15]*d[21]*d[25] - 4*d[9]*d[23]*d[25] + 2*d[10]*d[24]*d[25] + 2*d[12]*d[24]*d[25] - 2*d[9]*std::pow(d[25],2) + 8*d[13]*d[18]*d[26] - 4*d[10]*d[19]*d[26] - 4*d[12]*d[19]*d[26] - 4*d[10]*d[21]*d[26] - 4*d[12]*d[21]*d[26] + 8*d[9]*d[22]*d[26];
    coeffs[19] = -2*std::pow(d[20],2)*d[22] + 2*d[19]*d[20]*d[23] + 2*d[20]*d[21]*d[23] - 2*d[18]*std::pow(d[23],2) - 4*d[20]*d[22]*d[24] + 2*d[19]*d[23]*d[24] + 2*d[21]*d[23]*d[24] - 2*d[22]*std::pow(d[24],2) + 2*d[19]*d[20]*d[25] + 2*d[20]*d[21]*d[25] - 4*d[18]*d[23]*d[25] + 2*d[19]*d[24]*d[25] + 2*d[21]*d[24]*d[25] - 2*d[18]*std::pow(d[25],2) - 2*std::pow(d[19],2)*d[26] - 4*d[19]*d[21]*d[26] - 2*std::pow(d[21],2)*d[26] + 8*d[18]*d[22]*d[26];



	// Setup elimination template
	static const int coeffs0_ind[] = { 0, 10, 1, 0, 10, 11, 2, 1, 0, 10, 11, 12, 3, 2, 1, 11, 12, 13, 3, 2, 12, 13, 4, 0, 10, 14, 5, 4,
        14, 1, 10, 11, 0, 15, 6, 5, 4, 14, 15, 2, 11, 12, 1, 16, 7, 4, 14, 0, 10, 17, 8, 7, 17, 5, 14, 15, 1, 11, 4, 18, 6, 5, 15, 16, 3,
        12, 13, 2, 3, 13 };
	static const int coeffs1_ind[] = { 9, 19, 9, 19, 7, 17, 9, 7, 17, 4, 14, 19, 9, 19, 8, 17, 18, 5, 15, 7, 8, 7, 17, 18, 6, 15, 16, 2, 
        12, 5, 19, 8, 18, 9, 9, 19, 18, 6, 16, 8, 8, 18, 16, 3, 13, 6, 6, 16, 13, 3 };
	static const int C0_ind[] = { 0, 11, 12, 13, 16, 23, 24, 25, 26, 27, 28, 35, 36, 37, 38, 39, 40, 47, 49, 50, 51, 52, 60, 65, 67, 71, 
        72, 73, 76, 77, 78, 79, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 94, 95, 96, 101, 103, 104, 105, 107, 108, 109, 112, 113, 114, 115, 
        116, 117, 118, 119, 121, 122, 123, 124, 125, 126, 127, 130, 134, 135} ;
	static const int C1_ind[] = { 8, 9, 17, 19, 20, 21, 24, 29, 31, 32, 33, 35, 37, 40, 41, 42, 43, 44, 45, 46, 49, 50, 51, 52, 53, 54, 55, 
        56, 57, 58, 66, 68, 69, 70, 74, 75, 78, 80, 81, 82, 86, 87, 90, 92, 93, 94, 98, 99, 102, 106 };

    Eigen::MatrixXd C0 = Eigen::MatrixXd::Zero(12, 12); 
    Eigen::MatrixXd C1 = Eigen::MatrixXd::Zero(12, 9); 
	for (int i = 0; i < 70; i++) { 
        C0(C0_ind[i]) = coeffs(coeffs0_ind[i]); 
    }
	for (int i = 0; i < 50; i++) {
         C1(C1_ind[i]) = coeffs(coeffs1_ind[i]); 
    } 
    
	Eigen::MatrixXd C12 = C0.partialPivLu().solve(C1);

	// Setup action matrix
	Eigen::Matrix<double,12, 9> RR;
	RR << -C12.bottomRows(3), Eigen::Matrix<double,9,9>::Identity(9, 9);

	static const int AM_ind[] = { 8,6,0,7,1,9,10,11,2 };
	Eigen::Matrix<double, 9, 9> AM;
	for (int i = 0; i < 9; i++) {
		AM.row(i) = RR.row(AM_ind[i]);
	}

    // Eigen::Matrix<std::complex<double>, 2, 9> sols;
	sols.setZero();

	// Solve eigenvalue problem
    double p[1 + 9];
    Eigen::Matrix<double, 9, 9> AMp = AM;
    charpoly_danilevsky_piv(AMp, p);
    double roots[9];
    int nroots;
    find_real_roots_sturm(p, 9, roots, &nroots, 8, 0);
    fast_eigenvector_solver(roots, nroots, AM, sols);

	// EigenSolver<Matrix<double, 9, 9> > es(AM);
	// Eigen::ArrayXcd D = es.eigenvalues();	
	// Eigen::ArrayXXcd V = es.eigenvectors();
    // V = (V / V.row(0).array().replicate(9, 1)).eval();


    // sols.row(0) = V.row(1).array();
    // sols.row(1) = D.transpose().array();

	return nroots; //sols;
}

namespace poselib {
    int relpose_6pt_planar(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
                       std::vector<Eigen::Matrix3d> *fundamental_matrices) {

    // Compute nullspace to epipolar constraints
    Eigen::Matrix<double, 9, 6> epipolar_constraints;
    for (size_t i = 0; i < 6; ++i) {
        epipolar_constraints.col(i) << x1[i](0) * x2[i], x1[i](1) * x2[i], x1[i](2) * x2[i];
    }
    Eigen::Matrix<double, 9, 9> Q = epipolar_constraints.fullPivHouseholderQr().matrixQ();
    Eigen::Matrix<double, 9, 3> N = Q.rightCols(3);

    Eigen::VectorXd B(Eigen::Map<Eigen::VectorXd>(N.data(), N.cols() * N.rows()));

    Eigen::Matrix<std::complex<double>, 2, 9> sols;
    int n_sols = solver_6pt_planar(B, sols);
    
    fundamental_matrices->clear();
    fundamental_matrices->reserve(n_sols);

    for (int i = 0; i < n_sols; i++) {
        Eigen::Vector<double, 9> fundamental_matrix_vector = sols(0, i).real() * N.col(0) + sols(1, i).real() * N.col(1) + sols(2, i).real() * N.col(2) + N.col(3);
		fundamental_matrix_vector.normalize();
        Eigen::Matrix3d fundamental_matrix = Eigen::Map<Eigen::Matrix3d>(fundamental_matrix_vector.data());
        fundamental_matrices->push_back(fundamental_matrix);
    }

    return n_sols;
}

// int relpose_6pt_planar(const std::vector<Eigen::Vector3d> &x1, const std::vector<Eigen::Vector3d> &x2,
//                        std::vector<CameraPose> *output) {
//     std::vector<Eigen::Matrix3d> fundamental_matrices;
	
//     int n_sols = relpose_6pt_planar(x1, x2, &fundamental_matrices);

//     output->clear();
//     output->reserve(n_sols);
//     // for (int i = 0; i < n_sols; ++i) {
//     //     motion_from_essential(essential_matrices[i], x1[0], x2[0], output);
//     // }

//     return output->size();
// }

} // namespace poselib