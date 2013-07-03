struct extra_struct
{
	pointer	cor;		// Correlation matrix of the modes (for the cross term)
	pointer var_ortho;	// Variance of orthogonal modes nomalized at D/r0 = 1
	pointer cww;		// Covariance matrix of the measurements
	pointer cMat;		// Command matrix of the system
	pointer iMat;		// Interaction matrix for an "infinite" basis of modes 
	pointer modes;		// 3D matrix containing the modes used to do the iMat
	pointer cbact;		// Circular buffer of the coefficient of the modes
	pointer ipupil;		// Pupil
	pointer	imav;		// PSF to obtain
	pointer	alias		// Aliased modes in the DM space
	pointer	Dphi_ortho;	// Structure function of the orthogonal phase
	int		nact;		// Number of corrected modes
	float	normJ;		// Normalization factor for the criterion
	float	normModes;	// In case of the modes basis is not normalized to 1  i.e. (Mi . Mi) != 1
	float	teldiam;	// Telescope diameter
	float	cobs;		// Central obstruction in fraction of tel diam
	int		niter;		// Number of iteration
	float	lambdawfs;	// Wavelength of the WFS in micron
	float	lambdaim;	// Wavelength of the observation in micron
}

func correlate(f, g, dir=) 
/* DOCUMENT
 
	Compute the correlation with FFTs of 2 array.
	Auto-correlation f = g
 
 */
{
	dir = dir;
	if (is_void(dir)) dir = array(1, dimsof(f)(1));
	
	if (allof(f == g) == 1) {									// if auto-correlation
		cor		= fft(abs(fft(f, dir))^2, -dir);
	} else {
		ft_f	= fft(f, dir);
		ft_g	= fft(g, dir);
		conv	= conj(ft_f) * ft_g;
		cor		= fft(conv, -dir);								// caution with the normalization
	}
	
	return cor;
}

func calc_Uij(mode1, mode2, pupil, otf=) {						//
	/* DOCUMENT
	 
		Compute Uij modes
	 
	 */
	p	= pupil;
	
	if (is_void(otf)) otf = correlate(p, p).re;
	
	mask		= (otf > max(otf) * 1.e-7);
	pix			= where(mask);
	uij			= mode1 * 0.;
	mi			= mode1;
	mj			= mode2;
	otf			= otf(pix);
	
	tmp1		= -2. * correlate(mi * p, mj * p);
	tmp2		= 2. * correlate(mi * mj * p, p);
	
	uij(pix)	= (tmp1 + tmp2).re(pix) / otf;
	
	return uij;
}

func calc_Viif(ftMode1, ftMode2, mode1, mode2, den, conjftpup)
/* DOCUMENT
 
	Compute Vii modes fast
 
 */
{
	mask     = den > max(den) * 1.e-7;
	pix      = where(mask);
	Vii      = mode1 * 0.;
	Vii(pix) = ((-2. * fft((ftMode1 * conj(ftMode2)).re, -1) + \
				 2. * fft((fft(mode1 * mode2, 1) * conjftpup).re, -1)).re)(pix) / (den(pix));
	
	return Vii;
}

func calc_dphif(phase, pup, den, conjftpup)
/* DOCUMENT 

	Compute the structure function of the phase
 
*/ 
{
	phase	= phase(:pupd, :pupd);
	pup		= pup(:pupd, :pupd);
	
	npix	= dimsof(phase)(2);
	mi		= p = dphi = array(float, [2, 2 * npix, 2 * npix]);
	
	mi(1:npix, 1:npix)	= phase;
	p(1:npix, 1:npix)	= pup;
	
	mask	= den > max(den) * 1.e-7;
	pix		= where(mask);
	
    dphi(pix)	= fft(fft(mi^2 * p, 1) * conjftpup + fft(p, 1) * conj(fft(mi^2 * p, 1)) - \
					2. * fft(mi * p, 1) * conj(fft(mi * p, 1)), -1).re(pix) / den(pix);

	return dphi;
}



func get_param(cov, d_r0, extra)
/* DOCUMENT
 
	Get the parameters vector for minimization using J1J2_diag()
	The function does not take into account the off-diagonal covariance
 
*/

{
	par	= [];
	grow, par, diag(cov)(:extra.nact), d_r0;
	
	return par;
}

func map_indices(matind, d_r0, extra)
/* DOCUMENT
 
	Tool to find the position of a parameter in the covariance matrix
 
*/
{
	ind	= [];
	dim	= extra.nact;
	
	grow, ind, diag(matind)(:dim);
	
	for (j=1 ; j<=dim ; j++) {				// X = [var, cor, r0]
		for (i=1+j ; i<=dim ; ++i) {
			grow, ind, matind(i, j);
		}
	}
	
	if (is_void(d_r0)) {
		write, "/!\\ no D/r0 spefified, set to 30";
		d_r0 = 30.;
	}
	grow, ind, d_r0;
	
	return ind;
}

func cov2param(cov, d_r0, extra)
/* DOCUMENT
 
	Build the parameters vector from the covariance and D_r0.
	Parameters are sorted such as : X = [variances, covariances, D_r0]
 
 */
{
	X	= [];
	dim	= extra.nact;
	
	grow, X, diag(cov)(:dim);
	
	cor	= cov2cor(cov, var2);
	
	for (j=1 ; j<=dim ; j++) {
		for (i=1+j ; i<=dim ; ++i) {
			grow, X, cor(i, j);
		}
	}
	
	if (is_void(d_r0)) {
		write, "/!\\ no D/r0 spefified, set to 30";
		d_r0 = 30.;
	}
	grow, X, d_r0;
	
	return double(X);
}

func param2cov(X, extra)
/* DOCUMENT
 
	Compute the covariance matrix from the parameters.
 
 */
{
	var		= [];
	dim		= extra.nact;
	cov		= *extra.cor;							// getting cross-ortho correlation terms
	cormir	= param2cor_mir(X, extra);				// building mirror covariance
	
	cov(:dim, :dim)	= cormir;						// adding mirror correlation matrix to correlation matrix
	
	d_r0	= X(0);
	var		= X(:extra.nact);
	
	grow, var, (*extra.var_ortho) * d_r0^(5. / 3);
	
	nmodes	= numberof(var);
	
	for (i=1 ; i<=nmodes ; ++i) {
		cov(i, ) *= sqrt(var(i));
		cov(, i) *= sqrt(var(i));
	}
	
	return cov;
}

func param2cor_mir(X, extra)
/* DOCUMENT
 
	Compute the correlation matrix using the parameters in the mirror space.
 
 */
{
	dim		= extra.nact;
	
	cor		= xy = reform(indgen(1:dim*dim), dim, dim);			// mapping
	cor		*= 0.;
	xy		= long(map_indices(xy, -1000, extra));				// fake r0
	indcor	= xy(dim+1:numberof(xy)-1);
	inddiag	= xy(:dim);
	
	
	cor(*)(indcor)	= X(dim+1:numberof(X)-1);
	cor				+= transpose(cor);
	cor(*)(inddiag)	= unit(dim)(inddiag);						// diagonal
	
	return cor;
}


func J1J2_diag_nor0(X, &gx, extra)
/* DOCUMENT
 
	Criterion to ONLY fit the VARIANCE.
 
 */
{
	
	N		= niter;
	Nw		= nmes;
	dim		= extra.nact;
	dimDt	= dimsof(Dtot)(3);
	cww		= *extra.cww;
	grow, X, 35.;
	Ce		= cor2cov(extra, X);
	X=X(:dim);
	U		= Dtot(, +) * (Ce(, +) * Dtot(, +))(+, ) + cnn;
	Um1		= SVsolve(U(+, ) * U(+, ), transpose(U)); // inversion for highly singular matrices
	
	J1		= 0.5 * N * ln_det(U);
	J2		= 0.5 * N * trace(Um1( , +) * cww(+, ));
	
	gx		= grad_J1J2_diag_nor0(X, extra);
	gx		/= extra.normJ;
	
	J		= J1 + J2;
	J		/= extra.normJ;
	
	return J;
}

func grad_J1J2_diag_nor0(X, extra)
/* DOCUMENT
 
	Gradient computation for the J1J2_diag_nor0 function.
 
 */
{
	gx		= [];
	
	N		= niter;
	dim		= extra.nact;
	var		= X(:dim);
	dimDt	= dimsof(Dtot)(3);
	cww		= *extra.cww;
	grow, X, 35.;
	Ce		= cor2cov(extra, X);
	X=X(:dim);
	U		= (Dtot(, +) * Ce(+, ))(, +) * Dtot(, +) ;
	Um1		= LUsolve(U);
	Um1t	= transpose(Um1);
	
	g		= var * 0.;
	fact_wk = transpose(-(Um1t( , +) * \
						  (N * cww)(+, ))( , +) * Um1t(+, ));
	
	for (i=1 ; i<=dim ; ++i) {
		dc			= Dii(Ce, i);
		D_dc_Dt		= Dtot(, +) * (dc(, +) * Dtot(, +))(+, );
		dJ1			= 0.5 * N * trace(Um1(, +) * D_dc_Dt(+, ));
		dJ2			= 0.5 * trace(fact_wk(, +) * D_dc_Dt(+, ));
		g(i)		= dJ1 + dJ2;
	}
	
	gx	= g;
	
	return gx;
}

func J1J2_diag(X, &gx, extra)
/* DOCUMENT
 
	Criterion to fit the r0 and the VARIANCE only.
 
 */
{
	
	N		= extra.niter;
	Dr0		= X(0);
	dim		= extra.nact;
	Dtot	= *extra.iMat;
	dimDt	= dimsof(Dtot)(3);
	cww		= *extra.cww;
	Ce		= cor2cov(extra, X);
	U		= Dtot(, +) * (Ce(, +) * Dtot(, +))(+, ) + cnn;
	Um1		= SVsolve(U(+, ) * U(+, ), transpose(U)); // inversion for highly singular matrices
	
	J1		= 0.5 * N * ln_det(U);
	J2		= 0.5 * N * trace(Um1( , +) * cww(+, ));
	
	gx		= grad_J1J2_diag(X, extra);
	gx		/= extra.normJ;
	
	J		= J1 + J2;
	J		/= extra.normJ;
	
	return J;
}

func grad_J1J2_diag(X, extra)
/* DOCUMENT
 
	Gradient computation for the J1J2_diag() function.
 
 */
{
	gx		= [];
	
	N		= extra.niter;
	Dr0		= X(0);
	dim		= extra.nact;
	var		= X(:dim);
	Dtot	= *extra.iMat;
	dimDt	= dimsof(Dtot)(3);
	cww		= *extra.cww;
	Ce		= cor2cov(extra, X);
	U		= (Dtot(, +) * Ce(+, ))(, +) * Dtot(, +) ;
	Um1		= LUsolve(U);
	Um1t	= transpose(Um1);
	
	g		= var * 0.;
	fact_wk = transpose(-(Um1t( , +) * \
						  (N * cww)(+, ))( , +) * Um1t(+, ));
	
	for (i=1 ; i<=dim ; ++i) {
		dc			= Dii(Ce, i);
		D_dc_Dt		= Dtot(, +) * (dc(, +) * Dtot(, +))(+, );
		dJ1			= 0.5 * N * trace(Um1(, +) * D_dc_Dt(+, ));
		dJ2			= 0.5 * trace(fact_wk(, +) * D_dc_Dt(+, ));
		g(i)		= dJ1 + dJ2;
	}
	
	dr		= Rij(Ce, Dr0, extra);
	dr0dt	= (Dtot(, +) * dr(+, ))(, +) * Dtot(, +);
	dJ1r0	= 0.5 * N * trace(Um1(, +) * dr0dt(+, ));
	dJ2r0	= 0.5 * trace(fact_wk(, +) * dr0dt(+, ));
	
	gr0		= dJ1r0 + dJ2r0;
		
	grow, gx, g, gr0;
	
	return gx;
}




func J1J2(X, &gx, extra)
/* DOCUMENT
 
	Criterion to minimize.
 
 */
{
	
	N			= extra.niter;
	Dtot		= *extra.iMat;
	Dr0			= X(0);
	dim			= extra.nact;
	dimDt		= dimsof(Dtot)(3);
	cww			= *extra.cww
		
	Ce			= param2cov(X, extra);
	U			= Dtot(, +) * (Ce(, +) * Dtot(, +))(+, ) + cnn;
	Um1			= SVsolve(U(+, ) * U(+, ), transpose(U)); // inversion for highly singular matrices
	
	J1		= 0.5 * N * ln_det(U);
	J2		= 0.5 * N * trace(Um1( , +) * cww(+, ));
	
	gx		= grad_J1J2(X, extra);
	gx		/= extra.normJ;
		
	J		= J1 + J2;
	J		/= extra.normJ;
	
	return J;
}


func grad_J1J2(X, extra)
/* DOCUMENT
 
	Gradient computation for the J1J2() function.
 
 */
{
	gx		= [];
	
	N			= extra.niter;
	Dtot		= *extra.iMat;
	Dr0			= X(0);
	dim			= extra.nact;
	var			= X(:dim);
	dimDt		= dimsof(Dtot)(3);
	cww			= *extra.cww;
		
	Ce			= param2cov(X, extra);
	U			= (Dtot(, +) * Ce(+, ))(, +) * Dtot(, +) ;
	Um1			= SVsolve(U(+, ) * U(+, ), transpose(U));
	Um1t		= transpose(Um1);
	
	g_diag		= var * 0.;
	g_nodiag	= array(float, numberof(X)-dim-1);						// number of element - diag - d_r0
	
	fact_wk = transpose(-(Um1t( , +) * \
						  (N * cww)(+, ))( , +) * Um1t(+, ));
	
	for (i=1 ; i<=dim ; ++i) {
		dv			= Dii(Ce, i);										// derivative of variance
		D_dv_Dt		= Dtot(, +) * (dv(, +) * Dtot(, +))(+, );
		dJ1			= 0.5 * N * trace(Um1(, +) * D_dv_Dt(+, ));
		dJ2			= 0.5 * trace(fact_wk(, +) * D_dv_Dt(+, ));
		g_diag(i)	= dJ1 + dJ2;
	}
	
	k	= 0;
	for (i=1 ; i<dim ; ++i) {
		for (j=i ; j<=dim ; ++j) {
			if (i != j) {
				++k;
				dc			= Dij(Ce, i, j);							// derivative of covariance
				D_dc_Dt		= Dtot(, +) * (dc(, +) * Dtot(, +))(+, );
				dJ1			= 0.5 * N * trace(Um1(, +) * D_dc_Dt(+, ));
				dJ2			= 0.5 * trace(fact_wk(, +) * D_dc_Dt(+, ));
				g_nodiag(k)	= dJ1 + dJ2;								// element and its transpose vary
			}
		}
	}
	
	dr		= Rij(Ce, Dr0, extra);										// derivative wrt r0
	dr0dt	= (Dtot(, +) * dr(+, ))(, +) * Dtot(, +);
	dJ1r0	= 0.5 * N * trace(Um1(, +) * dr0dt(+, ));
	dJ2r0	= 0.5 * trace(fact_wk(, +) * dr0dt(+, ));
	
	gr0		= dJ1r0 + dJ2r0;
	
	grow, gx, g_diag, g_nodiag, gr0;
	
	return gx;
}

func Dii(Ce, i)
/* DOCUMENT
 
	Derivative of the covariance dC/dC_ii  i.e. wrt the variance 
 
 */
{
	cor		= cov2cor(Ce, var);
	D		= Ce * 0;
	D(, i)	= D(i, ) = 0.5 * Ce(i, ) / var(i);  
	D(i, i)	= 1;
	
	return D;
}

func Dij(Ce, i, j)
/* DOCUMENT
	
	Derivative of covariance dC/d(rho_ij)  i.e. wrt the correlation coefficent. 
 
 */
{
	cor		= cov2cor(Ce, var);
	m		= numberof(var);
	D		= Ce * 0;
	D(i, j)	= D(j, i) = sqrt(var(i)) * sqrt(var(j));
	
	return D;
}

func Rij(Ce, Dr0, extra)
/* DOCUMENT
 
	Derivative of covariance dC/d(D/r0)  i.e. wrt D/r0
 
 */
{
	cor		= cov2cor(Ce, var);
	
	si		= var(:extra.nact);
	sj		= *extra.var_ortho;
	
	grow, normvar, si, sj, 1;				// for D/r0 = 1 here because of cor2cov that take into account the r0
	
	nCe		= cor2cov(extra, normvar);		// normalized at D/r0 = 1
	
	nCe(:extra.nact, :extra.nact) *= 0.;	// mirror space not affected by r0
	
	dr0		= (5. / 6) * Dr0^(-1. / 6) * nCe;
	
	dr0(extra.nact+1:, extra.nact+1:)	= diag((5. / 3) * Dr0^(2. / 3) * sj); // on the diagonale
	
	return dr0;
}


func ln_det(a)
/* DOCUMENT
 
	Returns the absolute value of the logarithm determinant
	of matrix by using singular value decomposition
 
 */
{
	u		= SVdec(a);											//	Decomposition of a
	N		= numberof(u);										//	Number of singular values
	K		= 1. / min(u);										//	"Condition number"
	u		*= K;												//	Multiplying to avoid very small values
	ln_DET	= noneof(u == 0) ? log(abs(u))(sum) : 0.00;			//	product of SV -> determinant
	ln_DET	= ln_DET - N * log(K);								//	Correcting the multiplication
	
	return ln_DET;
}



func cor2cov(extra, X)
/* DOCUMENT
 
	Compute the covariance matrix given the correlation matrix and the variances
 
 */
{
	var		= [];
	cov		= *extra.cor;
	Dr0		= X(0);
	var		= X(:extra.nact);
	
	grow, var, (*extra.var_ortho) * Dr0^(5. / 3);
	
	nmodes	= numberof(var);
	
	for (i=1 ; i<=nmodes ; ++i) {
		cov(i, ) *= sqrt(var(i));
		cov(, i) *= sqrt(var(i));
	}
	
	return cov;
}


func cov2cor(cov, &var)
/* DOCUMENT
 
	Compute the correlation matrix given the covariance.
	The function returs the correlation matrix and the variance vector as a pointer.
 
 */
{
	cor		= cov;
	var		= diag(cov);
	
	nmodes	= numberof(var);
	
	for (i=1 ; i<=nmodes ; ++i) {
		cor(i, )	/= sqrt(var(i));
		cor(, i)	/= sqrt(var(i));
	}
	return cor;
}



func calc_cor_term(n, g, ipupil, eps_para, cMat, Dpinf, eps_ortho, phi_para)
/* DOCUMENT
 
	Iterative method
	to develop
 
 */
{
	niter	= dimsof(eps_ortho)(3);
	pix		= where(ipupil);
	
	dimp	= dimsof(cMat)(2);
	dimo	= dimsof(Dpinf)(3);
	
	c			= array(float, [2, dimp + dimo, niter]);
	
	write, ""; write, "Computing aliased commands"; write, "";
	com_ortho	= (cMat(, +) * (Dpinf(, +) * eps_ortho(+, ))(+, ));
	
	write, "Computing terms";
	t1	= t2 = t3 = array(float, [2, dimp, niter]);
	
	for (i=n+1 ; i<=niter ; ++i) {
		for (m=1 ; m<=n ; ++m) {
			t2(, i)	+= (1 - g)^(m - 1) * com_ortho(, i - m);
			t3(, i)	+= (1 - g)^(m - 1) * (phi_para(, i - m + 1) - phi_para(, i - m));
		}
	}
	t2	*= -g;
	
	
	write, "Computing cross terms";
	
	c(:dimp, )		= t2 + t3;
	c(dimp+1:, )	= eps_ortho;
	
	C_cross			= c(, n+1:)(, +) * c(, n+1:)(, +) / (niter - n);
	
	cor				= cov2cor(C_cross, var2);
	
	return cor;
}


func hist(t,over=,color=, pasZero=, marks=)
/* DOCUMENT 
 
	hist(tab,over=,color=, marks=, pasZero=)
	
	Displays the histogram of array 'tab'.
	See also histo, or histogramme.
 
	Red  = median value
	Blue = average value
 */
{
	local t;
	
	if (pasZero) {
		nn = where(t!=0);
		if (is_array(nn)) t = t(nn);
	}
	
	limits,square=0;
	if( is_void(over) ) fma;
	logxy,0,0;
	t = double(t);
	if (min(t)!=max(t)) {
		y = histogram(int(1.5+254.*(t-min(t))/(max(t)-min(t))));
		x = span(min(t),max(t),dimsof(y)(2));
		tab = where(y>(max(y)/2));
		plh, y, x, color=color, marks=marks;
	}
	
	stdev	= sqrt(((y - y(avg))^2)(avg)); // sigma = sqrt( variance ) 
	ind		= where(abs((y - stdev)) == min(abs(y - stdev)));
	sigma	= abs(x(ind))(1);
	
	w		= 1.;
	res		= moffat1d_fit(y, x, w);
	plg, moffat1d(x, res), x, marks=0, color="magenta";
	limits, -5 * sigma, 5 * sigma;
	write, "FWHM               : ", 2.*res(3) * sqrt(0.5^(-1/res(4)) - 1.);
	write, "Standard deviation : ", sigma;
}


func matCovAlias(nssp, obs, corr=, xy= )
/* DOCUMENT
 
	See le Rico parce que c'est chaud
 
 */
{
	local corr, xy;
	if( corr==[] )
		corr=-0.5;
	x = span(-1,1,nssp+1)(zcen)(,-:1:nssp);
	x = array(0., nssp+2, nssp+2);
	x(2:-1,2:-1) = span(-1,1,nssp+1)(zcen)(,-:1:nssp);
	y = transpose(x);
	r=sqrt(x^2+y^2);
	pup = r<0.95 & r>obs;
	nn = where(pup);
	N = numberof(nn);
	cov = array(0.,2*N,2*N);
	tmp = array(0., dimsof(x));
	for(i=1; i<=N; i++) {
		tmp *= 0;
		tmp(nn(i)) = 1;
		tmp(nn(i)+1) = tmp(nn(i)-1) = corr;
		cov(1:N,i) = tmp(nn);
		tmp *= 0;
		tmp(nn(i)) = 1;
		tmp(nn(i)+nssp+2) = tmp(nn(i)-nssp-2) = corr;
		cov(N+1:,i+N) = tmp(nn);
	}
	if( xy==1 )
		cov = roll(cov);
	return cov;
}

func circavg(a,center=,middle=)
/* DOCUMENT circavg
 *  
 * average=circavg(array[,center=,middle=])
 *
 * This routine returns the circular mean of an array. The routine dist is
 * used to compute an array of coordinate of points to be averaged.
 *
 * KEYWORDS :
 * a      (input) : The array of which we want the circular mean. It can be
 *                  a long, a float or a double. Complex are not supported
 * center (input) : An array of the form [x,y] which give the position of
 *                  the origin for the circular mean calculation
 * middle (input) : A flag to indicate that the origin for the circular mean
 *                  calculation is the middle of the array a
 *
 * SEE ALSO: dist, circavg2
 */ 
{
	s=dimsof(a);
	
	if (!is_set(middle)) middle=0;
	if (s(1) != 2) write,"error - invalid dimensions";
	if (s(3) != s(2)) write,"error - invalid dimensions";
	
	dim=s(2);
	
	if (center!=[]) {
		s=dimsof(center);
		if ((s(1) != 1) | (s(1) != 2)) \
			write,"error - center has invalid dimensions";
		
		center=long(center);
		
		if (middle) {
			center=long([0,0]);
			write,"error - center and middle are not compatible keywords";
		}
	} else { 
		if (middle) center=long([0,0]);
		else center=[dim/2,dim/2];
	}
	
	r=long(roll(long(dist(dim)+.5)+1,[center(1),center(2)]));
	j=long(max(r));
	n=array(long,j);
	sx=array(double,j);
	dim2=long(dim)*long(dim);
	
	for (i=1;i<=dim2;i++) {
		j=r(i);
		sx(j)=sx(j)+a(i);
		n(j)=n(j)+1;
	}
	
	return sx/n;
}
