require, "optimpack-mod.i"
require, "yao.i"
require, "yao_funcs.i"

func test(void)

{
	extern niter, nmes, nact, nmode, cww, cee, cnn, cwwnn, extra, D, Dtot, cbmes, act_all;
	
	aoread, "canary-TS-zernike.par";
	wfs(1).shmethod = 1;
	aoinit, clean=1, forcemat=1;
	
	cbmes		= double(fits_read("cbmes-150Hz_gain0.500000only200modes.fits"))(145:,);
	cberr		= double(fits_read("cberr-150Hz_gain0.500000only200modes.fits"));
	cbcom		= double(fits_read("cbcom-150Hz_gain0.500000only200modes.fits"));
	cbphase		= double(fits_read("cbphase-150Hz_gain0.500000only200modes.fits"));
	smesal		= double(fits_read("cbmes-150Hz_gain0.500000only200modes.fits"))(:72,);
	smesnn		= cbmes; // pour geometrique //
	//smesnn	= fits_read("smes_nonoise-300Hz_gain0.300000only200modes.fits")(, 2:);
	act_all		= double(fits_read("act_all-150Hz_gain0.500000only200modes.fits")(, 2:));
	act_ortho	= double(fits_read("act_ortho-150Hz_gain0.500000only200modes.fits")(, 2:));
	//kl			= fits_read("new_kl.fits");
	D			= double(fits_read("iMat_geom.fits"));
	Dtot		= double(fits_read("iMat_tot_geom.fits"));
	Dplus		= double(fits_read("cMat_geom.fits"));
	mc_cov		= fits_read("cov_mc_150.fits");
	
	// Dtot pas bonne unit phase micron || inutil en sh geometrique//
	//D		*= 0.912871;
	//Dtot	*= 0.912871;
	
	niter	= dimsof(act_all)(3);
	nmes	= dimsof(Dtot)(2);
	nact	= dimsof(D)(3);
	nmode	= dimsof(Dtot)(3);
	g		= loop.gain;
	
	cee		= act_all(, +) * act_all(, +) / niter;
	
	n		= cbmes - smesnn;
	cnn		= diag(diag(n(, +) * n(, +) / niter))(:nmes, :nmes);
	
	cww		= cbmes(, +) * cbmes(, +) / niter;
	cwwnn	= smesnn(, +) * smesnn(, +) / niter;
	
	crr		= smesal(, +) * smesal(, +) / niter;
	cer		= smesnn(, +) * smesal(, +) / niter;
	cre		= smesal(, +) * smesnn(, +) / niter;
	cerec	= (Dplus(, +) * (cwwnn+crr-cer-cre)(+, ))(, +) * Dplus(, +); // eq3.11 these damien
	
	/*
	epsmes	= Dplus(, +) * cbmes(+, );				// modes via mesures bruitees
	epsn	= Dplus(, +) * n(+, );					// modes via bruit
	epsr	= Dplus(, +) * smesal(+, );				// modes via aliasing
	epsmnn	= Dplus(, +) * smesnn(+, );				// modes via mesures non bruitees
	
	epsp	= epsmnn - epsr;						//eq3.9 these damien
	
	cerec	= epsp(, +) * epsp(, +) / niter;
	*/
	
	extra			= extra_struct();
	extra.cor		= &(cov2cor(cee,var));
	extra.var_ortho	= &(diag(cee)(nact + 1:) / 35.^(5. / 3)); // at D/r0 =1
	extra.cww		= &(cww);
	extra.nact		= nact;
	extra.normJ		= 1.e7;
	extra.teldiam	= 4.2;
	//extra(:nact, :nact)	-= diag(diag(cee)(:nact));
	//tmpt	= extra(31:,31:)-diag(diag(extra(31:,31:)));
	//extra(31:,31:)	-= tmpt;
	
	
	window,1,dpi=140;
	
	i = 20;
	
	plg, act_all(,i),marks=0, width=3, color="red";
	tmp = cbphase(,i);
	tmp(:nact) += cbcom(2*nact+1:,i -1);
	plg, tmp,marks=0;
	
	typeReturn;
	
	tmp2 = cbphase(,i);
	tmp2(:nact) += (cbcom(2*nact+1:,i - 2) - g * cberr(2 * nact + 1:,i - 1));
	plg, tmp2, marks=0, color="green";
	
	typeReturn;
	
	fma;
	plg,cberr(2*nact+1:,i), marks=0, width=3;
	plg,cMat(:30,:72)(,+)*(Dtot(,:30)(,+)*act_all(:30,i)(+))(+,) \
		+cMat(:30,:72)(,+)*(Dtot(,31:)(,+)*act_all(31:,i)(+))(+,), color="green", marks=0;
	
	typeReturn;
	
	/*fma;
	plg,act_all(:nact, i),width=3,marks=0;
	plg, (1-g)*act_all(:30, i-1) + (cbphase(:30,i) - cbphase(:30, i-1)) \
		- g * (cMat(:30,:72)(,+)*(Dtot(,31:)(,+)*act_all(31:,i-1)(+))(+)), color="green", marks=0;
	
	pltitle, "recurrence n=1";
	
	typeReturn;
	
	fma;
	plg,act_all(:nact, i),width=3,marks=0;
	epsm1	= (1-g)*act_all(:30, i-2) + (cbphase(:30,i-1) - cbphase(:30, i-2)) \
		- g * (cMat(:30,:72)(,+)*(Dtot(,31:)(,+)*act_all(31:,i-2)(+))(+)); // eps para à i-1 = f(i-2)
	
	plg, (1-g)*(epsm1) + (cbphase(:30,i) - cbphase(:30, i-1)) \
	- g * (cMat(:30,:72)(,+)*(Dtot(,31:)(,+)*act_all(31:,i-1)(+))(+)), color="green", marks=0;
	
	pltitle, "recurrence n=2";
	typeReturn;
	
	
	fma;
	plg,act_all(:nact, i),width=3,marks=0;
	epsm2	= (1-g)*act_all(:30, i-3) + (cbphase(:30,i-2) - cbphase(:30, i-3)) \
	- g * (cMat(:30,:72)(,+)*(Dtot(,31:)(,+)*act_all(31:,i-3)(+))(+)); // eps para à i-1 = f(i-2), etc.
	
	epsm1	= (1-g)*epsm2 + (cbphase(:30,i-1) - cbphase(:30, i-2)) \
	- g * (cMat(:30,:72)(,+)*(Dtot(,31:)(,+)*act_all(31:,i-2)(+))(+));
	
	epsi	= (1-g)*epsm1 + (cbphase(:30,i) - cbphase(:30, i-1)) \
	- g * (cMat(:30,:72)(,+)*(Dtot(,31:)(,+)*act_all(31:,i-1)(+))(+));
	
	plg, epsi, color="green", marks=0;
	
	pltitle, "recurrence n=3";
	typeReturn;*/
	
	
	for (n=1 ; n<=100 ; n+=50){
		
		fma;
		//write, "Computing terms";
		t1	= t2 = t3 = array(float, nact);
	
		com	= (cMat(:nact,:nmes)(,+)*(Dtot(,nact+1:)(,+)*act_all(nact+1:,)(+, ))(+, ));
	
	
		t1		= (1 - g)^n * act_all(:nact, i-n);
	
		for (m=1 ; m<=n ; ++m) {
			t2		+= (1 - g)^(m - 1) * com(, i - m);
			t3		+= (1 - g)^(m - 1) * (cbphase(:nact, i-m+1) - cbphase(:nact, i - m));
		}
		t2		*= -g;
	
		plg, act_all(:nact, i), marks=0, width=3;
		pltitle, "n = " + strtrim(swrite(n));
		plg, t1 + t2 + t3, marks=0, color="green";
		typeReturn;
	}
	n	= [];
	error;
}






/*
func epara_recur(iter, n, g, Dtot, C, phip, epara, eortho)  // FONCTION RECURSIVE A DEBUG
{
	
	flag	= (n <= 0) ? 1 : 0;
	flag	= &flag;
	
	i		= iter;
	nact	= dimsof(phip)(2);
	Dinf	= Dtot(, nact + 1:);									// iMat pour modes > 30
	
	como	= C(,+) * (Dinf(, +) * eortho(+, i - n))(+);			// com ortho
	
	dphi	= phip(, i) - phip(, i - n);							// delta phi
	
	ep		= epara_recur(iter, n - 1, g, Dtot, C, phip, epara, eortho);
	ep		*= 1 - g;
	
	ep		+= dphi - g * como;
	
	if (*flag) return ep;
}
*/



func get_vec_diag(Ce)
/* DOCUMENT
 
 
 */
{
	X		= [];
	dimD	= dimsof(D)(3);
	
	X		= diag(Ce(:dimD, :dimD));
	
	return double(X);
}

func get_mat_diag(vec, extra)
/* DOCUMENT
 
 
 */
{
	if (is_void(extra)) extra = 0.;

	dimD				= dimsof(D)(3);
	dimDt				= dimsof(Dtot)(3);
	Ce					= array(float, [2, dimDt, dimDt]);
	
	Ce(:dimD, :dimD)	= diag(vec(:dimD));
	Ce					+= extra;
	
	return Ce;
}

func J1_diag(X, &gx, extra)
/* DOCUMENT
 */
{
	if (is_void(extra)) extra = 0.;
	
	N		= niter;
	Nw		= nmes;
	dim		= dimsof(X);
	dimD	= dimsof(D)(3);
	dimDt	= dimsof(Dtot)(3);
	Ce		= get_mat_diag(X, extra);
	
	U		= Dtot(, +) * (Ce(, +) * Dtot(, +))(+, ) + cnn;
	
	J1		= 0.5 * N * ln_det(U);
	
	return J1;
}


func J2_diag(X, &gx, extra)
{

	local J2;
	
	if (is_void(extra)) extra = 0.;
	
	N		= niter;
	Nw		= nmes;
	dim		= dimsof(X);
	dimDt	= dimsof(Dtot)(3);
	Ce		= get_mat_diag(X, extra);
	U		= Dtot(, +) * (Ce(, +) * Dtot(, +))(+, ) + cnn;
	Um1		= SVsolve(U(+, ) * U(+, ), transpose(U));
	//Um1		= LUsolve(U);
	
	J2		= 0.5 * N * trace(Um1( , +) * cww(+, ));
	
	/*J2		= 0.;
	for (i=1 ; i<=N ; i+=2) {
	
		J2	+= cbmes(, i)(+) * ( Um1(,+) * cbmes(, i)(+))(+);
		
	}
	J2		*= 0.5;
	//gx		= grad_J2_diag(X); 
	*/
	return J2;
}


//////////////////////////////////////////
//SVsolve
////////////////////////////////////////
func get_param(cov, d_r0, extra)
{
	par	= [];
	grow, par, diag(cov)(:extra.nact), d_r0;
	
	return par;
}

func map_indices(matind, d_r0, extra)
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
{
	X	= [];
	dim	= extra.nact;
	
	grow, X, diag(cov)(:dim);
	
	cor	= cov2cor(cov, var2);
	
	for (j=1 ; j<=dim ; j++) {				// X = [var, cor, r0]
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
{
	dim		= extra.nact;
	
	cor		= xy	= reform(indgen(1:dim*dim), dim, dim);		// mapping
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
	
	//if (is_void(gx)) write, "/!\\ gradient computaion disabled";
	
	J		= J1 + J2;
	J		/= extra.normJ;
	
	return J;
}

func grad_J1J2_diag_nor0(X, extra)
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
	//gx = get_vec_diag(g);
	
	gx	= g;
	
	return gx;
}


func J1J2_diag(X, &gx, extra)
{

	N		= niter;
	Nw		= nmes;
	Dr0		= X(0);
	dim		= extra.nact;
	dimDt	= dimsof(Dtot)(3);
	cww		= *extra.cww;
	Ce		= cor2cov(extra, X);
	U		= Dtot(, +) * (Ce(, +) * Dtot(, +))(+, ) + cnn;
	Um1		= SVsolve(U(+, ) * U(+, ), transpose(U)); // inversion for highly singular matrices
	
	J1		= 0.5 * N * ln_det(U);
	J2		= 0.5 * N * trace(Um1( , +) * cww(+, ));

	gx		= grad_J1J2_diag(X, extra);
	gx		/= extra.normJ;
	
	//if (is_void(gx)) write, "/!\\ gradient computaion disabled";
	
	J		= J1 + J2;
	J		/= extra.normJ;
	
	return J;
}

func grad_J1J2_diag(X, extra)
{
	gx		= [];
	
	N		= niter;
	Dr0		= X(0);
	dim		= extra.nact;
	var		= X(:dim);
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
		//D_Jij		= Dtot(, i)(-, ) * Dtot(, j);  // transpose(D Jij Dt) (si un seul element varie
		//dJ2			= 0.5 * trace(fact_wk(, +) * D_Jij(+, ));
		//dJ1			= 0.5 * N * trace(Um1(, +) * D_Jij(+, ));
		dc			= Dii(Ce, i);
		D_dc_Dt		= Dtot(, +) * (dc(, +) * Dtot(, +))(+, );
		dJ1			= 0.5 * N * trace(Um1(, +) * D_dc_Dt(+, ));
		dJ2			= 0.5 * trace(fact_wk(, +) * D_dc_Dt(+, ));
		//dJ			= (i==j) ? dJ1 + dJ2 : 2. * (dJ1 + dJ2);
		g(i)		= dJ1 + dJ2;
	}
	
	//fact_r0	= -(5. / 3) * extra.teldiam^(5. / 3) * CK * r0^(-8. / 3);	// prendre en compte Cross terms... done
	dr		= Rij(Ce, Dr0, extra);
	dr0dt	= (Dtot(, +) * dr(+, ))(, +) * Dtot(, +);
	dJ1r0	= 0.5 * N * trace(Um1(, +) * dr0dt(+, ));
	dJ2r0	= 0.5 * trace(fact_wk(, +) * dr0dt(+, ));
	
	gr0		= dJ1r0 + dJ2r0;
	//gx = get_vec_diag(g);
	
	grow, gx, g, gr0;
	
	return gx;
}




func J1J2(X, &gx, extra)
{
	
	N			= niter;
	Nw			= nmes;
	Dr0			= X(0);
	dim			= extra.nact;
	dimDt		= dimsof(Dtot)(3);
	cww			= *extra.cww
	
	//(*extra.cor)(:dim, :dim)	= cov2cor(param2cov_mir(X, extra), v);
	
	Ce			= param2cov(X, extra);
	U			= Dtot(, +) * (Ce(, +) * Dtot(, +))(+, ) + cnn;
	Um1			= SVsolve(U(+, ) * U(+, ), transpose(U)); // inversion for highly singular matrices
	
	J1		= 0.5 * N * ln_det(U);
	J2		= 0.5 * N * trace(Um1( , +) * cww(+, ));
	
	gx		= grad_J1J2(X, extra);
	gx		/= extra.normJ;
	
	//if (is_void(gx)) write, "/!\\ gradient computaion disabled";
	
	J		= J1 + J2;
	J		/= extra.normJ;
	
	return J;
}


func grad_J1J2(X, extra)
{
	gx		= [];
	
	N			= niter;
	Dr0			= X(0);
	dim			= extra.nact;
	var			= X(:dim);
	dimDt		= dimsof(Dtot)(3);
	cww			= *extra.cww;
	
	//(*extra.cor)(:dim, :dim)	= cov2cor(param2cov_mir(X, extra), v);
	
	Ce			= param2cov(X, extra);
	U			= (Dtot(, +) * Ce(+, ))(, +) * Dtot(, +) ;
	Um1			= SVsolve(U(+, ) * U(+, ), transpose(U));
	Um1t		= transpose(Um1);
	
	g_diag		= var * 0.;
	g_nodiag	= array(float, numberof(X)-dim-1); // number of element - diag - d_r0
	
	fact_wk = transpose(-(Um1t( , +) * \
						  (N * cww)(+, ))( , +) * Um1t(+, ));
	
	for (i=1 ; i<=dim ; ++i) {
		dv			= Dii(Ce, i);									// derivative of variance
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
				//fma;pli,dc(:30,:30);typeReturn;
				D_dc_Dt		= Dtot(, +) * (dc(, +) * Dtot(, +))(+, );
				dJ1			= 0.5 * N * trace(Um1(, +) * D_dc_Dt(+, ));
				dJ2			= 0.5 * trace(fact_wk(, +) * D_dc_Dt(+, ));
				g_nodiag(k)	= dJ1 + dJ2;								// element and its transpose vary
			}
		}
	}
	
	dr		= Rij(Ce, Dr0, extra);									// derivative wrt r0
	dr0dt	= (Dtot(, +) * dr(+, ))(, +) * Dtot(, +);
	dJ1r0	= 0.5 * N * trace(Um1(, +) * dr0dt(+, ));
	dJ2r0	= 0.5 * trace(fact_wk(, +) * dr0dt(+, ));
	
	gr0		= dJ1r0 + dJ2r0;

	
	//gx = get_vec_diag(g);
	
	grow, gx, g_diag, g_nodiag, gr0;
	
	//fma;plg,gx;redraw;
	return gx;
}

func Dii(Ce, i)
/* DOCUMENT
 derivative of covariance dC/dC_ii 
 
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
 derivative of covariance dC/d rho_ij 
 
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
 derivative of covariance dC/dr0 
 
 */
{
	cor		= cov2cor(Ce, var);
	
	si		= var(:extra.nact);
	sj		= *extra.var_ortho;
	
	grow, normvar, si, sj, 1; // for D/r0 = 1 here because of cor2cov that take into account the r0
	
	nCe		= cor2cov(extra, normvar); // normalized at D/r0 = 1
	
	nCe(:extra.nact, :extra.nact) *= 0.; // mirror space not affected by r0
		
	//dr0		= -(5. / 3) * extra.teldiam^(5. / 3) * nCe * r0^(-8. / 3); 

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
 
 */
{
	niter	= dimsof(eps_ortho)(3);
	pix		= where(ipupil);
	//kl		= kl(*, )(pix, );
	
	dimp	= dimsof(cMat)(2);
	dimo	= dimsof(Dpinf)(3);
	
	//c1		= c2 = c3 = array(float, [2, dimp + dimo, niter]);
	c			= array(float, [2, dimp + dimo, niter]);
	
	write, ""; write, "Computing aliased commands"; write, "";
	com_ortho	= (cMat(, +) * (Dpinf(, +) * eps_ortho(+, ))(+, ));
	
	write, "Computing terms";
	t1	= t2 = t3 = array(float, [2, dimp, niter]);
	
	for (i=n+1 ; i<=niter ; ++i) {
		//t1(, i)		= (1 - g)^n * 	eps_para(, i-n)
		
		for (m=1 ; m<=n ; ++m) {
			t2(, i)	+= (1 - g)^(m - 1) * com_ortho(, i - m)
			t3(, i)	+= (1 - g)^(m - 1) * (phi_para(, i - m + 1) - phi_para(, i - m));
		}
	}
	t2	*= -g;
	
	
	write, "Computing cross terms";
	//c1(:dimp, )		= t1; 
	//c2(:dimp, )		= t2;
	//c3(:dimp, )		= t3;
	//c1(dimp+1:, )	= c2(dimp+1:, ) = c3(dimp+1:, ) = eps_ortho;
	//C_cross1		= c1(, +) * c1(, +) / niter;
	//C_cross2		= c2(, +) * c2(, +) / niter;
	//C_cross3		= c3(, +) * c3(, +) / niter;
	//error, "in calc_cor_term";
	//write, "Computing correlations";
	//C_cross_cor1		= cov2cor(C_cross1);
	//C_cross_cor2		= -cov2cor(C_cross2);
	//C_cross_cor3		= cov2cor(C_cross3);
	
	c(:dimp, )		= t2 + t3;
	c(dimp+1:, )	= eps_ortho;
	
	C_cross			= c(, n+1:)(, +) * c(, n+1:)(, +) / niter;
	//C_cross	= [C_cross1, C_cross2, C_cross3];
	//error;
	
	return C_cross;
}


func hist(t,over=,color=, pasZero=, marks=)
/* DOCUMENT hist(tab,over=,color=, marks=, pasZero=)
 displays the histogram of array 'tab'.
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

func minim_manu(cov, extra, istart, istop, ntot, nmodes, cross=, n=, g=) 
{
	
	//extra	= extra_struct();
	
	if (cross && n && g) {
		cor		= cov2cor(cov, var);
		cross	= calc_cor_term(n, g, ipupil, act_all(:nmodes,), cMat, Dtot(,nmodes+1:), act_all(nmodes+1:,), cbphase(:nmodes,));
		
		cross	= cov2cor(cross, var2);
		//cross	= cor2cov(cross, var);
		//cross	= cross(nmodes + 1:, :nmodes);

		//extra(nmodes + 1:, :nmodes)	= cross;
		//extra(:nmodes, nmodes + 1:)	= transpose(cross);
		extra.cor	= &(cross);
		
		//error;
	} else {
		extra.cor	= &(cov2cor(cov, var));
	}
	
	//extra.var_ortho	= &(diag(cov)(nmodes+1:));
	//extra.nact		= nact;
	//extra.normJ		= 1.e7;
	//extra.teldiam	= 4.2;
	
	tmp	= [];
	window, 5;
	fma;
	x	= span(istart, istop, ntot);
	c	= cov;
	
	vextra	= c;
	vextra(:30,:30) -= diag(diag(c(:30,:30)));
	
	
	for(s=1 ; s<=nmodes ; s+=1){
		write, format="\rMode #%d", s;
		limits;
		t	= x * 0;
		v	= get_vec_diag(c);
		
		for (j=1 ; j<=numberof(x) ; ++j){
			v(s)		*= x(j);
			vextra(,s)	*= sqrt(x(j));
			vextra(s,)	*= sqrt(x(j));
			t(j)		= J1J2_diag(v,,extra);
			v(s)		/= x(j);
			vextra(s,)	/= sqrt(x(j));
			vextra(,s)	/= sqrt(x(j));
		} 
		
		plg, t, x, color="green", marks=0; 
		m=where(t==min(t));
		grow, tmp, x(m);
	}
	
	if (!window_exists(1)) {
		window, 6, dpi=120;
	} else window, 6;
	fma;
	plg, v, marks=0, width=3;
	plg, v * tmp, marks=0, color="red";
	xytitles, "modes #", "Variance [V^2^]";
	pltitle, "Estimation de la variance";
	write, "";
	
	return tmp;
}



func minim_auto(cee_vraie, extra, var_guess, nmodes, cross=, n=, g=, mc_cov=) 
{
	//extra	= extra_struct();
	
	if (cross && n && g) {
		cor		= cov2cor(cee_vraie, var);
		var		= var(:nmodes);
		cross	= calc_cor_term(n, g, ipupil, act_all(:nmodes,), cMat, Dtot(,nmodes+1:), act_all(nmodes+1:,), cbphase(:nmodes,));
		
		cross		= cov2cor(cross, var2);
		//cross	= cor2cov(cross, var);
		extra.cor	= &(cross);
	
		//error;
	} else {
		extra.cor	= &(cov2cor(cee_vraie, var));
	}
	
	v	= var(:nmodes);
	
	//extra.cor		= &(cov2cor(cee_vraie, var));
	//extra.var_ortho	= &(diag(cee_vraie)(nmodes+1:));
	//extra.nact		= nmodes;
	//extra.normJ		= 1.e7;
	//extra.teldiam	= 4.2;
	
	//op_check_gradient(J1J2_diag, grad_J1J2_diag, var_guess, extra=extra, tiny=0.1);
	
	var_guess	= v*1.1;
	if (is_void(var_guess)) var_guess	= get_param(mc_cov, 32., extra);

	//(*extra.cor)(:30,:31) *=0;
	//(*extra.cor)(31:,:30) *=0;
	
	
	tmp	= op_mnb(J1J2_diag_nor0, var_guess, fout, gout, extra=extra, xmin=1.e-3, verb=5, ftol=1.e-15, gtol=1.e-15);
	
	
	if (!window_exists(6)) {
		window, 6, dpi=120;
	} else window, 6;
	//fma;
	//plg, v, marks=0, width=3;
	plg, abs(v-tmp(:nmodes)) * (dm(1).unitpervolt)^2 * 1.e6, marks=0, color="red";			// multiplied to convert V^2 into microns^2 into nm^2
	if (!is_void(mc_cov)) plg, abs(v-diag(mc_cov)) * (dm(1).unitpervolt)^2 * 1.e6, marks=0, color="blue";
	limits,,,0,;
	xytitles, "modes #", "\|!s^2^_simu_ - !s^2^_esti_ \| [nm^2^]",[-0.01,0];
	//pltitle, "Estimation de la variance";
	//write, "r0 = ", extra.teldiam / tmp(0);
	write, "r0 not fitted";
	
	return tmp;
}






func minim_auto_all(cee_vraie, extra, var_guess, nmodes, cross=, n=, g=, mc_cov=) 
{
	//extra	= extra_struct();
	nact	= extra.nact

	if (cross && n && g) {
		cor		= cov2cor(cee_vraie, var);
		var		= var(:nmodes);
		cross	= calc_cor_term(n, g, ipupil, act_all(:nmodes,), cMat, Dtot(,nmodes+1:), act_all(nmodes+1:,), cbphase(:nmodes,));
		
		cross		= cov2cor(cross, var2);
		//cross	= cor2cov(cross, var);
		extra.cor	= &(cross);
		
		(*extra.cor)(:nmodes, :nmodes) *= 0.;
		//error;
	} else {
		extra.cor	= &(cov2cor(cee_vraie, var));
		(*extra.cor)(:nmodes, :nmodes) *= 0.;
	}
	
	v	= cov2param(cee, 35., extra);
	
	//extra.cor		= &(cov2cor(cee_vraie, var));
	//extra.var_ortho	= &(diag(cee_vraie)(nmodes+1:));
	//extra.nact		= nmodes;
	//extra.normJ		= 1.e7;
	//extra.teldiam	= 4.2;
	
	//op_check_gradient(J1J2_diag, grad_J1J2_diag, var_guess, extra=extra, tiny=0.1);
	
	//var_guess	= v;
	
	if (is_void(var_guess)) {
		var_guess	= cov2param(cee, 35., extra);
	}
		
	xmin	= xmax = var_guess * 0.;
	
	xmin(:nact)		= 1.e-3;				// positive variances
	xmin(nact+1:)	= -1.;					// min correlation coefficient 
	xmin(0)			= 0.9 * var_guess(0);	// min D/r0
	xmax(:nact)		= 1000.;				// max variances
	xmax(nact+1:)	= 1.;					// max correlation coefficient
	xmax(0)			= 1.1 * var_guess(0);	// max D/d0
	
	tmp	= op_mnb(J1J2, var_guess, extra=extra, verb=1, /*ftol=1.e-14, gtol=1.e-15,*/ xmin=xmin, xmax=xmax);
	
	
	if (!window_exists(6)) {
		window, 6, dpi=120;
	} else window, 6;
	//fma;
	plg, abs(v - tmp)(:numberof(v)-1), marks=0, color="red";
	//plg, tmp, marks=0, color="red";
	if (!is_void(mc_cov)) plg, abs(v - cov2param(mc_cov, 35., extra))(:numberof(v)-1), marks=0, color="blue";
	xytitles, "Parameters", "";
	pltitle, "Performance of the estimation";
	write, "r0 = ", extra.teldiam / tmp(0);
	write, "";
	
	return tmp;
}






func test_n(cov, nmodes, g, act_all, cMat, Dtot, nmin, nmax)
{
	
	for (n=nmin ; n<=nmax ; ++n) {

		winkill, n;
		//cor		= cov2cor(cov, var);
		//cor		= [];
		
		cross	= calc_cor_term(n, g, ipupil, act_all(:nmodes,), cMat, Dtot(,nmodes+1:), act_all(nmodes+1:,), cbphase(:nmodes,));
	
		cross	= cross(, , 1:3)(, , sum);
		//cross	= cov2cor(cross, var2);
		//cross	= cor2cov(cross, var)(nmodes + 1:, :nmodes);
	
		window, n, dpi=110;
		plp, cov(nact+1:,:nact)(*), cross(nact+1:,:nact)(*), size=0.1;
		plg, span(-9,9,2), span(-9,9,2), marks=0, color="red";
		pltitle, "n = " + strtrim(swrite(n));
	}
	
	return;
}



func test_all(n, g, cov, nmodes, istart, istop, ntot)
{
	cor		= cov2cor(cov, var);
	grow, var, 35.
	
	cross	= calc_cor_term(n, g, ipupil, act_all(:nmodes,), cMat, Dtot(,nmodes+1:), act_all(nmodes+1:,), cbphase(:nmodes,));
	
	cross		= cov2cor(cross, var2);
	extra.cor	= &(cross);
	cross	= cor2cov(extra, var);
	
	
	window, 1;
	fma;
	plp, cee(nmodes + 1:, :nmodes)(*), cross(nmodes + 1:, :nmodes)(*), size=0.1;
	plg, span(-9, 9, 2), span(-9, 9, 2), marks=0, color = "red";
	pltitle, "Performance de l\'estimation pour n=" + strtrim(swrite(n));
	xytitles, "C_!e!]!|\|_ (simulation)", "C_!e!]!|\|_ (estimation)";

	window, 2;
	fma; limits;
	hist((cee(nmodes + 1:, :nmodes)(*) - cross(nmodes + 1:,:nmodes)(*)), marks=0, over=1);

	tmp = minim_manu(cov, istart, istop, ntot, nmodes, cross=1, n=n, g=g)

	return;
}

func test_r0(n, g, cov, extra, nmodes, istart, istop, ntot, cross=)
{

	//extra	= extra_struct();
	
	if (cross && n && g) {
		cor		= cov2cor(cov, var);
		cross	= calc_cor_term(n, g, ipupil, act_all(:nmodes,), cMat, Dtot(,nmodes+1:), act_all(nmodes+1:,), cbphase(:nmodes,));
		
		cross	= cov2cor(cross, var2);
		//cross	= cor2cov(cross, var);
		//cross	= cross(nmodes + 1:, :nmodes);
		
		//extra(nmodes + 1:, :nmodes)	= cross;
		//extra(:nmodes, nmodes + 1:)	= transpose(cross);
		extra.cor	= &(cross);
		
		//error;
	} else {
		extra.cor	= &(cov2cor(cov, var));
	}
	
	//extra.var_ortho	= &(diag(cov)(nmodes+1:));
	//extra.nact		= nact;
	//extra.normJ		= 1.e7;
	//extra.teldiam	= 4.2;
	
	tmp	= [];
	window, 5;
	fma;
	x	= span(istart, istop, ntot);
	c	= cov;
	
	limits;
	t	= x * 0;
	
	for (j=1 ; j<=numberof(x) ; ++j) {
		v		= get_param(c, x(j), extra);
		t(j)	= J1J2_diag(v,,extra);
	} 
	
	plg, t, x, color="green", marks=0;
	write, "r0 = ", x(where(t == min(t)));
	return;
}



func test_cij(n, g, cov, extra, nmodes, istart, istop, ntot, coef, cross=, Dr0=)
{
	
	//extra	= extra_struct();
	//winkill, 11;
	
	if (cross && n && g) {
		cor		= cov2cor(cov, var);
		cross	= calc_cor_term(n, g, ipupil, act_all(:nmodes,), cMat, Dtot(,nmodes+1:), act_all(nmodes+1:,), cbphase(:nmodes,));
		
		cross	= cov2cor(cross, var2);
		//cross	= cor2cov(cross, var);
		//cross	= cross(nmodes + 1:, :nmodes);
		
		//extra(nmodes + 1:, :nmodes)	= cross;
		//extra(:nmodes, nmodes + 1:)	= transpose(cross);
		extra.cor	= &(cross);
		
		//error;
	} else {
		extra.cor	= &(cov2cor(cov, var));
	}
	
	//extra.var_ortho	= &(diag(cov)(nmodes+1:));
	//extra.nact		= nact;
	//extra.normJ		= 1.e7;
	//extra.teldiam	= 4.2;
	
	tmp	= [];

	//fma; limits;
	window, 5;
	animate, 1;

	x	= span(istart, istop, ntot);
	c	= cov;
	
	limits;
	t	= x * 0;
	
	if (is_void(Dr0)) Dr0	= 35.;
	
	
	for (k=1 ; k<=numberof(x) ; ++k) {
		//(*extra.cor)(i,j) *= x(k);
		//(*extra.cor)(j,i) *= x(k);
		v		= cov2param(c, Dr0, extra);
		v(coef)	*= x(k);
		fma; pli, param2cov(v,extra)(:30,:30);pltitle,strtrim(swrite(v(coef)));
		t(k)	= J1J2(v,,extra);
		//(*extra.cor)(i,j) /= x(k);
		//(*extra.cor)(j,i) /= x(k);
		v(coef)	/= x(k);
	} 
	
	window, 11; plg, t, x, color="red", marks=0;redraw;
	
	return;
}









func op_check_gradient(f, g, x0, number=, tiny=, dir=, extra=)
/* DOCUMENT op_check_gradient(f, g, x)
 
 This function a part of the optimpack-mod.i package
 
 Compare gradient function with gradient estimated by finite
 differences.  F is the function, G is the gradient and X the
 parameters.
 
 The number of gradient values to check may be specified by keyword
 NUMBER (the subset of parameters is randomly chosen). The default is
 to compute the finite difference gradient for all parameters.
 
 The (relative) finite difference step size cand be specified by
 keyword TINY (default 1e-5).
 
 Keyword DIR can be used to specify which kind of finite differences to
 use to estimate the gradient: forward (DIR>0), backward (DIR<0) or
 centered (DIR=0 or nil).  The default is to use centered finite
 differences which are more precise (of order TINY^3) but twice more
 expensive to compute.
 
 SEE ALSO: op_mnb */
{
	n = numberof(x0);
	if (is_void(number) || number >= n) {
		list = indgen(n);
		number = n;
	} else {
		list = long(1 + n*random(number));
		list = indgen(31:465);
	}
	if (is_void(tiny)) tiny = 1e-5;
	x1 = x0;
	f0 = f(x0, gx, extra);
	
	g0 = gx(list);
	
	g1 = array(double, number);
	if (dir) {
		/* Forward or backward differences. */
		tiny = abs(tiny)*sign(dir);
		for (i=1 ; i<=number ; ++i) {
			l = list(i);
			x = x0(l);
			x1(l) = x + (dx = tiny*max(1.0, abs(x)));
			g1(i) = (f(x1, , extra) - f0)/dx;
			x1(l) = x;
		}
	} else {
		/* Centered differences. */
		for (i=1 ; i<=number ; ++i) {
			l = list(i);
			x = x0(l);
			dx = tiny*max(1.0, abs(x));
			x1(l) = x + dx;
			f1 = f(x1, , extra);
			x1(l) = x - dx;
			f2 = f(x1, , extra);
			x1(l) = x;
			g1(i) = (f1 - f2)/dx;
			write, format="\r %d", i;
		}
		write, "";
		g1 *= 0.5;
	}
	abs_err = abs(g1 - g0);
	max_abs_err = max(abs_err);
	avg_abs_err = avg(abs_err);
	rms_abs_err = sqrt(avg(abs_err*abs_err));
	u = abs(g0) + abs(g1);
	rel_err = 2*(g1 != g0)*abs_err/(u + !u);
	max_rel_err = max(rel_err);
	avg_rel_err = avg(rel_err);
	rms_rel_err = sqrt(avg(rel_err*rel_err));
	my			= abs(g1-g0)/abs(g0 +!g0);
	max_my		= max(my);
	min_my		= min(my)
	avg_my		= avg(my);
	
	write, format="GRADIENT CHECK WITH: tiny=%.1e,  number=%d\n",
    double(tiny), number;
	write, format="ABSOLUTE ERROR: max=%.1e,  avg=%.1e,  rms=%.1e\n",
    max_abs_err, avg_abs_err, rms_abs_err;
	write, format="RELATIVE ERROR: max=%.1e,  avg=%.1e,  rms=%.1e\n",
    max_rel_err, avg_rel_err, rms_rel_err;
	write, format="RELATIVE ERROR: max=%.1e,  min=%.1e,  avg=%.1e\n",
    max_my, min_my, avg_my;
	winkill, 10;
	window,10,dpi=140;
	plg,g0,marks=0;
	plg,g1,marks=0,color="red";
	//logxy, 0, 1;
	write, "black : analitical gradient";
	write, "red   : finite difference";
	
	error, "check";
}


struct extra_struct
{
	pointer	cor;												// matrix of cross terms and orthogonal variance
	pointer var_ortho;
	pointer cww;
	int		nact;
	float	normJ;
	float	teldiam;
}
