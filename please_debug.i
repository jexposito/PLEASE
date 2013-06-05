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