require, "please_utils.i";
require, "please_debug.i";
require, "yao.i";
require, "yao_funcs.i";
require, "optimpack-mod.i";

func psf_please(cee, extra)
/*	DOCUMENT
 
 
 */
{
	/*	Inits	*/
	pupil			= *extra.ipupil;
	Mi				= *extra.modes;
	nmodes			= dimsof(cee)(2);
	lambdaim		= extra.lambdaim;
	lambdawfs		= extra.lambdawfs;
	otf				= correlate(pupil, pupil).re;
	dphi			= pupil * 0.;
	
	cee	*= extra.normModes;
	Mi	/= sqrt(extra.normModes);		// Mi . Mi = extra.normModes
	
	/*	Uij	*/
	/*for (i=1 ; i<=nmodes ; ++i) {
		for (j=i ; j<=nmodes ; ++j) {
			
			write, format="\rComputing U%d%d", i, j;
			
			uij	= calc_Uij(Mi(.., i), Mi(.., j), pupil);
			
			if (i != j) {
				dphi	+= (cee(i, j) + cee(j, i)) * uij;
			} else {
				dphi	+= cee(i, j) * uij;
			}
		}
	}
	/*
	for (i=nm+1 ; i<=nminf ; ++i) {
		write, format="\rComputing U%d%d", i, i;
		
		uij		= calc_Uij(Mi(.., i), Mi(.., i), pupil);
		dphi	+= cee(i, i) * uij;
	}*/
	
	
	/*	Reconstruction of the structure function of the phase using Vij algorithm	*/
	/*	Diagonalisation	*/
	l = SVdec(cee, u, vt, full=1);	// reconstruction: cee = (u(,+) * diag(l)(+,))(,+) * vt(+,)
	
	write, "\nVii computation\n";
	vii		= Mi(, , +) * vt(,+);
	ftmi	= complex(vii * 0.);
	
	/*	Computation of the FFT of the modes	*/
	for (i=1 ; i<=nmodes ; i++) {
		ftmi(, , i)	= fft(vii(, , i), 1);
		write, format=" \rComputing fft of mode %d", i;
	}	
	write, "";
	
	/*	Vii computation and structure function computation	*/
	for (i=1 ; i<=nmodes ; i++) {
		write, format=" \rComputing V%d%d mode", i;
		modei	= vii(..,i);
		
		tmp		= calc_Viif(ftmi(, , i), ftmi(, , i), modei, modei, otf, conj(fft(ipupil)));
		
		dphi	+= tmp * l(i);
	}
	write, "\n";
	
	/*	Reconstruction of the OTF	*/

	fact					= (2 * pi / lambdaim)^2;						// conversion lambda_measurements to lambda_images
	tmp						= exp(-0.5 * dphi * fact);
	
	/*	Support definition	*/
	mask					= otf > max(otf) * 1.e-7;
	mask					= 1 - mask;
	tmp(where(mask))		= 0.;	
	sz						= dimsof(tmp)(2);
	fto_turb				= array(float, [2, 512, 512]);
	fto_turb(1:sz, 1:sz)	= eclat(tmp);
	
	fto_turb				= eclat(roll(fto_turb, [512 / 2 - sz / 2, 512 / 2 - sz / 2]));
	
	/*	Various PSF : Telescope / on sky	*/
	fto_tel				= telfto(0.589, 0.01, 4.2, 0.25, .589/4.2/4.85/(float(512)/140), 512) // just for Canary/WHT
	psftel				= eclat(abs(fft(fto_tel, -1)));
	psf					= eclat(abs(fft(fto_tel * fto_turb, -1)));
	
	/*	PSF from yao	*/
	imav				= *extra.imav;
	psftest				= imav / sum(imav);
	
	/* re-normalization of the reconstructed PSF and the telescope PSF	*/
	psfrec				= psf / sum(psf);
	psftel				= psftel / sum(psftel);
	
	difract				= circavg(psftel, middle=1);
	
	/*	Display	*/
	window, 11, dpi=130; // PB after computation : WARNING Gist GdText plotter failed
	pause, 100;
	plg,[1]; pause, 100; redraw;
	fma; limits, , 10;
	plg, difract / max(difract), marks=0, width=3;
	plg, circavg(eclat(psftest)) / max(difract), marks=0, color="blue";
	plg, circavg(eclat(psfrec)) / max(difract), color="red", marks=0;
	xytitles, "Pixels", "Strehl ratio";
	pltitle, "PSF circ avg";
	
	error, "test dphi";
	
	return psf;
}

func covar_please(cee_guess, D_r0_guess, cnn_guess, extra)
/* DOCUMENT
 
 
 
*/
{
	
	cee	= cee_guess;//develop
	
	return cee;
}

func test_please(void)
/* DOCUMENT
 
 
 
 
 */
{
	extern smes_nonoise, act_all, cbphase;
	
	aoread, "canary-TS-zernike.par";
	aoinit, clean=1, forcemat=1;
	aoloop, disp=0, savecb=1;
	
	// Inits for simulation
	nzer	= 200; // First 200 - piston
	size	= sim._size;
	nmodes	= dm(1)._nact;
	mInf	= array(0.0f, size, size, nzer);
	cpt		= 0;
	
	tmp		= fits_read("new_kl.fits"); // 200 orthogonal KL modes not normalized to 1  i.e. (Mi . Mi) != 1
	
	// computing matrix of modes
	nm		= 1;
	n1		= dm(nm)._n1;
	n2		= dm(nm)._n2;
	for (cc=1 ; cc<=nzer ; cc++) {
		cpt++;
		mInf(n1:n2, n1:n2, cpt)	= tmp(n1:n2, n1:n2, cpt);
	}
	
	// Computing the projection matrix from phase screen to modes i.e. phase2modes
	// and init some variables
	valid_pix	= where(ipupil);
	tmp			= mInf(*, )(valid_pix, );
	matPass		= LUsolve(tmp(+, ) * tmp(+, ))(+, ) * tmp(*, )(, +);
	
	smes_nonoise	= cbmes;
	act_all			= array(0.0f, [2, nzer, loop.niter]);
	cbphase			= array(float, [2, nzer, loop.niter]);

	go, all=1;
	
	// inits for reconstruction
	extra	= extra_struct();
	
	Dtot	= fits_read("iMat_tot_geom.fits"); // implement computation of iMat with 200 modes
	g		= loop.gain;
	n		= int(ceil(-6. / log(1 - g)));
	
	cross	= calc_cor_term(n, g, ipupil, act_all(:nmodes,), cMat, Dtot(,nmodes+1:), act_all(nmodes+1:,), cbphase(:nmodes,));
	
	cee			= act_all(, +) * act_all(, +) / loop.niter;
	var_ortho	= diag(cee)(nmodes+1:) / atm.dr0at05mic^(5. / 3);
	
	/*	normalisation factor so that (Mi . Mi) = 1	*/
	m			= mInf(*, )(valid_pix, );
	normfact	= diag(m(+, ) * m(+, ))(avg);
	
	
	/*		*/
	extra.cor		= &(cross);
	extra.var_ortho	= &(var_ortho);
	extra.cww		= &(cbmes(, +) * cbmes(, +) / loop.niter);
	extra.cMat		= &cMat;
	extra.iMat		= &Dtot;
	extra.modes		= &(mInf);
	extra.cbact		= &act_all;
	extra.ipupil	= &ipupil;
	extra.imav		= &(imav(.., 1, 1));
	extra.nact		= dm(1)._nact;
	extra.teldiam	= tel.diam;
	extra.niter		= loop.niter;
	extra.normJ		= 1.e7;
	extra.normModes	= normfact;
	extra.lambdawfs	= wfs(1).lambda;
	extra.lambdaim	= (*target.lambda)(1);
	
	cnn	= 0.;
	// only fitting the variance and D/r0
	write, "";
	write, "Minimization of the criterion";
	
	tmp	= op_mnb(J1J2_diag, get_param(cee, atm.dr0at05mic, extra), fout, gout, extra=extra, verb=5);
	
	write, "";
	write, "r0 = ", extra.teldiam / tmp(0);
	write, "";
	
	cee	= cor2cov(tmp,extra);
	psf	= psf_please(cee, extra);
	
	return;	
}
