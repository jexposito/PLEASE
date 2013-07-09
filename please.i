require, "please_utils.i";
require, "please_debug.i";
require, "yao.i";
require, "yao_funcs.i";
require, "optimpack-mod.i";

func psf_please(cee, extra, MC=)
/*	DOCUMENT
 
 
 */
{
	/*	Inits	*/
	pupil			= *extra.ipupil;
	Mi				= *extra.modes;
	nmodes			= dimsof(cee)(2);
	nact			= extra.nact;
	lambdaim		= extra.lambdaim;
	lambdawfs		= extra.lambdawfs;
	otf				= correlate(pupil, pupil).re;
	dphi			= pupil * 0.;
	cMat			= *extra.cMat;
	Dplus			= *extra.Dplus
		
	/*	Reconstruction of the structure function of the phase using Vij algorithm	*/
	/*	Diagonalisation	*/
	l		= SVdec(cee, u, vt, full=1);	// reconstruction: cee = (u(,+) * diag(l)(+,))(,+) * vt(+,)
	
	write, "\nVii computation\n";
	Mii		= Mi(, , +) * vt(, +);	// diagonalisation
	ftmi	= complex(Mii * 0.);	//
	
	/*	Computation of the FFT of the modes	for fast computing */
	for (i=1 ; i<=nmodes ; i++) {
		ftmi(, , i)	= fft(Mii(, , i), 1);
		write, format=" \rComputing FFT of mode #%d", i;
	}	
	write, "";
	
	conjftpup	= conj(fft(ipupil));
	/*	Vii computation and structure function computation	*/
	for (i=1 ; i<=nmodes ; i++) {
		write, format=" \rComputing V%d%d mode", i;
		modei	= Mii(.., i);
		
		vii		= calc_Viif(ftmi(, , i), ftmi(, , i), modei, modei, otf, conjftpup);
		
		dphi	+= l(i) * vii;
	}
	write, "\n";
	
	if (MC) {
		//cee_p	= cMat(, +) * ((*extra.cww)(, +) * cMat(, +))(+, );
		cee_p	= Dplus(, +) * ((*extra.cww)(, +) * Dplus(, +))(+, )
		cee_mc	= cee_p + *extra.alias;
		l		= SVdec(cee_mc, u, vt, full=1);			// reconstruction: cee = (u(,+) * diag(l)(+,))(,+) * vt(+,)
		
		write, "\nVii computation for MC method\n";
		Mii		= Mi(, , :extra.nact)(, , +) * vt(,+);	// 
		ftmi	= complex(Mii * 0.);					//
		
		/*	Computation of the FFT of the modes	for fast computing */
		for (i=1 ; i<=nact ; i++) {
			ftmi(, , i)	= fft(Mii(, , i), 1);
			write, format=" \rComputing FFT of mode #%d", i;
		}	
		write, "";
		
		/*	Vii computation and structure function computation	*/
		dphi_mc	= dphi * 0.;
		for (i=1 ; i<=nact ; i++) {
			write, format=" \rComputing V%d%d mode for MC method", i;
			
			modei	= Mii(.., i);
			vii		= calc_Viif(ftmi(, , i), ftmi(, , i), modei, modei, otf, conjftpup);

			dphi_mc	+= l(i) * vii;
		}
		write, "\n";
	}
	
	/*	Reconstruction of the OTF	*/
	fact					= (2 * pi / lambdaim)^2;		// conversion µm^2 to rad^2
	tmp						= exp(-0.5 * dphi * fact);
	
	/*	Support definition	*/
	mask					= otf > max(otf) * 1.e-7;
	mask					= 1 - mask;
	tmp(where(mask))		= 0.;	
	sz						= dimsof(tmp)(2);
	otf_turb				= array(float, [2, 512, 512]);
	otf_turb(1:sz, 1:sz)	= eclat(tmp);
	
	otf_turb				= eclat(roll(otf_turb, [512 / 2 - sz / 2, 512 / 2 - sz / 2]));
		
	/*	Various PSF : Telescope / on sky	*/
	otf_tel				= telfto(lambdaim, 0.01, extra.teldiam, extra.cobs, lambdaim / extra.teldiam / 4.848 / (float(sim._size) / sim.pupildiam), sim._size) // just for Canary/WHT
	psftel				= eclat(abs(fft(otf_tel, -1)));
	psf					= eclat(abs(fft(otf_tel * otf_turb, -1)));
	
	
	if (MC)	{
		tmp					= exp(-0.5 * dphi_mc * fact) * exp(-0.5 * eclat((*extra.Dphi_ortho)) * fact);
		tmp(where(mask))	= 0.;
		otf_mc				= otf_turb * 0.;
		otf_mc(:sz, :sz)	= eclat(tmp);
		otf_mc				= eclat(roll(otf_mc, [512 / 2 - sz / 2, 512 / 2 - sz / 2]));
		psf_mc				= eclat(abs(fft(otf_tel * otf_mc, -1)));
	}
		
	/*	PSF from Yao	*/
	imav				= *extra.imav;
	psftest				= imav / sum(imav);
	
	/* re-normalization of the reconstructed PSF and the telescope PSF	*/
	psfrec				= psf / sum(psf);
	psftel				= psftel / sum(psftel);
	if (MC)	psf_mc		= psf_mc / sum(psf_mc)
	
	difract				= circavg(psftel, middle=1);
	
	/*	Display	*/
	if (!window_exists(11)) window, 11, dpi=130, style="boxed.gs"; // PB after computation : WARNING Gist GdText plotter failed
	pause, 100;
	plg,[1]; pause, 100; redraw;
	fma; limits, , 10;
	plg, difract / max(difract), marks=0, width=3, type="dash";
	plg, circavg(eclat(psftest)) / max(difract), marks=0, color="red", width=3;
	plg, circavg(eclat(psfrec)) / max(difract), color="blue", marks=0;
	if (MC) plg, circavg(eclat(psf_mc)) / max(difract), color="green", marks=0;
	xytitles, "Pixels", "Strehl ratio";
	pltitle, "PSF circ avg";
	
	write, "";
	write, "Dashed: perfect telescope\n"
	write, "Red   : Observation";
	write, "Blue  : MV reconstruction";
	if (MC) write, "Green : MC reconstruction";
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

func test_please(MC=)
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
	
	cross	= calc_cor_term(n, g, ipupil, act_all(:nmodes ,), cMat, Dtot(, nmodes+1:), act_all(nmodes+1:, ), cbphase(:nmodes, ));
	
	cee			= act_all(, +) * act_all(, +) / loop.niter;
	var_ortho	= diag(cee)(nmodes+1:) / atm.dr0at05mic^(5. / 3);
	
	/*	normalisation factor so that (Mi . Mi) = 1	*/
	m			= mInf(*, )(valid_pix, );
	normfact	= diag(m(+, ) * m(+, ))(avg);
	normfact	= sqrt(normfact);
	
	if (MC) {
		dphi_ortho	= fits_read("dphi_ortho_canary.fits") / 2000.; // not normalized by the number of iteration when recorded niter was 2000
		
		/*	Aliasing on the mirror for MC computation	*/
		alias		= cMat(, +) * (Dtot(, dm(1)._nact+1:)(, +) * act_all(dm(1)._nact+1:, )(+, ))(+, );
	}

	/*	Filling extra	*/
	extra.cor			= &(cross);
	extra.var_ortho		= &(var_ortho);
	extra.cww			= &(cbmes(, +) * cbmes(, +) / loop.niter);
	extra.cMat			= &cMat;
	extra.iMat			= &Dtot;
	extra.Dplus			= &(LUsolve(iMat(+, ) * iMat(+, ))(, +) * iMat(, +));
	extra.modes			= &(mInf);
	extra.cbact			= &act_all;
	extra.ipupil		= &ipupil;
	extra.imav			= &(imav(.., 1, 1));
	extra.nact			= dm(1)._nact;
	extra.teldiam		= tel.diam;
	extra.cobs			= tel.cobs;
	extra.niter			= loop.niter;
	extra.normJ			= 1.e7;
	extra.normModes		= normfact;
	extra.lambdawfs		= wfs(1).lambda;
	extra.lambdaim		= (*target.lambda)(1);
	
	if (MC) {
		extra.Dphi_ortho	= &(dphi_ortho * atm.dr0at05mic^(5./3));					// only for MC method
		extra.alias			= &(alias(, +) * alias(, +) / extra.niter);
	}
	
	cnn	= 0.;
	// only fitting the variance and D/r0
	write, "";
	write, "Minimization of the criterion";
	write, "";
	
	tmp	= op_mnb(J1J2_diag, get_param(cee, atm.dr0at05mic+0.1, extra), fout, gout, extra=extra, verb=5, xmin=1.e-3);
	
	write, "";
	write, "r0 = ", extra.teldiam / tmp(0);
	write, "";
	
	cee	= cor2cov(extra, tmp);			// Cee reconstruction with estimated parameters
	
	psf	= psf_please(cee, extra, MC=MC);
	
	return;	
}