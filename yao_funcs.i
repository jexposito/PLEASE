/*
 * yao_funcs.i
 *
 * modified yao functions for please.i
 *
*/

require, "iterkolmo.i";

func mult_wfs(iter, disp=)
/* DOCUMENT func mult_wfs(iter, disp=)
   Goes through all WFS and concatenate the resulting measurement vectors.
   SEE ALSO:
 */
{
	extern wfs;

	if (iter == 1) write, "Using the modified mult_wfs version";
	mes			= [];

	mircube2		= mircube;
	
	for (ns=1 ; ns <= nwfs ; ns++) {
		
		offsets		= wfs(ns).gspos;
		phase		= get_phase2d_from_optics(ns, "wfs");
		phase		+= get_turb_phase(iter, ns, "wfs");

		if (wfs(ns).correctUpTT) {
			phase =	correct_uplink_tt(phase,ns);
		}
		
		//////////////////////////////////////////////////
		//	Start of HACK : filter noise and aliasing	//
		//////////////////////////////////////////////////
		
		pix				= where(ipupil);

		phase			-= phase(*)(pix)(avg);
		cbphase(, iter)	= matPass(, +) * (phase * ipupil)(pix)(+);		// saving the phase at iteration
		phase			= cbphase(, iter)(+) * mInf(,,+);				// High-frequency-filtered turbulent phase only 200 modes

		phase_3			= phase;
		phase_3			-= phase_3(pix)(avg);
				
		mircube			*= 0.0f;
		mircube(, , ns)	= mircube2(, , ns);
		phase_3			+= get_phase2d_from_dms(ns, "wfs");				// remove DM phase
		phase_3			= float(phase_3 - phase_3(pix)(avg));			// remove piston mode

		vect			= matPass(, +) * (phase_3 * ipupil)(pix)(+);	// Projection on KL modes of the residual phase
		act_all(, iter)	= vect;											// Saving 200 KL modes
		

		// Noise free measurements
		nsave			= wfs(ns).noise;
		wfs(ns).noise	= 0;
		
		if (wfs(ns).type == "hartmann" ) {
			tmpmes_all	= sh_wfs(ipupil, phase_3 * ipupil, ns);
		}

		// Subtract the reference vector for this sensor:
		if (wfs(ns)._cyclecounter == 1) {
			tmpmes_all	= tmpmes_all - *wfs(ns)._refmes;
		}

		grow, smes_nonoise, tmpmes_all;			

		wfs(ns).noise	= nsave;
		
		//////////////////
		//	END OF HACK	//
		//////////////////
		
		smes	= sh_wfs(ipupil, phase_3 * ipupil, ns);

		// Subtract the reference vector for this sensor:
		if (wfs(ns)._cyclecounter == 1) {
			smes	= smes - *wfs(ns)._refmes;
		}
		
		wfs(ns)._tt(1)	= sum( smes * (*wfs(ns)._tiprefvn) );
		wfs(ns)._tt(2)	= sum( smes * (*wfs(ns)._tiltrefvn) );
		if (wfs(ns).filtertilt) {
			smes	= smes - wfs(ns)._tt(1) * (*wfs(ns)._tiprefv) \
			- wfs(ns)._tt(2) * (*wfs(ns)._tiltrefv);
		}
		if (wfs(ns)._cyclecounter == 1) {
			wfs(ns)._lastvalidtt	= wfs(ns)._tt;
		}
		
		grow, mes, smes;
	} // end of the loop on WFSs
	
	mircube = mircube2;

	return mes;
}

//----------------------------------------------------


struct screen_struct
{
	float   r0;    
	pointer A;    
	pointer B;   
	pointer istencil;
	float   deltax;    
	float   deltay;    
};

//----------------------------------------------------

func get_turb_phase_init(skipReadPhaseScreens=)
/* DOCUMENT get_turb_phase_init(skipReadPhaseScreens=)
 Initializes everything for get_turb_phase (see below), which
 returns the interpolated, integrated phase to the loop function
 for iteration number "iter". Returns the "ok" parameter which is
 set to 0 if this step is not to be used to compute statistics.
 AUTHOR: F.Rigaut, June 11, 2002.
 SEE ALSO: aoinit, aoloop, get_turb_phase.
 */
{
	extern pscreens,// contains screens, dim [dim Screen X,dim Screen Y, nscreens]
    nscreens,     // number of screens
    noptics,      // number of optics
    optphasemaps, // phase maps for each optics
    xposvec,      // contains X position vs iteration# of center of beam in
	// phase screen N, in pixel coords. dim [N iteracurwtions,nscreens]
    yposvec,      // same for Y. This should be cst in the current implementation
	// as the screens have been cut in Y, therefore are not periodic
	// and therefore we can not wrap.
    wfsxposcub,   // contains the position of the ray at which a given X pixel
	// intersects the Nth phase screen, for each WFS GS.
	// dimension [dimX outphase, nscreens, #GS]
    wfsyposcub,   // Same for Y
    gsxposcub,    // contains the position of the ray at which a given X pixel
	// intersects the Nth phase screen, for each star at which the perf
	// is evaluated dimension [dimX outphase, nscreens, #Star]
    gsyposcub,    // Same for Y
    dmwfsxposcub, // contains the position of the ray at which a given X pixel
	// intersects the Nth mirror, for each WFS GS.
	// dimension [dimX outphase, ndm, #GS]
    dmwfsyposcub, // Same for Y
    dmgsxposcub,  // contains the position of the ray at which a given X pixel
	// intersects the Nth mirror, for each star at which the perf
	// is evaluated dimension [dimX outphase, ndm, #Star]
    dmgsyposcub,  // Same for Y
	
    optwfsxposcub,// contains the position of the ray at which a given X pixel
	// intersects the Nth optics, for each WFS GS.
	// dimension [dimX outphase, noptics, #GS]
    optwfsyposcub,// Same for Y
    optgsxposcub, // contains the position of the ray at which a given X pixel
	// intersects the Nth optics, for each star at which the perf
	// is evaluated dimension [dimX outphase, noptics, #Star]
    optgsyposcub, // Same for Y
	
    xmargins,     // value of xposvec so that no pixel indice < 1 or > dimx(screen)
	// in any off-axis beam and altitude (low margin, up margin)
    ymargins,     // same for Y. We use that in get_turb_phase to determine when to
	// wrap.
    statsokvec,   // 1 if it is ok to collect stats at this iteration.
	// 0 if not, e.g. we just did a jump. dim [iteration]
    inithistory,  // 1 if init has been done.
    screendim;    // [phase screen X dim, Y dim] before they are extended for safe wrapping
	
	extern currentScreenNorm; // current screen normalization. Used when swaping screen.
	extern screens;
	
	// Define a few variables:
	
	nscreens = numberof(*atm.screen);
	if (opt!=[]) noptics = numberof(opt.phasemaps); else noptics=0;
	
	wfsxposcub = wfsyposcub = array(float,[3,_n,nscreens,nwfs]);
	gsxposcub = gsyposcub = array(float,[3,_n,nscreens,target._ntarget]);
	
	dmwfsxposcub = dmwfsyposcub = array(float,[3,_n,ndm,nwfs]);
	dmgsxposcub = dmgsyposcub = array(float,[3,_n,ndm,target._ntarget]);
	
	if (noptics) {
		optwfsxposcub = optwfsyposcub = array(float,[3,_n,noptics,nwfs]);
		optgsxposcub = optgsyposcub = array(float,[3,_n,noptics,target._ntarget]);
	}
	
	//=======================================================
	// READS AND NORMALIZE THE PHASE SCREENS
	// the phase screens now (v2.4) are normalized in microns
	//=======================================================
	
	if (!is_set(skipReadPhaseScreens)) {
		// Compute normalization factor:
		(*atm.layerfrac) = (*atm.layerfrac)/sum(*atm.layerfrac);
		
		weight = float(sqrt(*atm.layerfrac)*(atm.dr0at05mic/
											 cos(gs.zenithangle*dtor)^0.6/sim.pupildiam)^(5./6.));
		weight = weight *0. + 1.;
		
		// above: in radian at 0.5 microns
		weight = weight * float(0.5/(2*pi));
		// ... and now in microns.
		
		
		psize  = tel.diam/sim.pupildiam;  // pixel in meter
		r0tot = (tel.diam/atm.dr0at05mic)/psize;
		
		/*=====================================================
		 How to relate r0(layer) and atm.layerfrac ?
		 r0(i)   = r0 of layer i
		 r0tot   = total r0
		 f(i)    = "fraction" in layer i ( = (*atm.layerfrac)(i) )
		 we have:
		 weight(i) = sqrt(f(i)) * (D/r0tot)^(5/6.) = (D/r0(i))^(5/6.)
		 thus
		 f(i) = (r0tot/r0(i))^(5./3)
		 
		 inversely, we have:
		 r0(i) = r0tot / f(i)^(3./5)
		 =====================================================*/
		
		pscreens = array(float,[3,2*_n,2*_n,nscreens]);
		screen   = screen_struct();
		screens  = array(screen,nscreens);
		psize    = tel.diam/sim.pupildiam;
		// Build the position vector vs iteration by phase screen
		deltax   = sim.pupildiam/tel.diam*(*atm.layerspeed)*cos(dtor*gs.zenithangle)*loop.ittime;
		// Stuff it
		for (i=1;i<=nscreens;i++) {
			if (sim.verbose>=1) {
				write,format="Creating phase screen %d\n",i;
			}
			AB, 2*_n, A, B, istencil;
			screens(i).A        = &A;
			screens(i).B        = &B;
			screens(i).istencil = &istencil;
			// here we need to express r0 in pixels
			screens(i).r0 = r0tot / ((*atm.layerfrac)(i))^(3./5);
			for (j=1;j<=2*_n+1;j++) pscreens(,,i) = extrude(pscreens(,,i),screens(i).r0,A,B,istencil);
			screens(i).deltax   = deltax(i)*cos(dtor*(*atm.winddir)(i));
			screens(i).deltay   = deltax(i)*sin(dtor*(*atm.winddir)(i));
		}
		
		// apply weights to each phase screens (normalize):
		// the screens are expressed in microns
		//pscreens = pscreens*weight(-,-,);
		currentScreenNorm = float(weight);
		
		//=============================================
		// READ THE OPTICS PHASE MAPS
		// the dimension of all the optical phase maps
		// should be the same. These phase maps should
		// be big enough so that the interpolation will
		// not go outside of the provided map. That is,
		// the indices in opt**xposcub (and y), with **
		// being wfs and gs, should not define points
		// outside of the provided phase maps.
		// THIS IS LEFT TO THE RESPONSIBILITY OF THE
		// USER. I do however a simple check below.
		// Normaly, this should not be a problem:
		// if some indices are outside the maps, it means
		// either that the maps are not big enough, in
		// which case we should not use it, or that
		// there was an error in the definition of the
		// altitude.
		//=============================================
		if (noptics) {
			if (sim.verbose) {
				psize  = tel.diam/sim.pupildiam;
				write,"Reading optics. I am expecting phase maps with a pixel";
				write,format=" size of %fm/pixel, as projected in the entrance\n",psize;
				write,"pupil plane. It is also your responsibility to provide phase";
				write,"maps of adequate dimension, i.e. large enough";
			}
			if (sim.verbose) {
				write,format="Reading phase map for optics \"%s\"\n",opt(1).phasemaps;
			}
			tmp = yao_fitsread(YAO_SAVEPATH+opt(1).phasemaps);
			optdims      = dimsof(tmp);
			optdimx      = optdims(2);
			optdimy      = optdims(3);
			
			if (optdimx != optdimy)
				error,"Optics phase maps should be square arrays";
			
			optphasemaps = array(float,[3,optdimx,optdimy,noptics]);
			// Stuff it
			optphasemaps(,,1) = float(tmp);
			
			for (i=2;i<=noptics;i++) {
				if (sim.verbose>=1) {
					write,format="Reading phase map for optics \"%s\"\n",opt(i).phasemaps;
				}
				tmp = yao_fitsread(YAO_SAVEPATH+opt(i).phasemaps);
				if (anyof(dimsof(tmp) != optdims))
					error,"All optics phase maps should have the same dimensions";
				optphasemaps(,,i) = float(tmp);
			}
			opt._cent = optdimx/2.+(sim._cent-sim._size/2.);
			// (sim._cent-sim._size/2.) =  0. or 0.5
		}
	}
	
	//====================================
	// PREPARE GEOMETRY FOR INTERPOLATION
	//====================================
	
	// Build a vector of position (integer position in equivalent
	// iterations:
	iposvec = indgen(loop.niter);
	
	// build a vector of the iterations at which statistics
	// should be accumulated. 1 if ok to accumulate statistics, 0 if not.
	statsokvec   = (indgen(loop.niter)-1) % loop.skipevery;
	statsokvec   = (statsokvec >= loop.startskip);
	
	// has to start away from edge as this is the coordinate of the beam center
	/// will modified later on when we know the beams geometry. see few lines below.
	//xposvec = (1.+iposvec*deltax(-,)*cos(dtor*(*atm.winddir)(-,)));
	//yposvec = (1.+iposvec*deltax(-,)*sin(dtor*(*atm.winddir)(-,)));
	
	xposvec = (1.+deltax*cos(dtor*(*atm.winddir)));
	yposvec = (1.+deltax*sin(dtor*(*atm.winddir)));
	//===========================================================
	// PRE-COMPUTATION OF THE INTERSECT POSITIONS FOR EACH WFS GS
	//===========================================================
	
	// zero-centered vector of position (our reference)
	xref = indgen(_n)-(_n+1)/2.;
	yref = indgen(_n)-(_n+1)/2.;
	
	// loop on WFS
	for (n=1;n<=nwfs;n++) {
		
		// loop on screen
		for (ns=1;ns<=nscreens;ns++) {
			
			// offsets of the center of beam on screen NS
			xoff = (wfs(n).gspos)(1)*4.848e-6*(*atm._layeralt)(ns)/psize;
			yoff = (wfs(n).gspos)(2)*4.848e-6*(*atm._layeralt)(ns)/psize;
			
			// if we are dealing with LGS, there is a geometric
			// factor on the beam size
			if (wfs(n)._gsalt != 0) {
				fact = (wfs(n)._gsalt-(*atm._layeralt)(ns))/wfs(n)._gsalt;
			} else {
				fact = 1.;
			}
			
			// Compute and stuff the array
			wfsxposcub(,ns,n) = xref*fact + xoff;
			wfsyposcub(,ns,n) = yref*fact + yoff;
		}
		
		// loop on DMs
		for (ns=1;ns<=ndm;ns++) {
			
			// offsets of the center of beam on DM NS
			xoff = (wfs(n).gspos)(1)*4.848e-6*(dm.alt)(ns)/psize;
			yoff = (wfs(n).gspos)(2)*4.848e-6*(dm.alt)(ns)/psize;
			
			// if we are dealing with LGS, there is a geometric
			// factor on the beam size
			if (wfs(n)._gsalt != 0) {
				fact = (wfs(n)._gsalt-(dm.alt)(ns))/wfs(n)._gsalt;
			} else {
				fact = 1.;
			}
			
			// Compute and stuff the array
			dmwfsxposcub(,ns,n) = xref*fact + xoff;
			dmwfsyposcub(,ns,n) = yref*fact + yoff;
		}
		
		// loop on optics
		for (ns=1;ns<=noptics;ns++) {
			
			// offsets of the center of beam on optics #NS
			xoff = (wfs(n).gspos)(1)*4.848e-6*opt(ns).alt/psize;
			yoff = (wfs(n).gspos)(2)*4.848e-6*opt(ns).alt/psize;
			
			// if we are dealing with LGS, there is a geometric
			// factor on the beam size
			if (wfs(n)._gsalt != 0) {
				fact = (wfs(n)._gsalt-opt(ns).alt)/wfs(n)._gsalt;
			} else {
				fact = 1.;
			}
			
			// Compute and stuff the array
			optwfsxposcub(,ns,n) = xref*fact + xoff;
			optwfsyposcub(,ns,n) = yref*fact + yoff;
		}
		
	}
	
	// type conversion:
	wfsxposcub = float(wfsxposcub);
	wfsyposcub = float(wfsyposcub);
	dmwfsxposcub = float(dmwfsxposcub);
	dmwfsyposcub = float(dmwfsyposcub);
	optwfsxposcub = float(optwfsxposcub);
	optwfsyposcub = float(optwfsyposcub);
	
	
	//=============================================================
	// PRE-COMPUTATION OF THE INTERSECT POSITION FOR EACH PERF STAR
	//=============================================================
	
	// loop on target
	for (n=1;n<=target._ntarget;n++) {
		
		// loop on screen
		for (ns=1;ns<=nscreens;ns++) {
			
			// offsets of the center of beam on screen NS
			xoff = (*target.xposition)(n)*4.848e-6*(*atm._layeralt)(ns)/psize;
			yoff = (*target.yposition)(n)*4.848e-6*(*atm._layeralt)(ns)/psize;
			
			// note: that's target. we can't be dealing with LGS here (no fact)
			// Compute and stuff the array
			gsxposcub(,ns,n) = xref + xoff;
			gsyposcub(,ns,n) = yref + yoff;
		}
		
		// loop on DMs
		for (ns=1;ns<=ndm;ns++) {
			
			// offsets of the center of beam on DM NS
			xoff = (*target.xposition)(n)*4.848e-6*(dm.alt)(ns)/psize;
			yoff = (*target.yposition)(n)*4.848e-6*(dm.alt)(ns)/psize;
			
			// note: that's target. we can't be dealing with LGS here (no fact)
			// Compute and stuff the array
			dmgsxposcub(,ns,n) = xref + xoff;
			dmgsyposcub(,ns,n) = yref + yoff;
		}
		
		// loop on optics
		for (ns=1;ns<=noptics;ns++) {
			
			// offsets of the center of beam on DM NS
			xoff = (*target.xposition)(n)*4.848e-6*opt(ns).alt/psize;
			yoff = (*target.yposition)(n)*4.848e-6*opt(ns).alt/psize;
			
			// note: that's target. we can't be dealing with LGS here (no fact)
			// Compute and stuff the array
			optgsxposcub(,ns,n) = xref + xoff;
			optgsyposcub(,ns,n) = yref + yoff;
		}
		
	}
	// type conversion:
	gsxposcub = float(gsxposcub);
	gsyposcub = float(gsyposcub);
	dmgsxposcub = float(dmgsxposcub);
	dmgsyposcub = float(dmgsyposcub);
	optgsxposcub = float(optgsxposcub);
	optgsyposcub = float(optgsyposcub);
	
	
	//======================================
	// SOME CHECKS TO AVOID INDICES OVERFLOW
	//======================================
	
	// Now we can modify xposvec and yposvec to make sure we are not going
	// out of the phase screen arrays
	xmargins = abs([min(_(wfsxposcub(*),gsxposcub(*))),max(_(wfsxposcub(*),gsxposcub(*)))]);
	ymargins = abs([min(_(wfsyposcub(*),gsyposcub(*))),max(_(wfsyposcub(*),gsyposcub(*)))]);
	//xposvec = xposvec - min(xposvec) + xmargins(1) +1;
	//yposvec = yposvec - min(yposvec) + ymargins(1) +1;
	xposvec = xposvec * 0. + xmargins(1) + 1;
	yposvec = yposvec * 0. + ymargins(1) + 1;
    
	// type conversion:
	xposvec = float(xposvec);
	yposvec = float(yposvec);
	
	// so now we have everything initiliazed, and we will just have to
	// interpolate the phase screen at points
	// xposvec(iteration,screen#) + wfsxposcub(,screen#,wfs#)
	// and corresponding for Y, and integrate on screen#
	// to get the phase for a given WFS
	// this integration is done by the C routine _get2dPhase
	// which take, in addition to the screens and output phase parameters,
	// only a set of X and Y positions as input [the one we just talked
	// about, xposvec(iteration) + wfsxposcub(,,wfs#) ].
	
	inithistory = 1;
	return 1;
}

//----------------------------------------------------

func get_turb_phase(iter,nn,type)
/* DOCUMENT get_turb_phase(iter,n,type)
 Returns the interpolated, integrated phase along the turbulent phase
 screen data cube to the loop function for iteration number "iter".
 You have to call get_turb_phase_init to initialize prior using.
 SEE ALSO: aoinit, aoloop, get_turb_phase_init.
 */
{
	extern xpos,ypos;
	
	if (!inithistory) {error,"get_turb_phase has not been initialized !";}
	
	sphase = array(float,_n,_n);
	bphase = array(float,sim._size,sim._size);
	
	// Now we can call the C interpolation routine and get the integrated
	// phase for this star
	// there are a few things to do to get ready
	psnx = dimsof(pscreens)(2);
	psny = dimsof(pscreens)(3);
	nscreens = dimsof(pscreens)(4);
	if ((nn == 1) && (type == "wfs")) {
		if (iter > 1) {
			for (cc = 1;cc <= nscreens;cc++) {
				deltax = screens(cc).deltax;
				tmpx = floor(xpos(cc)+deltax) - floor(xpos(cc));
				xpos(cc) += deltax;
				if (abs(tmpx) > 0) { // we jumped of at least of 1 pixel so
					// we need to extrude the phase screen and remove the number
					// of pixels we extruded from xposvec
					if (tmpx > 0) rev = 0;
					else rev = 1;
					for (j=1;j<=abs(tmpx);j++) growScreenX,cc,reverse=rev;
					xpos(cc) -= floor(xpos(cc));
				}
				
				deltay = screens(cc).deltay;
				tmpy = floor(ypos(cc)+deltay) - floor(ypos(cc));
				ypos(cc) += deltay;
				if (abs(tmpy) > 0) { // we jumped of at least of 1 pixel so
					// we need to extrude the phase screen and remove the number
					// of pixels we extruded from xposvec
					if (tmpy > 0) rev = 0;
					else rev = 1;
					for (j=1;j<=abs(tmpy);j++) growScreenY,cc,reverse=rev;
					ypos(cc) -= floor(ypos(cc));
				}
			} 
		} else {
			xpos = xposvec * 0.;
			ypos = yposvec * 0.;
		}
	}
	//xposvec+=0.5;
	//yposvec+=0.5;
	
	// here we have a branch to be able to process wfs and targets with the same
	// subroutine, this one.
	if (type == "wfs") {
		// mod 2011jan19 w/ DG to get rid of screens above LGS
		if (wfs(nn).gsalt>0) nscreens = long(sum(*atm.layeralt < wfs(nn).gsalt));
		// stuff xshifts with fractionnal offsets, add xposvec for each screen
		xshifts = wfsxposcub(,,nn)+(xposvec+xpos)(-,);
		yshifts = wfsyposcub(,,nn)+(yposvec+ypos)(-,);
		//xshifts = wfsxposcub(,,nn)+(xposvec(1,)+xposvec(iter,)-long(xposvec(iter,)))(-,);
		//yshifts = wfsyposcub(,,nn)+(yposvec(1,)+yposvec(iter,)-long(yposvec(iter,)))(-,);
	} else if ( type == "target") {
		// stuff xshifts with fractionnal offsets, add xposvec for each screen
		xshifts = gsxposcub(,,nn)+(xposvec+xpos)(-,);
		yshifts = gsyposcub(,,nn)+(yposvec+ypos)(-,);
	}
	
	ishifts = int(xshifts);  xshifts = xshifts - ishifts;
	jshifts = int(yshifts);  yshifts = yshifts - jshifts;
	
	err = _get2dPhase(&(pscreens* currentScreenNorm(-,-,)),psnx,psny,nscreens,
					  &sphase,_n,_n,&ishifts,&float(xshifts),&jshifts,&float(yshifts));
	
	if (err != 0) {error,"Error in get_turb_phase";}
	
	bphase(_n1:_n2,_n1:_n2) = sphase;
	
	return bphase;
}

//----------------------------------------------------

func growScreenX(ns,reverse=)
{
	extern pscreens;
	
	if (reverse == []) reverse = 0;
	
	A = *(screens(ns).A);
	B = *(screens(ns).B);
	istencil = *(screens(ns).istencil);
	r0 = screens(ns).r0;
	
	if (reverse) scr = pscreens(,,ns)(::-1,);
	else scr = pscreens(,,ns);
	scr = extrude(scr,r0,A,B,istencil);
	if (reverse) pscreens(,,ns) = scr(::-1,);
	else pscreens(,,ns) = scr;
}

func growScreenY(ns,reverse=)
{
	extern pscreens;
	
	if (reverse == []) reverse = 0;
	
	A = *(screens(ns).A);
	B = *(screens(ns).B);
	istencil = *(screens(ns).istencil);
	r0 = screens(ns).r0;
	
	pscreens(,,ns) = transpose(pscreens(,,ns));
	if (reverse) scr = pscreens(,,ns)(::-1,);
	else scr = pscreens(,,ns);
	scr = extrude(scr,r0,A,B,istencil);
	if (reverse) pscreens(,,ns) = scr(::-1,);
	else pscreens(,,ns) = scr;
	pscreens(,,ns) = transpose(pscreens(,,ns));
}

//----------------------------------------------------

func make_kl_dm(nm,&def,disp=)

  /* DOCUMENT function make_kl_dm,dm_number,ActIF,disp=
   */
{
  require,"yaokl.i";
  dim   = dm(nm)._n2-dm(nm)._n1+1;
  //  nkllow = dm(nm).nklfiltered;
  nkllow = 1;
  nkl = dm(nm).nkl;
  cent  = sim._cent;
  psize = tel.diam/sim.pupildiam;
  gsdist = sqrt((abs(wfs.gspos)^2.)(sum,));
  patchDiam = long(ceil((sim.pupildiam+2*max(gsdist)*
                         4.848e-6*abs(dm(nm).alt)/psize)/2)*2);

  //  prepzernike,dim,patchDiam,sim._cent-dm(nm)._n1+1,sim._cent-dm(nm)._n1+1;

  if (dm(nm).alt==0) {
    // enforce pupil to be system pupil.
    i1 = sim._size/2 - sim.pupildiam/2+1;
    i2 = sim._size/2 + sim.pupildiam/2;
    outpup = ipupil(i1:i2,i1:i2);
    //patchDiam should be good.
    write,format="KL: PatchDiam = %d, sim.pupildiam=%d\n",patchDiam,sim.pupildiam;
    cobs = tel.cobs;
  } else {
    outpup = [];
    cobs = 0.;
    patchDiam += 2; // margin
  }

  //kl = float(make_kl(nkl,patchDiam,varkl,outbas,outpup,oc=cobs,nr=128));

  // order them in a similar order as zernike:
	//vt = fits_read("vt.fits");
	//kl = order_kls(kl,patchDiam,upto=20);
	write, "HACKING of KLs";
	kl	= fits_read("new_kl_dm.fits")(.., :nkl) ;
	//write, "HACKING of KLs";

  def = array(float,dim,dim,nkl-nkllow+1);

  n1 = dim/2-patchDiam/2+1;
  n2 = n1+patchDiam-1;
  for (i=nkllow;i<=nkl;i++) {
    def(n1:n2,n1:n2,i-nkllow+1) = kl(,,i);
    if (disp == 1) {fma; pli,def(,,i-nkllow+1);}
  }
  if (sim.verbose>=1) {write,format="Number of KL :%d\n",nkl-nkllow+1;}

  // the KL are normalized so that rms over surface = 1 unit
  // meaning the TT go from -2 to 2.

  // we'll use the same normalization as zernike, just for consistency

  // current: over tel.diam, we have 4 units, and a rms of 1 we want a
  // rms of 957nm = 0.957microns

  //def = def * 0.957f * float(dm(nm).unitpervolt);

  dm(nm)._nact = (dimsof(def))(4);
  dm(nm)._def = &def;

  return def;
}
