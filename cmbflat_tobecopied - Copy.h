
        


	/************************************************************************************************************************************/
	/*		All the calculation of the CMBFLAT is done in this function. This is the main function of this program. Some of the Cl will  	*/
	/*		be calculated in  this program other Cl's will be manupulated by fitting an spline															*/
	/************************************************************************************************************************************/

                                        
                
	void CMBflat()
	{

 // Tempurary variable
					//	double tauf[ntfmax+1];



		indtau[ntau + ntfmax + 1] = 0;



		int nriend, nk, nk1, nk2, n10;


		double akchange;
		double ak0[nk0 + 1], ak10, t10,  ar; // ak0[] <- 

		int itau, ind; // itau is a flag

		double phi1, akt, ho, xf, b0; //,bo
		int closestknum;
		double tol1; // tolerance vaue

		double zsanpass[4], phi, stpt, ysanpass[3], akdone, akgf, phihk;

		double c[25]; // dvark intigrator communication parameter


		double dk, dk0, dlnk1;


		double de2, xlim, xlmax1, xlmax2, tmin, tmax, h2, a2, b2;
		double ddt1, ddt2, ddp1, ddp2, ddk1, ddk2, ddt, ddp, ddk;
		double ddt1adb, ddt2adb, ddp1adb, ddp2adb, ddk1adb, ddk2adb, ddtadb, ddpadb, ddkadb;
		double ddt1iso, ddt2iso, ddp1iso, ddp2iso, ddk1iso, ddk2iso, ddtiso, ddpiso, ddkiso;

		double xlimmin = 10.0;
		double dtau3, xi, x;
		double ajl0;
		double ckj, cpkj, cckj, ckkkj, ctkkj;
		double dtdt1, dtdt2, dtde1, dtdb1, dtde2, dtdb2, dtdt, dtdb, dtde;
		double ctkj, ctekj, ctbkj, ctckj, apowert;
		double cllo, clhi;
		double ctnorm, xl[lmax];
		double clint, cplint, cclint, ctklint; //,ckklont
		double clintadb, cplintadb, cclintadb;
		double clintiso, cplintiso, cclintiso;
		double clintcross, cplintcross, cclintcross;
		double clinttotal, cplinttotal, cclinttotal;
		double ctlint, ctelint, ctblint, ctclint, ckklint;

		int nstart1, nstop1, nstop1a, nstart2, m2, nstop2, nstop2a;
		int no, no1, nko;
		int nkt, nstps, nstpt;
		int khi, klo;
		int mxx[nstep0 + 1], m1;
		int llo, lhi;


      
                /*
		double dtl[lmax + 1], dtl2[lmax + 1], dtl3[lmax + 1];
		double dtel[lmax + 1], dtel2[lmax + 1], dtel3[lmax + 1];
		double dtbl[lmax + 1], dtbl2[lmax + 1], dtbl3[lmax + 1];
		double st2[nstep0 + 1], ste2[nstep0 + 1];
		double stb2[nstep0 + 1];

		double dl[lmax + 1] = {0}, dl2[lmax + 1] = {0}, dl3[lmax + 1] = {0};
		double dpl[lmax + 1] = {0}, dpl2[lmax + 1] = {0}, dpl3[lmax + 1] = {0};
		double dkl[lmax + 1] = {0}, dkl2[lmax + 1] = {0}, dkl3[lmax + 1] = {0};
		double s2[nstep0 + 1], sp2[nstep0 + 1], sk2[nstep0 + 1];

		double dladb[lmax + 1] = {0}, dl2adb[lmax + 1] = {0}, dl3adb[lmax + 1] = {0};
		double dpladb[lmax + 1] = {0}, dpl2adb[lmax + 1] = {0}, dpl3adb[lmax + 1] = {0};
		double dkladb[lmax + 1] = {0}, dkl2adb[lmax + 1] = {0}, dkl3adb[lmax + 1] = {0};
		double s2adb[nstep0 + 1], sp2adb[nstep0 + 1], sk2adb[nstep0 + 1];

		double dliso[lmax + 1] = {0}, dl2iso[lmax + 1] = {0}, dl3iso[lmax + 1] = {0};
		double dpliso[lmax + 1] = {0}, dpl2iso[lmax + 1] = {0}, dpl3iso[lmax + 1] = {0};
		double dkliso[lmax + 1] = {0}, dkl2iso[lmax + 1] = {0}, dkl3iso[lmax + 1] = {0};
		double s2iso[nstep0 + 1], sp2iso[nstep0 + 1], sk2iso[nstep0 + 1];

                
                
	double d[nk0 + 1][nstep0 + 1];	//	Source terms 												/memory1/
	double dpr[nk0 + 1][nstep0 + 1];  //	Scaler														/memory1/
	double dp[nk0 + 1][nstep0 + 1];   //																	/memory1/
	double dppr[nk0 + 1][nstep0 + 1]; //																	/memory1/
	double dkk[nk0 + 1][nstep0 + 1];  //																	/memory1/
	double dpkr[nk0 + 1][nstep0 + 1]; //																	/memory1/

	double dkpradb[nk0 + 1][nstep0 + 1]; //																	/memory4/
	double dkpriso[nk0 + 1][nstep0 + 1]; //																	/memory4/

	double dadb[nk0 + 1][nstep0 + 1];   //	Adiabatic perturbation									/memory3/
	double dpadb[nk0 + 1][nstep0 + 1];  //																	/memory3/
	double dkkadb[nk0 + 1][nstep0 + 1]; //																	/memory3/
	double diso[nk0 + 1][nstep0 + 1];   //	Isocurvature perturbation								/memory3/
	double dpiso[nk0 + 1][nstep0 + 1];  //																	/memory3/
	double dkkiso[nk0 + 1][nstep0 + 1]; //																	/memory3/

	double dpradb[nk0 + 1][nstep0 + 1];   //	Adiabatic perturbation									/memory3/
	double dppradb[nk0 + 1][nstep0 + 1];  //																	/memory3/
	double dkkpradb[nk0 + 1][nstep0 + 1]; //																	/memory3/
	double dpriso[nk0 + 1][nstep0 + 1];   //	Isocurvature perturbation								/memory3/
	double dppriso[nk0 + 1][nstep0 + 1];  //																	/memory3/
	double dkkpriso[nk0 + 1][nstep0 + 1]; //																	/memory3/

	double dt[nk0 + 1][nstep0 + 1];	//	Source terms												/memory2/
	double dtpr[nk0 + 1][nstep0 + 1];  //	Tensor														/memory2/
	double dte[nk0 + 1][nstep0 + 1];   //																	/memory2/
	double dtepr[nk0 + 1][nstep0 + 1]; //																	/memory2/
	double dtb[nk0 + 1][nstep0 + 1];   //																	/memory2/
	double dtbpr[nk0 + 1][nstep0 + 1]; //																	/memory2/
                
                

    	double ctl[lmax + 1][nnmax + 1];	//
	double ctel[lmax + 1][nnmax + 1];   //
	double ctbl[lmax + 1][nnmax + 1];   //
	double ctcl[lmax + 1][nnmax + 1];   //

  	double clpr[lmax + 1][nnmax + 1] = {{0}};   //
	double cplpr[lmax + 1][nnmax + 1] = {{0}};  //
	double cclpr[lmax + 1][nnmax + 1] = {{0}};  //
	double ckklpr[lmax + 1][nnmax + 1] = {{0}}; //
	double ctklpr[lmax + 1][nnmax + 1] = {{0}}; //
 	
	double clpradb[lmax + 1][nnmax + 1] = {{0}};  //
	double cplpradb[lmax + 1][nnmax + 1] = {{0}}; //
	double cclpradb[lmax + 1][nnmax + 1] = {{0}}; // 
  
	double clpriso[lmax + 1][nnmax + 1] = {{0}};  //
	double cplpriso[lmax + 1][nnmax + 1] = {{0}}; //
	double cclpriso[lmax + 1][nnmax + 1] = {{0}}; //
 
  	double clprcross[lmax + 1][nnmax + 1] = {{0}};  // 
	double cplprcross[lmax + 1][nnmax + 1] = {{0}}; //
	double cclprcross[lmax + 1][nnmax + 1] = {{0}}; //
 
	double clprtotal[lmax + 1][nnmax + 1] = {{0}};  //  
	double cplprtotal[lmax + 1][nnmax + 1] = {{0}}; //  
  	double cclprtotal[lmax + 1][nnmax + 1] = {{0}}; //
  
 	double ctlpr[lmax + 1][nnmax + 1];  //	Tensor Cl's 
	double ctelpr[lmax + 1][nnmax + 1]; //
  	double ctblpr[lmax + 1][nnmax + 1]; //
  	double ctclpr[lmax + 1][nnmax + 1]; //

                 
                dl2[lmax] = 0;
		dl3[lmax] = 0;
		dpl2[lmax] = 0;
		dpl3[lmax] = 0;
		dkl2[lmax] = 0;
		dkl3[lmax] = 0;
                
                */                
                
		numericx Numericx;
		others Others;
		
		
		neutrino Neutrino;
		doubleinflation Doubleinflation;
	
		
                
		nw = nvar0;


                
                
                














					for (int j = 1; j <= Variables.l0; j++)
					{
						dtl[j] = 0.0;
						dtl2[j] = 0.0;
						dtl3[j] = 0.0;
						dtel2[j] = 0.0;
						dtbl2[j] = 0.0;
					}









			}	 // K LOOP ENDS HERE


        
        
        
			/***************************************************************************************************************************/
			/*     Tensor Case																																			*/
			/***************************************************************************************************************************/
			if (Variables.itflag != 0)
			{
				/************************************************************************************************************************/
				/*    Normalization																																		*/
				/************************************************************************************************************************/
				for (int j = 1; j <= Variables.l0; j++)
				{
					ctnorm = Variables.l[j] * (Variables.l[j] - 1.0) * (Variables.l[j] + 1) * (Variables.l[j] + 2);
					for (int in = 1; in <= Variables.nn; in++)
					{
						ctl[j][in] = ctnorm * ctl[j][in] * (Variables.l[j] * (Variables.l[j] + 1)) / 8.0;
						ctel[j][in] = ctel[j][in] * (Variables.l[j] * (Variables.l[j] + 1)) / 8.0;
						ctbl[j][in] = ctbl[j][in] * (Variables.l[j] * (Variables.l[j] + 1)) / 8.0;
						ctcl[j][in] = sqrt(ctnorm) * ctcl[j][in] * Variables.l[j] * (Variables.l[j] + 1) / 8.0;
					}
				}

				/************************************************************************************************************************/
				/*		Making the interpolation tables to get other l-values.																				*/
				/************************************************************************************************************************/
				cllo = 1.0e30;
				clhi = 1.0e30;

				for (int in = 1; in <= Variables.nn; in++)
				{
					for (int jj = 1; jj <= lmax + 1; jj++)
					{
						bpctl[jj] = ctl[jj][in];
						bpctel[jj] = ctel[jj][in];
						bpctbl[jj] = ctbl[jj][in];
						bpctcl[jj] = ctcl[jj][in];
					}

					Numericx.spline(xl, bpctl, Variables.l0, cllo, clhi, bpctlpr);
					Numericx.spline(xl, bpctel, Variables.l0, cllo, clhi, bpctelpr);
					Numericx.spline(xl, bpctbl, Variables.l0, cllo, clhi, bpctblpr);
					Numericx.spline(xl, bpctcl, Variables.l0, cllo, clhi, bpctclpr);

					for (int jj = 1; jj <= lmax + 1; jj++)
					{
						ctlpr[jj][in] = bpctlpr[jj];
						ctelpr[jj][in] = bpctelpr[jj];
						ctblpr[jj][in] = bpctblpr[jj];
						ctclpr[jj][in] = bpctclpr[jj];
					}
				}
			} // IF Variables.itflag!=1 ENDS HERE

			/************************************************************************************************************************/
			/*     Calculating Cls for every l.																													*/
			/************************************************************************************************************************/

			// (P20) : Interpolate the Cls for all the l's
			for (int in = 1; in <= Variables.nn; in++)
			{
				llo = 1;
				for (int il = 2; il <= Variables.l[Variables.l0]; il++)
				{
					xi = il;

					if ((xi > xl[llo + 1]) && (llo < Variables.l0))
						llo = llo + 1;

					lhi = llo + 1;
					ho = xl[lhi] - xl[llo];
					a0 = (xl[lhi] - xi) / ho;
					b0 = (xi - xl[llo]) / ho;

					if (Variables.itflag != 2)
					{
						if (Variables.twofield == 0)
						{
							clint = a0 * cl[llo][in] + b0 * cl[lhi][in] + ((pow(a0, 3) - a0) * clpr[llo][in] + (pow(b0, 3) - b0) * clpr[lhi][in]) * ho * ho / 6.0;
							cplint = a0 * cpl[llo][in] + b0 * cpl[lhi][in] + ((pow(a0, 3) - a0) * cplpr[llo][in] + (pow(b0, 3) - b0) * cplpr[lhi][in]) * ho * ho / 6.0;
							cclint = a0 * ccl[llo][in] + b0 * ccl[lhi][in] + ((pow(a0, 3) - a0) * cclpr[llo][in] + (pow(b0, 3) - b0) * cclpr[lhi][in]) * ho * ho / 6.0;
							ckklint = a0 * ckkl[llo][in] + b0 * ckkl[lhi][in] + ((pow(a0, 3) - a0) * ckklpr[llo][in] + (pow(b0, 3) - b0) * ckklpr[lhi][in]) * ho * ho / 6.0;
							ctklint = a0 * ctkl[llo][in] + b0 * ctkl[lhi][in] + ((pow(a0, 3) - a0) * ctklpr[llo][in] + (pow(b0, 3) - b0) * ctklpr[lhi][in]) * ho * ho / 6.0;
						}
						if (Variables.twofield == 1)
						{
							clintadb = a0 * cladb[llo][in] + b0 * cladb[lhi][in] + ((pow(a0, 3) - a0) * clpradb[llo][in] + (pow(b0, 3) - b0) * clpradb[lhi][in]) * ho * ho / 6.0;
							cplintadb = a0 * cpladb[llo][in] + b0 * cpladb[lhi][in] + ((pow(a0, 3) - a0) * cplpradb[llo][in] + (pow(b0, 3) - b0) * cplpradb[lhi][in]) * ho * ho / 6.0;
							cclintadb = a0 * ccladb[llo][in] + b0 * ccladb[lhi][in] + ((pow(a0, 3) - a0) * cclpradb[llo][in] + (pow(b0, 3) - b0) * cclpradb[lhi][in]) * ho * ho / 6.0;

							clintiso = a0 * cliso[llo][in] + b0 * cliso[lhi][in] + ((pow(a0, 3) - a0) * clpriso[llo][in] + (pow(b0, 3) - b0) * clpriso[lhi][in]) * ho * ho / 6.0;
							cplintiso = a0 * cpliso[llo][in] + b0 * cpliso[lhi][in] + ((pow(a0, 3) - a0) * cplpriso[llo][in] + (pow(b0, 3) - b0) * cplpriso[lhi][in]) * ho * ho / 6.0;
							cclintiso = a0 * ccliso[llo][in] + b0 * ccliso[lhi][in] + ((pow(a0, 3) - a0) * cclpriso[llo][in] + (pow(b0, 3) - b0) * cclpriso[lhi][in]) * ho * ho / 6.0;

							clintcross = a0 * clcross[llo][in] + b0 * clcross[lhi][in] + ((pow(a0, 3) - a0) * clprcross[llo][in] + (pow(b0, 3) - b0) * clprcross[lhi][in]) * ho * ho / 6.0;
							cplintcross = a0 * cplcross[llo][in] + b0 * cplcross[lhi][in] + ((pow(a0, 3) - a0) * cplprcross[llo][in] + (pow(b0, 3) - b0) * cplprcross[lhi][in]) * ho * ho / 6.0;
							cclintcross = a0 * cclcross[llo][in] + b0 * cclcross[lhi][in] + ((pow(a0, 3) - a0) * cclprcross[llo][in] + (pow(b0, 3) - b0) * cclprcross[lhi][in]) * ho * ho / 6.0;

							clinttotal = a0 * cltotal[llo][in] + b0 * cltotal[lhi][in] + ((pow(a0, 3) - a0) * clprtotal[llo][in] + (pow(b0, 3) - b0) * clprtotal[lhi][in]) * ho * ho / 6.0;
							cplinttotal = a0 * cpltotal[llo][in] + b0 * cpltotal[lhi][in] + ((pow(a0, 3) - a0) * cplprtotal[llo][in] + (pow(b0, 3) - b0) * cplprtotal[lhi][in]) * ho * ho / 6.0;
							cclinttotal = a0 * ccltotal[llo][in] + b0 * ccltotal[lhi][in] + ((pow(a0, 3) - a0) * cclprtotal[llo][in] + (pow(b0, 3) - b0) * cclprtotal[lhi][in]) * ho * ho / 6.0;
						}
					}

					if (Variables.itflag != 0)
					{
						ctlint = a0 * ctl[llo][in] + b0 * ctl[lhi][in] + ((pow(a0, 3) - a0) * ctlpr[llo][in] + (pow(b0, 3) - b0) * ctlpr[lhi][in]) * ho * ho / 6.0;
						ctelint = a0 * ctel[llo][in] + b0 * ctel[lhi][in] + ((pow(a0, 3) - a0) * ctelpr[llo][in] + (pow(b0, 3) - b0) * ctelpr[lhi][in]) * ho * ho / 6.0;
						ctblint = a0 * ctbl[llo][in] + b0 * ctbl[lhi][in] + ((pow(a0, 3) - a0) * ctblpr[llo][in] + (pow(b0, 3) - b0) * ctblpr[lhi][in]) * ho * ho / 6.0;
						ctclint = a0 * ctcl[llo][in] + b0 * ctcl[lhi][in] + ((pow(a0, 3) - a0) * ctclpr[llo][in] + (pow(b0, 3) - b0) * ctclpr[lhi][in]) * ho * ho / 6.0;
					}

					if (Variables.itflag != 2)
					{
						if (Variables.twofield == 0)
						{
							GlobalArray.clts[il][in] = clint / cl[1][in];
							GlobalArray.cles[il][in] = cplint / cl[1][in];
							GlobalArray.clcs[il][in] = cclint / cl[1][in];
							GlobalArray.clkk[il][in] = ckklint / cl[1][in];
							GlobalArray.cltk[il][in] = ctklint / cl[1][in];
						}
						if (Variables.twofield == 1)
						{
							/*
		         	cltsadb[il][in]=clintadb/cladb[1][in];
   	           	clesadb[il][in]=cplintadb/cladb[1][in];
      	        	clcsadb[il][in]=cclintadb/cladb[1][in];

		         	cltsiso[il][in]=clintiso/cliso[1][in];
   	           	clesiso[il][in]=cplintiso/cliso[1][in];
      	        	clcsiso[il][in]=cclintiso/cliso[1][in];

		         	cltscross[il][in]=clintcross/clcross[1][in];
   	           	clescross[il][in]=cplintcross/clcross[1][in];
      	        	clcscross[il][in]=cclintcross/clcross[1][in];
						*/
							cltsadb[il][in] = clintadb / cltotal[1][in];
							clesadb[il][in] = cplintadb / cltotal[1][in];
							clcsadb[il][in] = cclintadb / cltotal[1][in];

							cltsiso[il][in] = clintiso / cltotal[1][in];
							clesiso[il][in] = cplintiso / cltotal[1][in];
							clcsiso[il][in] = cclintiso / cltotal[1][in];

							cltscross[il][in] = clintcross / cltotal[1][in];
							clescross[il][in] = cplintcross / cltotal[1][in];
							clcscross[il][in] = cclintcross / cltotal[1][in];

							cltstotal[il][in] = clinttotal / cltotal[1][in];
							clestotal[il][in] = cplinttotal / cltotal[1][in];
							clcstotal[il][in] = cclinttotal / cltotal[1][in];
						}
					}

					if (Variables.itflag != 0)
					{
						cltt[il][in] = ctlint / ctl[1][in];
						clet[il][in] = ctelint / ctl[1][in];
						clbt[il][in] = ctblint / ctl[1][in];
						clct[il][in] = ctclint / ctl[1][in];
					}
				}

				if (Variables.itflag != 2)
				{
					if (Variables.twofield == 0)
						GlobalArray.clts[1][in] = cl[1][in];
					else if (Variables.twofield == 1)
					{
						cltsadb[1][in] = cltotal[1][in];
						cltsiso[1][in] = cltotal[1][in];
						cltscross[1][in] = cltotal[1][in];
						cltstotal[1][in] = cltotal[1][in];
					}
				}

				if (Variables.itflag != 0)
					cltt[1][in] = ctl[1][in];
			}

 		}
	}
        }
        
	void COBEnormalizetwofield()
	{

/*		const double fourpi = 4.0 * 3.14159265;
		double xlog10, h, xlnh, hc, curv, r;
		double c10, d1, d2, d3, d4, d5, d6, d7, x1, x2, x3, x4, x5, x6, x7;
		double sy, s, sx, sxy, sxx, delt, d1pr, d1ppr;

		double d2norm; // Needed for lensing calculation

		xlog10 = log(10.0);
		h = Variables.h0 / 100.0;
		xlnh = log(h);

		// Curvature radius
		if (fabs(Variables.omegak) > 1.0e-3)
		{
			hc = 2.998e5 / Variables.h0;			// h*c. c = 2.338*e5 km/s. h = km/sec/pc
			curv = -Variables.omegak / (hc * hc); // Gaussian curvature
			r = 1.0 / sqrt(fabs(curv)); // Radious of curvature
		}

		// COBE normalization
		// fit the spectrum to a quadratic around C_10 with equal weights in logl
		/*
	for(int in=1;in<=Variables.nn;in++)
	{
		c10=cltsadb[10][in];
		printf("cl10 = %e cltt[10][in] = %e\n",c10,cltsadb[10][in]);
      d1=(cltsadb[Variables.l[2]][in])/c10-1.0;
      d2=(cltsadb[Variables.l[3]][in])/c10-1.0;
      d3=(cltsadb[Variables.l[5]][in])/c10-1.0;
      d4=(cltsadb[Variables.l[7]][in])/c10-1.0;
      d5=(cltsadb[Variables.l[10]][in])/c10-1.0;
      d6=(cltsadb[Variables.l[11]][in])/c10-1.0;
      d7=(cltsadb[Variables.l[12]][in])/c10-1.0;

      x1=log(1.0*Variables.l[2])/xlog10-1.0;
      x2=log(1.0*Variables.l[3])/xlog10-1.0;
      x3=log(1.0*Variables.l[5])/xlog10-1.0;
      x4=log(1.0*Variables.l[7])/xlog10-1.0;
      x5=log(1.0*Variables.l[10])/xlog10-1.0;
      x6=log(1.0*Variables.l[11])/xlog10-1.0;
      x7=log(1.0*Variables.l[12])/xlog10-1.0;

		// Quadratic List Square fit
      sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7;
      s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7;
      sx=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4+x5*x5*x5+x6*x6*x6+x7*x7*x7;
      sxy=x1*x1*d1+x2*x2*d2+x3*x3*d3+x4*x4*d4+x5*x5*d5+x6*x6*d6+x7*x7*d7;
      sxx=pow(x1,4)+pow(x2,4)+pow(x3,4)+pow(x4,4)+pow(x5,4)+pow(x6,4)+pow(x7,4);

      delt=s*sxx-sx*sx;
      d1pr=(sxx*sy-sx*sxy)/delt;
      d1ppr=2.0*(s*sxy-sx*sy)/delt;

		printf("%e %e %e\n",delt,d1pr,d1ppr);

		// Bunn and White fitting formula
      c10=(0.64575+0.02282*d1pr+0.01391*d1pr*d1pr-0.01819*d1ppr-0.00646*d1pr*d1ppr +0.00103*d1ppr*d1ppr)/c10;

		// C_l normalization; output l(l+1)C_l/twopi

		c10=(c10*2.2e-9)/fourpi;
		printf("cl10 = %e \n",c10);
		for(int il=2;il<=Variables.l[Variables.l0];il++)
		{
			printf("clts = %e \n",cltsadb[il][in]);
			cltsadb[il][in]=cltsadb[il][in]*c10;
	      clesadb[il][in]=clesadb[il][in]*c10;
   	   clcsadb[il][in]=clcsadb[il][in]*c10;
   	   }
		}
	*/
		/*
	for(int in=1;in<=Variables.nn;in++)
	{
		c10=cltsiso[10][in];
		printf("cl10 = %e cltt[10][in] = %e\n",c10,clts[10][in]);
      d1=(cltsiso[Variables.l[2]][in])/c10-1.0;
      d2=(cltsiso[Variables.l[3]][in])/c10-1.0;
      d3=(cltsiso[Variables.l[5]][in])/c10-1.0;
      d4=(cltsiso[Variables.l[7]][in])/c10-1.0;
      d5=(cltsiso[Variables.l[10]][in])/c10-1.0;
      d6=(cltsiso[Variables.l[11]][in])/c10-1.0;
      d7=(cltsiso[Variables.l[12]][in])/c10-1.0;

      x1=log(1.0*Variables.l[2])/xlog10-1.0;
      x2=log(1.0*Variables.l[3])/xlog10-1.0;
      x3=log(1.0*Variables.l[5])/xlog10-1.0;
      x4=log(1.0*Variables.l[7])/xlog10-1.0;
      x5=log(1.0*Variables.l[10])/xlog10-1.0;
      x6=log(1.0*Variables.l[11])/xlog10-1.0;
      x7=log(1.0*Variables.l[12])/xlog10-1.0;

		// Quadratic List Square fit
      sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7;
      s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7;
      sx=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4+x5*x5*x5+x6*x6*x6+x7*x7*x7;
      sxy=x1*x1*d1+x2*x2*d2+x3*x3*d3+x4*x4*d4+x5*x5*d5+x6*x6*d6+x7*x7*d7;
      sxx=pow(x1,4)+pow(x2,4)+pow(x3,4)+pow(x4,4)+pow(x5,4)+pow(x6,4)+pow(x7,4);

      delt=s*sxx-sx*sx;
      d1pr=(sxx*sy-sx*sxy)/delt;
      d1ppr=2.0*(s*sxy-sx*sy)/delt;

		printf("%e %e %e\n",delt,d1pr,d1ppr);

		// Bunn and White fitting formula
      c10=(0.64575+0.02282*d1pr+0.01391*d1pr*d1pr-0.01819*d1ppr-0.00646*d1pr*d1ppr +0.00103*d1ppr*d1ppr)/c10;

		// C_l normalization; output l(l+1)C_l/twopi

		c10=(c10*2.2e-9)/fourpi;
		printf("cl10 = %e \n",c10);
		for(int il=2;il<=Variables.l[Variables.l0];il++)
		{
			printf("clts = %e \n",cltsiso[il][in]);
			cltsiso[il][in]=cltsiso[il][in]*c10;
	      clesiso[il][in]=clesiso[il][in]*c10;
   	   clcsiso[il][in]=clcsiso[il][in]*c10;
   	   }
		}
	*/
		/*
	for(int in=1;in<=Variables.nn;in++)
	{
		c10=cltscross[10][in];
		printf("cl10 = %e cltt[10][in] = %e\n",c10,clts[10][in]);
      d1=(cltscross[Variables.l[2]][in])/c10-1.0;
      d2=(cltscross[Variables.l[3]][in])/c10-1.0;
      d3=(cltscross[Variables.l[5]][in])/c10-1.0;
      d4=(cltscross[Variables.l[7]][in])/c10-1.0;
      d5=(cltscross[Variables.l[10]][in])/c10-1.0;
      d6=(cltscross[Variables.l[11]][in])/c10-1.0;
      d7=(cltscross[Variables.l[12]][in])/c10-1.0;

      x1=log(1.0*Variables.l[2])/xlog10-1.0;
      x2=log(1.0*Variables.l[3])/xlog10-1.0;
      x3=log(1.0*Variables.l[5])/xlog10-1.0;
      x4=log(1.0*Variables.l[7])/xlog10-1.0;
      x5=log(1.0*Variables.l[10])/xlog10-1.0;
      x6=log(1.0*Variables.l[11])/xlog10-1.0;
      x7=log(1.0*Variables.l[12])/xlog10-1.0;

		// Quadratic List Square fit
      sy=x1*d1+x2*d2+x3*d3+x4*d4+x5*d5+x6*d6+x7*d7;
      s=x1*x1+x2*x2+x3*x3+x4*x4+x5*x5+x6*x6+x7*x7;
      sx=x1*x1*x1+x2*x2*x2+x3*x3*x3+x4*x4*x4+x5*x5*x5+x6*x6*x6+x7*x7*x7;
      sxy=x1*x1*d1+x2*x2*d2+x3*x3*d3+x4*x4*d4+x5*x5*d5+x6*x6*d6+x7*x7*d7;
      sxx=pow(x1,4)+pow(x2,4)+pow(x3,4)+pow(x4,4)+pow(x5,4)+pow(x6,4)+pow(x7,4);

      delt=s*sxx-sx*sx;
      d1pr=(sxx*sy-sx*sxy)/delt;
      d1ppr=2.0*(s*sxy-sx*sy)/delt;

		printf("%e %e %e\n",delt,d1pr,d1ppr);

		// Bunn and White fitting formula
      c10=(0.64575+0.02282*d1pr+0.01391*d1pr*d1pr-0.01819*d1ppr-0.00646*d1pr*d1ppr +0.00103*d1ppr*d1ppr)/c10;

		// C_l normalization; output l(l+1)C_l/twopi

		c10=(c10*2.2e-9)/fourpi;
		printf("cl10 = %e \n",c10);
		for(int il=2;il<=Variables.l[Variables.l0];il++)
		{
			printf("clts = %e \n",clts[il][in]);
			cltscross[il][in]=cltscross[il][in]*c10;
	      clescross[il][in]=clescross[il][in]*c10;
   	   clcscross[il][in]=clcscross[il][in]*c10;
   	   }
		}
	*/
/*		for (int in = 1; in <= Variables.nn; in++)
		{
			c10 = cltstotal[10][in];
			printf("cl10 = %e cltt[10][in] = %e\n", c10, GlobalArray.clts[10][in]);
			d1 = (cltstotal[Variables.l[2]][in]) / c10 - 1.0;
			d2 = (cltstotal[Variables.l[3]][in]) / c10 - 1.0;
			d3 = (cltstotal[Variables.l[5]][in]) / c10 - 1.0;
			d4 = (cltstotal[Variables.l[7]][in]) / c10 - 1.0;
			d5 = (cltstotal[Variables.l[10]][in]) / c10 - 1.0;
			d6 = (cltstotal[Variables.l[11]][in]) / c10 - 1.0;
			d7 = (cltstotal[Variables.l[12]][in]) / c10 - 1.0;

			x1 = log(1.0 * Variables.l[2]) / xlog10 - 1.0;
			x2 = log(1.0 * Variables.l[3]) / xlog10 - 1.0;
			x3 = log(1.0 * Variables.l[5]) / xlog10 - 1.0;
			x4 = log(1.0 * Variables.l[7]) / xlog10 - 1.0;
			x5 = log(1.0 * Variables.l[10]) / xlog10 - 1.0;
			x6 = log(1.0 * Variables.l[11]) / xlog10 - 1.0;
			x7 = log(1.0 * Variables.l[12]) / xlog10 - 1.0;

			// Quadratic List Square fit
			sy = x1 * d1 + x2 * d2 + x3 * d3 + x4 * d4 + x5 * d5 + x6 * d6 + x7 * d7;
			s = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5 + x6 * x6 + x7 * x7;
			sx = x1 * x1 * x1 + x2 * x2 * x2 + x3 * x3 * x3 + x4 * x4 * x4 + x5 * x5 * x5 + x6 * x6 * x6 + x7 * x7 * x7;
			sxy = x1 * x1 * d1 + x2 * x2 * d2 + x3 * x3 * d3 + x4 * x4 * d4 + x5 * x5 * d5 + x6 * x6 * d6 + x7 * x7 * d7;
			sxx = pow(x1, 4) + pow(x2, 4) + pow(x3, 4) + pow(x4, 4) + pow(x5, 4) + pow(x6, 4) + pow(x7, 4);

			delt = s * sxx - sx * sx;
			d1pr = (sxx * sy - sx * sxy) / delt;
			d1ppr = 2.0 * (s * sxy - sx * sy) / delt;

			printf("%e %e %e\n", delt, d1pr, d1ppr);

			// Bunn and White fitting formula
			c10 = (0.64575 + 0.02282 * d1pr + 0.01391 * d1pr * d1pr - 0.01819 * d1ppr - 0.00646 * d1pr * d1ppr + 0.00103 * d1ppr * d1ppr) / c10;

			// C_l normalization; output l(l+1)C_l/twopi

			c10 = (c10 * 2.2e-9) / fourpi;
			printf("cl10 = %e \n", c10);
			for (int il = 2; il <= Variables.l[Variables.l0]; il++)
			{
				printf("clts = %e \n", cltstotal[il][in]);

				cltsadb[il][in] = (1 / 9.0) * cltsadb[il][in] * c10;
				clesadb[il][in] = (1 / 9.0) * clesadb[il][in] * c10;
				clcsadb[il][in] = (1 / 9.0) * clcsadb[il][in] * c10;

				cltsiso[il][in] = (4.0 / 25) * cltsiso[il][in] * c10;
				clesiso[il][in] = (4.0 / 25) * clesiso[il][in] * c10;
				clcsiso[il][in] = (4.0 / 25) * clcsiso[il][in] * c10;

				cltscross[il][in] = -(4.0 / 15) * cltscross[il][in] * c10;
				clescross[il][in] = -(4.0 / 15) * clescross[il][in] * c10;
				clcscross[il][in] = -(4.0 / 15) * clcscross[il][in] * c10;

				cltstotal[il][in] = cltstotal[il][in] * c10;
				clestotal[il][in] = clestotal[il][in] * c10;
				clcstotal[il][in] = clcstotal[il][in] * c10;
			}
		} */
	}

	/*********************************************************************************/
	/* Normalization is requiredfor matching the data with the observed data set     */
	/* Otherwise different theory will give different amplitude of the spectrum      */
	/* and it would be difficult to find out the accuricy of the theory              */
	/*                                                                               */
	/* This part is based on the paper: 													      */
	/*      The Four-Year COBE Normalization and Large-Scale Structure			      */
	/*      Emory F. Bunn  and  Martin White												      */
	/*********************************************************************************/

	void COBEnormalize()
	{
/*
		const double fourpi = 4.0 * 3.14159265;
		double xlog10, h, xlnh, hc, curv, r;
		double c10, d1, d2, d3, d4, d5, d6, d7, x1, x2, x3, x4, x5, x6, x7;
		double sy, s, sx, sxy, sxx, delt, d1pr, d1ppr;

		double d2norm; // Needed for lensing calculation

		xlog10 = log(10.0);
		h = Variables.h0 / 100.0;
		xlnh = log(h);

		// Curvature radius
		if (fabs(Variables.omegak) > 1.0e-3)
		{
			hc = 2.998e5 / Variables.h0;			// h*c. c = 2.338*e5 km/s. h = km/sec/pc
			curv = -Variables.omegak / (hc * hc); // Gaussian curvature
			r = 1.0 / sqrt(fabs(curv)); // Radious of curvature
		}

		// Ensuring scalar to tensor ratio
		if (Variables.itflag == 1)
			for (int in = 1; in <= Variables.nn; in++)
				for (int il = 2; il <= Variables.l[Variables.l0]; il++)
				{
					printf("rat =  %e\n", Variables.rat[il]);
					cltt[il][in] = cltt[il][in] * Variables.rat[in];
					clet[il][in] = clet[il][in] * Variables.rat[in];
					clbt[il][in] = clbt[il][in] * Variables.rat[in];
					clct[il][in] = clct[il][in] * Variables.rat[in];
					GlobalArray.clkk[il][in] = GlobalArray.clkk[il][in] * Variables.rat[in];
					GlobalArray.cltk[il][in] = GlobalArray.cltk[il][in] * Variables.rat[in];
				}

		// COBE normalization
		// fit the spectrum to a quadratic around C_10 with equal weights in logl

		for (int in = 1; in <= Variables.nn; in++)
		{
			c10 = GlobalArray.clts[10][in] + cltt[Variables.l[2]][in];
			printf("cl10 = %e cltt[10][in] = %e\n", c10, GlobalArray.clts[10][in]);
			d1 = (GlobalArray.clts[Variables.l[2]][in] + cltt[Variables.l[2]][in]) / c10 - 1.0;
			d2 = (GlobalArray.clts[Variables.l[3]][in] + cltt[Variables.l[3]][in]) / c10 - 1.0;
			d3 = (GlobalArray.clts[Variables.l[5]][in] + cltt[Variables.l[5]][in]) / c10 - 1.0;
			d4 = (GlobalArray.clts[Variables.l[7]][in] + cltt[Variables.l[7]][in]) / c10 - 1.0;
			d5 = (GlobalArray.clts[Variables.l[10]][in] + cltt[Variables.l[10]][in]) / c10 - 1.0;
			d6 = (GlobalArray.clts[Variables.l[11]][in] + cltt[Variables.l[11]][in]) / c10 - 1.0;
			d7 = (GlobalArray.clts[Variables.l[12]][in] + cltt[Variables.l[12]][in]) / c10 - 1.0;

			x1 = log(1.0 * Variables.l[2]) / xlog10 - 1.0;
			x2 = log(1.0 * Variables.l[3]) / xlog10 - 1.0;
			x3 = log(1.0 * Variables.l[5]) / xlog10 - 1.0;
			x4 = log(1.0 * Variables.l[7]) / xlog10 - 1.0;
			x5 = log(1.0 * Variables.l[10]) / xlog10 - 1.0;
			x6 = log(1.0 * Variables.l[11]) / xlog10 - 1.0;
			x7 = log(1.0 * Variables.l[12]) / xlog10 - 1.0;

			// Quadratic List Square fit
			sy = x1 * d1 + x2 * d2 + x3 * d3 + x4 * d4 + x5 * d5 + x6 * d6 + x7 * d7;
			s = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5 + x6 * x6 + x7 * x7;
			sx = x1 * x1 * x1 + x2 * x2 * x2 + x3 * x3 * x3 + x4 * x4 * x4 + x5 * x5 * x5 + x6 * x6 * x6 + x7 * x7 * x7;
			sxy = x1 * x1 * d1 + x2 * x2 * d2 + x3 * x3 * d3 + x4 * x4 * d4 + x5 * x5 * d5 + x6 * x6 * d6 + x7 * x7 * d7;
			sxx = pow(x1, 4) + pow(x2, 4) + pow(x3, 4) + pow(x4, 4) + pow(x5, 4) + pow(x6, 4) + pow(x7, 4);

			delt = s * sxx - sx * sx;
			d1pr = (sxx * sy - sx * sxy) / delt;
			d1ppr = 2.0 * (s * sxy - sx * sy) / delt;

			printf("%e %e %e\n", delt, d1pr, d1ppr);

			// Bunn and White fitting formula
			c10 = (0.64575 + 0.02282 * d1pr + 0.01391 * d1pr * d1pr - 0.01819 * d1ppr - 0.00646 * d1pr * d1ppr + 0.00103 * d1ppr * d1ppr) / c10;

			// Calculations for lensing

			d2norm = c10 * (1.1e-9) / GlobalArray.clts[1][in] / (2 * 3.14159265);

			if (Variables.ict == 2)
			{
				for (int itf = 1; itf <= Variables.ntf; itf++)
				{
					d2norm = d2norm * exp(xlnh * 4.0);
					Variables.anorm[itf][in] = d2norm / ((4 * 3.14159265));
				}
			}

			// C_l normalization; output l(l+1)C_l/twopi

			c10 = (c10 * 2.2e-9) / fourpi;
			printf("cl10 = %e \n", c10);
			for (int il = 2; il <= Variables.l[Variables.l0]; il++)
			{
				printf("clts = %e \n", GlobalArray.clts[il][in]);
				GlobalArray.clts[il][in] = GlobalArray.clts[il][in] * c10;
				GlobalArray.cles[il][in] = GlobalArray.cles[il][in] * c10;
				GlobalArray.clcs[il][in] = GlobalArray.clcs[il][in] * c10;
				cltt[il][in] = cltt[il][in] * c10;
				clet[il][in] = clet[il][in] * c10;
				clbt[il][in] = clbt[il][in] * c10;
				clct[il][in] = clct[il][in] * c10;
				GlobalArray.clkk[il][in] = GlobalArray.clkk[il][in] * c10;
				GlobalArray.cltk[il][in] = GlobalArray.cltk[il][in] * c10;
			}
		}

   */  
 	}

	//void CMBflat() {
	//	cmbflat CMBflat;
	//	CMBflat.cmbflat();
	//}        
     