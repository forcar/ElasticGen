package org.clas.lib;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.jlab.clas.physics.Particle;

public class Elasgen {
	
	MoTsai elib = new MoTsai();
	GenEvent ge;
    Writer writer = new Writer();
	
	double f,d,tl,tr,vx,vy,vz,be,emn,smn,smx,c;
	int    g_,p_,M,opt_strag_0,opt_strag_1,opt_fiduc,idum_off;
	
	int    dismc[][] = new int[6][10];	
	int    npart, id[]=new int[5], pdgid[]=new int[5];//mc.inc
	double v[][]=new double[5][3], cs[][] = new double [5][3], p[][] = new double [5][4], q[]=new double[5], rm[] = new double [5];//mc.inc
	
	double es,ep,th0,ppx,ppy,ppz,eprot,delta,egam,egamx,egamy,egamz,cstk,phik,q2; //bos_events
	
    double ebeam,epmin0,ps,rs,tl_cm,tl_rl,tr_cm,tr_rl,th_min,th_max,fmcall; //parset
    int    nmax,nprint,elas_inc,id1,id2,id3,iext1; //parset    
    double vertex_x,vertex_y,vertex_z,v_z; //new_parset
    
	double mel,ucos0,ucrng,erng;
	double egdiv,epin0,cdiv1,cdiv2,cdiv3,cdiv4,csrng,crng,elast_inc;
	int    intreg;
	double uek,epk;
	
	double theta_pr, phi_pr, phi_r;
	
	double sigr1,sigr_all,sigr_max=0,sig_tot=0,wmax,cst00,cstmin;
	
	int    ntries,nevent,ntold,mcall,mcall_max;
	long   startTime;

    double  pi = 3.141592654;
    double  mp = 0.93827;
    double rme = mel = 0.00051099895069;
    double alpha = 1./137.01;
    double hbarc = 389.37966;
    double  bfac = 4./3.;
    double hydrogen_rad = 865.;
    double[] del_radcor = new double[1000];
    
    public Elasgen() {

    } 
    
    public void getFile(String file) {
    	readFile(file);
		init();    	
    }
    
    public void run() {
		elast_gen();    	
    }
    
    public void readFile(String file) {
    	assignInputParameters(parseInput(getInput(file)));
    }    
    
    public void assignInputParameters(String[] parts) {    	
        f   = Double.parseDouble(parts[0]);         // fraction of events to sample photon energies from 0 to delta
        d   = Double.parseDouble(parts[1]);         // egam < d considered elastic (delta cut for Mo-Tsai)
        g_  = Integer.parseInt(parts[2]);           // number of simulated particles 2: (e-,p+) 3: (e-,p+,gamma)
        tl  = Double.parseDouble(parts[3]);         // target cell length (cm)
        tr  = Double.parseDouble(parts[4]);         // target cell radius (cm)
        vx  = Double.parseDouble(parts[5]);         // target vertex x (cm)
        vy  = Double.parseDouble(parts[6]);         // target vertex y (cm)
        vz  = Double.parseDouble(parts[7]);         // target vertex z (cm)
        be  = Float.parseFloat(parts[8]);           // beam energy (GeV)
        emn = Double.parseDouble(parts[9]);         // minimum scattered e- energy (GeV)
        smn = Double.parseDouble(parts[10]);        // minimum e- scattering angle (deg)
        smx = Double.parseDouble(parts[11]);        // maximum e- scattering angle (deg)
        p_  = Integer.parseInt(parts[12]);          // 0=include 1=exclude elastic peak
        M   = Integer.parseInt(parts[13]);          // simulate M events
        c   = Double.parseDouble(parts[14]);        // maximum cross section (for MC sampling)
        opt_strag_0 = Integer.parseInt(parts[15]);  // straggling of beam e- (1=on, 0=off)
        opt_strag_1 = Integer.parseInt(parts[16]);  // straggling of scattered e- (1=on, 0=off)
        opt_fiduc   = Integer.parseInt(parts[17]);  // use CLAS6 type fiducial cut (no solenoid) (1=on 0=off)
        idum_off    = Integer.parseInt(parts[18]);  // random number offset (obsolete)
    }
    
    public void init() {
    	
        double s,cst0,a1,a2,a3,ep0,epmin;
        double es_save=0,ep_save=0,th0_save=0;
    	
    	csrng = .04;
        egdiv = f;
        delta = d;
        elast_inc = p_;
        npart = g_;
        
        ge = npart==2 ? new gen_e_p() : new gen_e_p_g();

        cdiv1 	= 0.25;
        cdiv2 	= 0.25;
        cdiv3 	= 0.05;
        cdiv4 	= 0.15;
        cdiv2	= cdiv1+cdiv2;
        cdiv3	= cdiv2+cdiv3;
        cdiv4	= cdiv3+cdiv4; //cdiv4 must be >= 0.95
        
        sigr_all = 0;
        wmax = 4.;
     
        par_setup();   
        
        String timeStamp = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format(new java.util.Date());
        output("Starting ElastGen: "+timeStamp,writer.outWriter,writer.sumWriter);
        
        for (int i=0; i<6; i++) {
        	for (int j=0; j<10; j++) {
        		dismc[i][j]=0;
        	}
        }
        
        cst00  = Math.cos(Math.toRadians(th_min));
        cstmin = Math.cos(Math.toRadians(th_max));
        
        ucos0 	= 1./(1.-cstmin);
        ucrng 	= 1./(1-cst00)-ucos0;
        	      
        th0	= Math.toRadians(th_min);
        	      
        s		= Math.pow(Math.sin(th0/2.),2);
        cst0 	= Math.cos(th0);

        a1	= Math.pow(mel*mel+mp*es,2)+Math.pow(mel*ps*cst0,2);
        a2	= (es+mp)*(mel*mel+mp*es);
        a3	= Math.pow(es+mp,2)-Math.pow(ps*cst0,2);
        ep0	= (a2+Math.sqrt(a2*a2-a1*a3))/a3;

        epmin = (mp*mp+2*mp*es-wmax*wmax)/(4*es*s+2*mp);
        	      
        if (epmin<epmin0) epmin=epmin0;
        	      
        erng	= ep0-epmin; 
        
        System.out.println("Integration range for photon energy = "+erng+" GeV");
        output("Minimum photon energy = "+delta,writer.outWriter);
        output("Maximum photon energy = "+erng, writer.outWriter);
       
        fmcall = 1;
        
        System.out.println("EP0,EPMIN = "+(float)ep0+" "+(float)epmin);
        System.out.println("Finding maximum of radiative cross section integrand");
        
        for (int ntries=1; ntries<50001; ntries++) {
        	double sigr1 = sigma_calc();
        	if (sigr1 > sigr_max) {
        		sigr_max = sigr1;
        		es_save = es;
        		ep_save = ep;
        		th0_save = th0;
        	}
        }
        
        System.out.println("Kinematics of sigr_max:");
        System.out.println("ES="+(float)es_save+" EP="+(float)ep_save+" TH="+(float)Math.toDegrees(th0_save)+" SIG="+(float)sigr_max);
        
        sigr_max = fmcall*sigr_max; 
        
        System.out.println("sigr_max changed to "+(float)sigr_max);        
        
    }

    public void int_out() {
    	
    	double space   = 4*pi * ucrng * 2*pi * erng;
    	double sig_int = (sigr_max/ntries) * nevent * space;
    	double sig_sum = (sig_tot /ntries)          * space;
    	
    	System.out.println(" ");
    	System.out.println(" ntries, nevent, mcall_max: "+ntries+" "+nevent+" "+mcall_max);	
    	System.out.println(" sig_int, sig_sum, sig_tot = "+(float)sig_int+" "+(float)sig_sum+" "+(float)sig_tot);
    	System.out.println(" photon phase space = "+(int)space);
    	System.out.println(" Number of seconds at Lum=1E34 is "+(float)((1.E-4)*nevent/sig_int));
    	System.out.println(" Elapsed CPU time = "+(float)cpuTime()+" seconds");
    	
    	output(" ntries,nevent, mcall_max: "+ntries+" "+nevent+" "+mcall_max,writer.outWriter);
    	output(" sig_int = "+sig_int+" sig_sum= "+sig_sum,                   writer.outWriter);
    	output(" Number of seconds at Lum=1E34 is "+(1.E-4)*nevent/sig_int,  writer.outWriter);
    	output(" Elapsed CPU time = "+cpuTime()+" seconds",                  writer.outWriter);
    }

    public void par_setup() {
    	
        tl_cm = tl;
        tr_cm = tr;
        tl_rl = bfac * tl / hydrogen_rad;
        tr_rl = bfac * tr / hydrogen_rad;

        vertex_x = vx;
        vertex_y = vy;
        v_z      = vz;

        ebeam = be;
	    epmin0 = emn;

        System.out.println("Incident electron energy = "+ebeam+" GeV");
        System.out.println("Minimum allowed scattered electron energy = "+epmin0+" GeV");

        es = ebeam;

        // minimum and maximum scattering angles
        th_min = smn;
        th_max = smx;
        	      
        // include the elastic peak in the n-tuple
        elas_inc = p_; //1=yes 0=no
        	      
        // number of events desired in the output file.
        nmax = M;
        nprint = nmax/25;

        sigr_max = c;

        // incident e- momentum
        ps = es-0.5*mel*mel/es;
        rs = ps/es;

        // Calculate a table of internal radiative corrections vs Q2 based on beam energy es and resolution cut delta 
        // for the elastic peak. For photon energies greater than delta, the reaction is considered non-elastic (real photons radiated).
        elas_radcor(es,delta);   
        
        output(" Incident electron energy = "+ebeam+" GeV",                  writer.outWriter,writer.sumWriter);
        output(" Minimum allowed scattered electron energy = "+epmin0+" GeV",writer.outWriter,writer.sumWriter);
        output(" Target thickness = "+tl_rl+" r.l.",                         writer.outWriter,writer.sumWriter);
        output(" Number of output particles = "+npart,                       writer.outWriter,writer.sumWriter);
    	output(" e-: e_min, th_min, th_max = "+epmin0+" "+th_min+" "+th_max, writer.outWriter,writer.sumWriter);
    	output(" Elastic peak "+(elas_inc==0 ? "included":"excluded"),       writer.outWriter,writer.sumWriter);
    }
    
    public void elast_gen() {

        ntries=0; nevent=0; ntold=0; mcall=0; mcall_max=0;
        boolean fail = false;
        startTime = System.nanoTime();
        
        System.out.println("START");
        do {fail=process(writer);} while (fail || nevent <= nmax); 
        System.out.println("STOP");
        
        System.out.println("Elapsed CPU time: " + cpuTime() + " seconds");
        
    	output(" Elapsed CPU time = "+cpuTime()+" seconds", writer.outWriter);
    	output(" Elapsed CPU time = "+cpuTime()+" seconds", writer.sumWriter);

    	double frac = 0; int intreg=0;
        for (int[] row : dismc) {intreg++;  
        	if (intreg==1) frac = cdiv1;
            if (intreg==2) frac = cdiv2-cdiv1;
            if (intreg==3) frac = cdiv3-cdiv2;
            if (intreg==4) frac = cdiv4-cdiv3;
            if (intreg==5) frac = 1.-cdiv4;            
                output(" distribution of mcall values in the inelastic tail region "+intreg,writer.sumWriter);
                output(" fraction of tail calculation spent here "+frac,writer.sumWriter);
                output(filter(Arrays.toString(row)),writer.sumWriter);
            if (intreg==6){frac = egdiv;
                output(" distribution of mcall values in the elastic peak ",writer.sumWriter);        
                output(" fraction of total integration spent here "+frac,writer.sumWriter);
                output(filter(Arrays.toString(row)),writer.sumWriter);            
            }
        }
        
        close(writer.lundWriter, writer.outWriter, writer.sumWriter); //88,12,14

    }
    
    public boolean process(Writer writer) {
    	
    	double rand, temp, targs, targp, sigr, rtest, stest, ep_sav, mm2, phi_h=0, w, th0d=0;

    	ntries++;
    	
    	// Calculate the energy of the electron at the scattering point after making its way through the target.  
    	// First, randomly choose the interaction point.
    	    
    	rand = Math.random();
    	
        // Change into proper coordinate system

    	vertex_z 	= v_z + tl_cm*(rand-0.5);
    		      
    	// Calculate the incident beam radiation loss
        // Note when f=1 and delta=0 always turn off incoming straggling: opt_strag_0=0

    	es	= ebeam;
    	   
        if (opt_strag_0==1) es = stragl(es, tl_rl*rand);
        
    	sigr1 = sigma_calc(); 
    	
    	if (sigr1==0) return true;
    	
    	// Calculate the distance from interaction to target wall for ep. 
    	// This works for cylindrical target only.  Use GSIM/GEMC for full 
        // target geometry and/or more precise straggling correction.

    	targp	= tr_rl/Math.sin(th0);    	
        temp	= tl_rl*(1-rand)/Math.cos(th0);
    	if (temp < targp) targp = temp;
 
     	sig_tot = sig_tot + sigr1;
    	sigr    = sigr1/sigr_max;
    	
//    	if((int)sigr>10)   System.out.println(ep);
//    	if((int)sigr>1000) System.out.println("ES="+(float)es+" EP="+(float)ep+" TH="+(float)Math.toDegrees(th0)+" SIG="+(float)sigr1);
    	
    	// Calculate the number of times, mcall, to call the routine used to calculate kinematic quantities for the n-tuple.   	
    	
    	rtest = Math.random();
    	mcall = (int) sigr;
    	stest = sigr-mcall;
    	
        if (stest>rtest) mcall++;
        if (mcall>mcall_max) mcall_max=mcall;
        	      
        // If mcall>0 generate mcall n-tuple events.

        if (mcall<1) return true;
        
        if (mcall<10) {
          dismc[intreg-1][mcall-1]++;
        } else {
          dismc[intreg-1][9]++;
        }
        
        ep_sav	= ep;
      
        mm2 = missm();
        
        for (int j=1; j<mcall+1; j++ ) {
        	
        	th0d = Math.toDegrees(th0); 
        	ep = ep_sav;
        	
        	// Calculate the radiation loss for the electron leaving the target
        	// Choose opt_strag_1=0 if this step is performed in GEMC
        	
        	if(opt_strag_1==1) ep = stragl(ep,targp);
        	
            if (ep >= epmin0) {
            
              do {	
            	  phi_r = 2*pi*Math.random();
            	  phi_h = newphi(Math.toDegrees(phi_r));}  
                  while (opt_fiduc == 1 && !accvb(th0d,phi_h,smn,smx));

            
              if (opt_fiduc==2) phi_r = 0.;
              if (opt_strag_1==1) mm2 = missm();

              // Calculate W as if there was no incoming photon radiated.

              w = Math.sqrt(mp*mp+2*mp*(ebeam-ep)-4.*ebeam*ep*Math.pow(Math.sin(th0/2),2));   
            
              double ntp10[] = new double[21];

              ntp10[0]=es;
              ntp10[1]=ep;
              ntp10[2]=th0d;
              ntp10[3]=w;
              ntp10[4]=ppx;
              ntp10[5]=ppy;
              ntp10[6]=ppz;
              ntp10[7]=eprot;
              ntp10[8]=Math.toDegrees(theta_pr);
              ntp10[9]=Math.toDegrees(phi_pr);
              ntp10[10]=mm2;
              ntp10[11]=cstk;
              ntp10[12]=Math.toDegrees(phik);
              ntp10[13]=egam;
              ntp10[14]=egamx;
              ntp10[15]=egamy;
              ntp10[16]=egamz;
              ntp10[17]=vertex_z;
              ntp10[18]=q2;
              ntp10[20]=phi_h;
            
              nevent++;
              
              // Write to LUND file                
          	  try {
				for (String line : ge.getev()) {
					if (line!=null) {
					   writer.lundWriter.write(line);
					   writer.lundWriter.newLine();
					}					
				}	
        	  }  
        	  catch (IOException e)
              {
        		System.out.println("lundWriter write exception!");
              }
            }
        }
        
        //  Talk to the user every now and then.
        
        int ntell = nevent/nprint-ntold;
        if(ntell>0) {
        	int_out();
        	ntold++;
        }
             
        return false;
    }

    public double stragl(double e1, double thik) {
        double gxs, eloss;
    	do {eloss = Math.pow(Math.random(),1./Math.abs(thik));
    	      gxs = 1.-eloss;
        } while (Math.random()>gxs);
        return e1*gxs;
    }
 
    public boolean accvb(double theta, double phi, double thmin, double thmax) {
        double delta_phi = 25*(1.-Math.exp(-0.2*(theta-thmin)));
        return (Math.abs(phi)<delta_phi && theta<thmax);      
    }
    
    public double newphi(double phi) {
        if (phi>330) {
           return phi-360;
        } else if (phi>=0  &&  phi<=30) {
           return phi;
        } else if (phi>0   &&  phi<=90) {
           return phi-60.;
        } else if (phi>90  && phi<=150) {
           return phi-120.;
        } else if (phi>150 && phi<=210) {
           return phi-180.;
        } else if (phi>210 && phi<=270) {
           return phi-240.;
        } else if (phi>270 && phi<=330) {
           return phi-300.;
        }
        return 0;
    }
    
    public double missm() {
    	    	
        // Choose the phi angle of the photon randomly.
        // Choose the hadronic decay angles randomly and calculate the missing mass, the proton momenta and pion momenta.	
    	
    	double csthe,snthe,nu,pp,qx,qz,qvec,sntk,csphk,snphk,cstq,sntq,q2,ps,pps,pprot,pbeam;
    	double csthp,snthp,csphip,snphip,q_dot_pp;
   
        csthe   = Math.cos(th0);
        snthe   = Math.sin(th0);
        nu	    = es-ep;
        ps      = Math.sqrt(es*es-mel*mel);
        pp      = Math.sqrt(ep*ep-mel*mel);
        q2      = 2.*es*ep-2.*ps*pp*csthe-2.*mel*mel;

        // Get lab components of the q vector using es
        // By definition es_x > 0 and qx < 0

        qx      =   -pp*snthe;
        qz      = ps-pp*csthe;
        qvec    = Math.sqrt(qx*qx+qz*qz);

        // Get components of the photon vector in the lab frame

        sntk	= Math.sqrt(1.-cstk*cstk);

        csphk	= Math.cos(phik);
        snphk	= Math.sin(phik);

        cstq	= qz/qvec;
        sntq	= Math.sqrt(1.-cstq*cstq);
        egamx	= egam*(sntk*csphk*cstq-cstk*sntq);
        egamy	= egam*sntk*snphk;
        egamz	= egam*(cstk*cstq+sntk*csphk*sntq);

        // Calculate proton momentum vector in lab frame

        ppz     = qz-egamz;
        ppx     = qx-egamx;
        ppy     =   -egamy;
        pps     = ppx*ppx+ppy*ppy+ppz*ppz;
        pprot   = Math.sqrt(pps);
        eprot   = Math.sqrt(pprot*pprot+mp*mp);
        csthp   = ppz/pprot;
        snthp   = Math.sqrt(1.-csthp*csthp);
        theta_pr= Math.acos(csthp);
        csphip  = ppx/pprot/snthp;
        snphip  = ppy/pprot/snthp;
        phi_pr  = Math.atan2(snphip,csphip);

        // Calculate the square of the missing mass, associated with the proton momentum components and the electron beam energy:

        nu      = ebeam-ep;
        pbeam   = Math.sqrt(ebeam*ebeam-mel*mel);
        q2      = 2.*ebeam*ep-2.*pbeam*pp*csthe-2.*mel*mel;
        	      
        // Get lab components of the q vector using ebeam

        qx      = -pp*snthe;
        qz      = pbeam-pp*csthe;
        q_dot_pp= qx*ppx+qz*ppz;

        return -q2+2*mp*mp+2*mp*(nu-eprot)-2*nu*eprot+2*q_dot_pp;
    
    }

    public void elas_radcor(double es, double delta) {
    	
    	double qs,s2,theta,cst1,eel,epr,gamma4,beta4,e1,e3,e4,eta;
    	double[] deltac = new double[28];

        for (int iq=0; iq<1000; iq++) {
            qs	= iq+1;
            qs	= qs/50.; //valid for 0.02<q2<20 GeV/c
            s2	= qs/4./es/(es-qs/2./mp);
            s2	= Math.sqrt(s2);
            theta	= 2.*Math.asin(s2);

            cst1 = 1.-Math.cos(theta);
            eel	 = es/(1.+es/mp*cst1);
            epr	 = es+mp-eel;
            gamma4	= epr/mp;
            beta4	= Math.sqrt(1.-1/(gamma4*gamma4));
            e1	= es;
            e3	= eel;
            e4	= epr;
            eta	= es/eel;
        
            deltac = elib.del_mo(1,mp,qs,delta,eta,beta4,e1,e3,e4);
            
        	double del_mo  = 0; 
        	for (int idel=0; idel<deltac.length; idel++) del_mo = del_mo + deltac[idel];
        	
        	del_mo = -alpha*del_mo/pi;
        	del_radcor[iq] = Math.exp(del_mo);
        }

    }
    
    public double elas_cor(double es, double theta) {

        int iq;
        double cst1,eel,qs;
        double th0d = Math.toDegrees(theta);
        cst1 = 1.-Math.cos(theta);
        eel	 = es/(1.+es/mp*cst1);
        qs	 = 2.*es*eel*cst1;
        iq	 = (int) qs*50;
        iq	 = iq+1;
    
        if (iq>1000) {
            System.out.println(" Q2 = "+qs+" exceeds programmed range");
            return 0.;
        }
        
        return	elib.elas(es,th0d)*del_radcor[iq];
    }
	
	public double sigma_calc() {
		
		double pp,ps,ucos,cst0,jacob,s,a1,a2,a3,ep0,sig_el;		
		double den,cstk1,cstk2;
		double nu,qvec,mpfac,mcfac;
		double csrnge,csrngb,cran,cran2;
		double sigr1;
		
		double th0d;

        sigr1	= 0;

        ps	= Math.sqrt(es*es-mel*mel);
    
//      Calculate scattering angle

        ucos	= ucos0 + ucrng*Math.random();
        cst0	= 1.-1./ucos;
        jacob	= 1./(ucos*ucos);
    
        th0	    = Math.acos(cst0); th0d = Math.toDegrees(th0);
        s		= Math.pow(Math.sin(th0/2.),2);
    
//      Calculate the energy of the scattered electron in the absence of radiation

        a1	= Math.pow(mel*mel + mp*es,2) + Math.pow(mel*ps*cst0,2);
        a2	= (es+mp)*(mel*mel + mp*es);
        a3	= Math.pow(es+mp,2) - Math.pow(ps*cst0,2);
        ep0	= (a2 + Math.sqrt(a2*a2 - a1*a3))/a3;

//      lcs 1/1/14: If delta=0 then purely elastic with no RC, regardless of egdiv (f)

        if (delta==0.00) { 
           sigr1 = elib.elas(es,th0d);
           ep    = ep0;
           sigr1 = sigr1*jacob/erng/4./pi;
           return sigr1;
        }
   
//      Randomly choose the energy of the photon and scattered electron and calculate Q**2
    
        uek = egdiv>0 ? Math.random():1.0;
    
        if (uek < egdiv) { //event is in elastic peak
           if (elast_inc==1) return sigr1;
           epk	= delta*Math.random();
           ep	= ep0-epk;
           if (ep < epmin0) return sigr1;
           pp	= ep - 0.5*mel*mel/ep;
           jacob = jacob*(delta/erng/egdiv);
        } else {          //event is in elastic tail
           epk	= delta+(erng-delta)*Math.random();
           ep	= ep0-epk;
           if (ep < epmin0) return sigr1;
           pp	= ep-0.5*mel*mel/ep;
           jacob = jacob*((1.-delta/erng)/(1.-egdiv));
        }

        if(uek < egdiv) { //event is in elastic peak      
 	
            intreg	= 6;
            sig_el =  elas_cor(es,th0);
    
//      Make an effective cross section differential in ek,csthk,phik

            sigr1 = sig_el*jacob/delta/4./pi;

            den   = Math.sqrt(ps*ps+pp*pp-2.*ps*pp*cst0);
            cstk1 = (ps-pp*cst0)/den;
            cstk2 = (ps*cst0-pp)/den;
        
            cstk = Math.random()<0.5 ? cstk1:cstk2;    
            phik = 0.;
    
//      The following calculation of the photon energy is approximate

            q2	= 4.*es*ep*s;
            nu	= es-ep;
            qvec = Math.sqrt(q2+nu*nu);
            egam = (2.*mp*nu-q2)/2./(mp+nu-qvec*cstk);
    
            return sigr1;
        }
        
        q2    = 4.*es*ep*s;
        den   = Math.sqrt(ps*ps+pp*pp-2.*ps*pp*cst0);
        cstk1 = (ps-pp*cst0)/den;
        cstk2 = (ps*cst0-pp)/den;
        csrnge = csrng;
    
        if ((1.-cstk1) < csrnge) csrnge=1.-cstk1;
        if (Math.abs(cstk1-cstk2) < 2.*csrnge) csrnge=Math.abs(cstk1-cstk2)/2.;
    
        csrngb	= csrng/40.;
        if (csrngb > csrnge/5.) csrngb=csrnge/5.;
    
       cran	 = Math.random();
       cran2 = Math.random();
    
       cran2 = cran2<0.5 ? -1:1;
       
       if (cran < cdiv1) {
          cstk	= cstk1+csrngb*(2.*Math.random()-1.);
          mcfac	= csrngb/cdiv1;
          phik	= (Math.random()-.5)*pi/9.;
          mpfac	= 1/18.;
	     intreg	= 1;
       } else if (cran < cdiv2) {
          cstk	= cstk1+cran2*(csrngb+(csrnge-csrngb)*Math.random());
          mcfac	= (csrnge-csrngb)/(cdiv2-cdiv1);
          phik	= (Math.random()-.5)*pi/9.;
          mpfac	= 1/18.;
	     intreg = 2;
       } else if (cran < cdiv3) {
         cstk	= cstk2+csrngb*(2.*Math.random()-1.);
         mcfac	= csrngb/(cdiv3-cdiv2);
         phik	= (Math.random()-.5)*pi/9.;
         mpfac	= 1/18.;
	    intreg	= 3;
       } else if (cran < cdiv4) {
           cstk	= cstk2+cran2*(csrngb+(csrnge-csrngb)*Math.random());
          mcfac	= (csrnge-csrngb)/(cdiv4-cdiv3);
           phik	= (Math.random()-.5)*pi/9.;
          mpfac	= 1/18.;
	    intreg	= 4;
       } else {
           boolean test1,test2;   	   
    	   do {
             cstk	= 2.*Math.random()-1.;
             phik	= 2.*pi*(Math.random()-0.5);
             test1 = Math.abs(cstk-cstk1) < csrnge || Math.abs(cstk-cstk2) < csrnge;
             test2 = Math.abs(phik)< pi/18.;
    	   } while (!test1 && !test2);
    	   
//     The following line corrected on Jan. 15, 1999
          mcfac = (1.-2.*csrnge)/(1.-cdiv4);
          mpfac	= 17./18.;
	     intreg = 5;
       }

       sigr1 = elib.radtail(es,ep,th0d,cstk,phik);
       
       egam = elib.egam;
   
       return sigr1 > 0 ? jacob*mpfac*mcfac*sigr1 : 0; //jacob normalizes to the integration region

	}

    public List<String> getInput(String path) { 
    	
    	List<String> out = new ArrayList();
    	String strCurrentLine;
       
    	try (BufferedReader br = new BufferedReader(new FileReader(path))) {
    	while ((strCurrentLine = br.readLine()) != null) out.add(strCurrentLine);
    	} catch (IOException e) {
    		e.printStackTrace();    		
    	}    	
    	return out;
    }    
    
    public String[] parseInput(List<String> list) {
    	
    	ArrayList<String> out = new ArrayList();
    	for (String item : list) {
			String split[] = item.split(" ");
			for (int i=0; i<split.length; i++) out.add(split[i]);
    	}
    	return out.toArray(new String[0]);
    }
    
	abstract class GenEvent {		
		public abstract List<String> getev();		
	}
	
	public class gen_e_p extends GenEvent {	
		List<String> out = new ArrayList<String>();
		public List<String> getev() {
			out.clear();
			out.add(getHeader(2,nevent)); 
			out.add(getString(getParticle(11),1));
			out.add(getString(getParticle(2212),2));
			return out;			
		}		
	}
		
	public class gen_e_p_g extends GenEvent {	
		List<String> out = new ArrayList<String>();
		public List<String> getev() {
			out.clear();
			out.add(getHeader(3,nevent)); 
			out.add(getString(getParticle(11),1));
			out.add(getString(getParticle(2212),2));
            out.add(getString(getParticle(22),3));
			return out;			
		}		
	}
	
	public String getHeader(int np, int n) {
		return np+" "+n+" 0. 0. 0. 0. 0. 0. 0. 0.";		
	}
	
	public String getString(Particle p, int np) {
		return np+" "+p.charge()+" 1 "+p.pid()+" 0 0 "+String.format("%.6f",p.px())+" "
                                                      +String.format("%.6f",p.py())+" "
                                                      +String.format("%.5f",p.pz())+" "
                                                      +"1 "+"1"+" "
                                                      +String.format("%.4f",p.vx())+" "
                                                      +String.format("%.4f",p.vy())+" "
                                                      +String.format("%.4f",p.vz())+" "
                                                      +String.format("%.4f",-3.);			    
	}
	
	public Particle getParticle(int pid) {
		
		double px=0,py=0,pz=0,vx=0,vy=0,vz=0;
				
		switch (pid) {
		case 11: 
			px = Math.sqrt(ep*ep-mel*mel)*Math.sin(th0);
			py = 0;
			pz = Math.sqrt(ep*ep-mel*mel)*Math.cos(th0);
			break;
		case 2212: 
			px = ppx;
			py = ppy;
			pz = ppz;
			break;
		case 22:
			px = (egam<delta) ? 0:egamx;
			py = (egam<delta) ? 0:egamy;
			pz = (egam<delta) ? 0:egamz;
		}
		
		vx = vertex_x; vy = vertex_y; vz= vertex_z;
		
		// Rotate all momenta by a random angle (phi_r) around the z axis	
		double sin_phir = Math.sin(phi_r);
		double cos_phir = Math.cos(phi_r);
		double pxtmp = px, pytmp = py;
		px = pxtmp*cos_phir + pytmp*sin_phir;
		py = pytmp*cos_phir - pxtmp*sin_phir;
		
		return new Particle(pid,px,py,pz,vx,vy,vz);
	}
	
	public void output(String str, BufferedWriter... writer)  {
		for (BufferedWriter w : writer) {
			try {w.write(str); w.newLine();} catch (IOException e){ }
	    }
	}
	
	public void close(BufferedWriter... writer) {
		for (BufferedWriter w : writer) {
			try {w.close();} catch (IOException e){ }
		}
	}
	
	public String filter(String in) {
		return in.replaceAll("\\[", "").replaceAll("\\]", "").replaceAll("\\,", " ");
	}
	
	public double cpuTime() {
	  return (System.nanoTime() - startTime)/1e9;
    }
   
    public class Writer {
    	
    	private File lundFile = new File("elastgen.lund"); 
    	private File  outFile = new File("elastgen.out");
    	private File  sumFile = new File("elastgen.sum");
    	
    	public BufferedWriter lundWriter,outWriter,sumWriter; //88,12,14 in elast_gen.F
    	
    	public Writer() {
    		try {    		
        		lundWriter = new BufferedWriter(new FileWriter(lundFile.getAbsoluteFile()));
        		outWriter  = new BufferedWriter(new FileWriter(outFile.getAbsoluteFile()));
        		sumWriter  = new BufferedWriter(new FileWriter(sumFile.getAbsoluteFile()));
    		}
    		catch (IOException e)
    		{
                System.out.println("Writer Exception ");       
            }
    		
    	}
    }
    	    
    public static void main(String[] args) { 
    	
    	Elasgen gen = new Elasgen();
    	gen.getFile("/Users/colesmith/clas12/ElasticGen/inp/elas_clas12.inp");
    	gen.run();

    }
}
