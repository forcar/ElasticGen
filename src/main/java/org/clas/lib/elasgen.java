package org.clas.lib;

public class elasgen {
	
	MoTsai elib = new MoTsai();

	double mel,es,ep,ucos0,ucrng,erng,delta;
	double egdiv,epin0,cdiv1,cdiv2,cdiv3,cdiv4,crng,elast_inc,intreg;
	double uek,epk;
	double epmin0;
    double  pi = 3.141592654;
    double  mp = 0.93827;
    double rme = 0.00051099895069;
    double alpha = 1./137.01;
    double hbarc = 389.37966;
    double[] del_radcor = new double[1000];
    
    public elasgen() {
		
    } 
    
    public void elas_radcor(double es, double delta) {
    	
    	double qs,s2,theta,snth,cst1,eel,epr,gamma4,beta4,e1,e3,e4,eta;
    	double[] deltac = new double[28];

        for (int iq=0; iq<1000; iq++) {
            qs	= iq+1;
            qs	= qs/50.;
            s2	= qs/4./es/(es-qs/2./mp);
            s2	= Math.sqrt(s2);
            theta	= 2.*Math.asin(s2);

            snth = Math.sin(theta);
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

        cst1 = 1.-Math.cos(theta);
        eel	 = es/(1.+es/mp*cst1);
        qs	 = 2.*es*eel*cst1;
        iq	 = (int) qs*50;
        iq	 = iq+1;
    
        if (iq>1000) {
            System.out.println(" Q2 = "+qs+" exceeds programmed range");
            return 0.;
        }

        return	elib.elas(es,theta)*del_radcor[iq];
    }
	
	public double sigma_calc(double es, double delta, double csrng) {
		
		double pp,ps,ucos,cst0,jacob,snt0,s,a1,a2,a3,ep0,sig_el;		
		double den,cstk1,cstk2,Tk;
		double nu,qvec,mpfac,mcfac;
		double csrnge,csrngb,cran,cran2;
		double sigr1;
		
		double ep,th0,th0d,ppx,ppy,ppz,eprot,egam,egamx,egamy,egamz,cstk,phik,q2;

        sigr1	= 0;

        ps	= Math.sqrt(es*es-rme*rme);
    
//      Calculate scattering angle

        ucos	= ucos0 + ucrng*Math.random();
        cst0	= 1.-1./ucos;
        jacob	= 1./(ucos*ucos);
    
        th0	    = Math.acos(cst0); th0d = Math.toRadians(th0);
        snt0	= Math.sin(th0);
        s		= Math.pow(Math.sin(th0/2.),2);
    
//      Calculate the energy of the scattered electron in the absence of radiation

        a1	= Math.pow(mel*mel + mp*es,2) + Math.pow(mel*ps*cst0,2);
        a2	= (es+mp)*(mel*mel + mp*es);
        a3	= Math.pow(es+mp,2) - Math.pow(ps*cst0,2);
        ep0	= (a2 + Math.sqrt(a2*a2 - a1*a3))/a3;

// lcs 1/1/14: If delta=0 then purely elastic with no RC, regardless of egdiv (f)

        if (delta==0.00) { 
           sigr1 = elib.elas(es,th0d);
           ep    = ep0;
           sigr1 = sigr1*jacob/erng/4./pi;
           return sigr1;
        }
   
//      Randomly choose the energy and momentum of the scattered electron and calculate Q**2
    
        uek = egdiv>0 ? Math.random():1.0;
    
        boolean inElasticPeak = uek < egdiv ? true:false;
        
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
           if (ep<epmin0) return sigr1;
           pp	= ep-0.5*mel*mel/ep;
           jacob = jacob*((1.-delta/erng)/(1.-egdiv));
//           go to 41
        }

        if(uek < egdiv) { //event is in elastic peak
        	
            if (elast_inc==1) return sigr1;
            
            epk	= delta*Math.random();
            ep	= ep0-epk;
            
            if (ep < epmin0) return sigr1;
            
            pp	= ep - 0.5*mel*mel/ep;
            jacob = jacob*(delta/erng/egdiv);

            intreg	= 6;
            sig_el =  elas_cor(es,th0d);
    
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
    
            mpfac	= 1.;
            mcfac	= 1.;
        
            return mpfac*mcfac*sigr1;
        }
        
//        go to 48
  
     
        // event is in elastic tail

        epk	  = delta+(erng-delta)*Math.random();
        ep	  = ep0-epk;
        if (ep<epmin0) return sigr1;
        
        pp	  = ep-0.5*mel*mel/ep;
        jacob = jacob*((1.-delta/erng)/(1.-egdiv));
        
      //41 
        
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

       Tk = Math.acos(cstk);

       elib.radtail(es,ep,th0d,cstk,phik);

//48     
       return sigr1 > 0 ? jacob*mpfac*mcfac*sigr1 : 0; //jacob normalizes to the integration region

	}
}
