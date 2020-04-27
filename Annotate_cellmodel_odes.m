%% 
% This file implements the right hand side of the ODE.



function dydt = Annotate_cellmodel_odes(t, y, rates, parameters)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% UNPACK COMPONENTS
    
    % Unpack paramters (see _Driver)
	s0 = parameters(1);
	ns = parameters(2);
	nr = parameters(3);
	nx = parameters(4);
	gmax = parameters(5);
	Kp = parameters(6);
	vet = parameters(7);
	Ket = parameters(8);
	vem = parameters(9);
	Kem = parameters(10);
	wr = parameters(11);
	we = parameters(12);
	wq = parameters(13);
	wp = parameters(14);
	thetar = parameters(15);
	thetax = parameters(16);
	Kq = parameters(17);
	nq = parameters(18);
	M = parameters(19);
    
    % Unpack rates
	dm = rates(1);
	kb = rates(2);
	ku = rates(3);
    
    % Populate variables
	r = y(1);
	et = y(2);
	em = y(3);
	q = y(4);
	p = y(5);     
	mr = y(6);
	mt = y(7);
	mm = y(8);
	mq = y(9);
	mp = y(10);	        
    cr = y(11);
	ct = y(12);
	cm = y(13);
	cq = y(14);
	cp = y(15);    
    si = y(16);
	a = y(17);
    
    % Useful to make variables for the sum of mRNA/complex
    m_all = mr + mt + mm + mq + mp;
    c_all = cr + ct + cm + cq + cp;
       
    
    
    %% CALCULATE RATES
    
    % Translation elongation thresh., ~=7 (Kp defined in _Driver)
	Kgamma = gmax/Kp;    
    % Translation elongation rate
	gamma = gmax*a/(Kgamma + a);    
    % Total translation rate = nx*(gamma/nx*rmx) = gamma*rmx
	ttrate = gamma*c_all;    
    % Dilution rate = growth rate
	lam = ttrate/M;    
    % Rate of import of nutrients
    nuimp = et*vet*s0/(Ket + s0);
    % Rate of metabolism into energy
    nucat = em*vem*si/(Kem + si);
    
    
    
    %% CREATE ODEs
        
    dydt(size(y,1),1) = 0; % Setup ODE structure
    
    % r (free ribo) = transl - dilu + transl_all - bind_all + unbind_all
    %%% Split transl_all due to different nr/nx values
    dydt(1) = gamma*cr/nr - lam*r + (gamma*cr/nr + gamma*(c_all-cr)/nx) - kb*r*m_all + ku*c_all;
    % et (free transporters) = transl - dilution
    dydt(2) = gamma*ct/nx - lam*et;	
    % em (metabolic base) = ""
	dydt(3) = gamma*cm/nx - lam*em;    
    % q (free housekeeping) = ""
    dydt(4) = gamma*cq/nx - lam*q;	
    % p (free placeholder) = ""
    dydt(5) = gamma*cp/nx - lam*p;
    
    % mr (mRNA, ribo) = transcr - bind + unbind + transl - degra - dilution
    dydt(6) = wr*a/(thetar + a) - kb*r*mr + ku*cr + gamma*cr/nr - dm*mr - lam*mr;	
    % mt (mRNA, transporters) = ""
    dydt(7) = we*a/(thetax + a) - kb*r*mt + ku*ct + gamma*ct/nx - dm*mt - lam*mt;	
    % mm (mRNA, metabolic) = ""
    dydt(8) = we*a/(thetax + a) - kb*r*mm + ku*cm + gamma*cm/nx - dm*mm - lam*mm;	
    % mq (mRNA, housekeeping) = transcr*negautoreg + ""
    dydt(9) = (wq*a/(thetax + a))*(1/(1+(q/Kq)^nq)) - kb*r*mq + ku*cq + gamma*cq/nx - dm*mq - lam*mq;	
    % mp (mRNA, PLACEHOLDER) = ""
    dydt(10) = wp*a/(thetax + a) - kb*r*mp + ku*cp + gamma*cp/nx - dm*mp - lam*mp;	
    
    % cr (complex, ribo) = bind - unbind - transl - dilution    
	dydt(11) = kb*r*mr - ku*cr - gamma*cr/nr - lam*cr;    
    % ct (complex, transporters) = ""
    dydt(12) = kb*r*mt - ku*ct - gamma*ct/nx - lam*ct;	
    % cm (complex, metabolic) = ""
    dydt(13) = kb*r*mm - ku*cm - gamma*cm/nx - lam*cm;	
    % cq (complex, housekeeping) = ""  
    dydt(14) = kb*r*mq - ku*cq - gamma*cq/nx - lam*cq;	
    % cp (complex, PLACEHOLDER) = ""  
	dydt(15) = kb*r*mp - ku*cp - gamma*cp/nx - lam*cp;	
    
    % si (intracellular nutrient) = import - metabolism - dilution
    dydt(16) = nuimp - nucat - lam*si;	
    % a (energy) = metabolism - translation - dilution
    dydt(17) = ns*nucat - ttrate - lam*a;