%% 
% Supplementary material 
%
% From Weisse et al., "A mechanistic link between cellular trade-offs, 
% gene expression and growth", PNAS, 2015
%
% This file implements the right hand side of the ODE.

%%
function dydt= cellmodel_odes(t, y, rates, parameters)

	b= rates(1);
	dm= rates(2);
	kb= rates(3);
	ku= rates(4);
	f= rates(5);

	thetar= parameters(1);
	k_cm= parameters(2);
	s0= parameters(3);
	gmax= parameters(4);
	cl= parameters(5);
	thetax= parameters(6);
	Kt= parameters(7);
	M= parameters(8);
	we= parameters(9);
	Km= parameters(10);
	vm= parameters(11);
	nx= parameters(12);
	Kq= parameters(13);
	Kp= parameters(14);
	vt= parameters(15);
	wr= parameters(16);
	wq= parameters(17);
	wp= parameters(18);
	nq= parameters(19);
	nr= parameters(20);
	ns= parameters(21);

	rmr= y(1);
	em= y(2);
	rmp= y(3);
	rmq= y(4);
	rmt= y(5);
	et= y(6);
	rmm= y(7);
	zmm= y(8);
	zmr= y(9);
	zmp= y(10);
	zmq= y(11);
	zmt= y(12);
	mt= y(13);
	mm= y(14);
	q= y(15);
	p= y(16);
	si= y(17);
	mq= y(18);
	mp= y(19);
	mr= y(20);
	r= y(21);
	a= y(22);

	Kgamma= gmax/Kp;
	gamma= gmax*a/(Kgamma + a);
	ttrate= (rmq + rmr + rmp + rmt + rmm)*gamma;
	lam= ttrate/M;
	fr= nr*(r + rmr + rmp + rmt + rmm + rmq + zmr + zmp + zmt + zmm + zmq) / ( nr*(r + rmr + rmp + rmt + rmm + rmq + zmr + zmp + zmt + zmm + zmq) + nx * (p + q + et + em));
	nucat= em*vm*si/(Km + si);

	dydt(size(y,1),1)= 0;
	dydt(1)= +kb*r*mr+b*zmr-ku*rmr-gamma/nr*rmr-f*rmr-lam*rmr;
	dydt(2)= +gamma/nx*rmm-lam*em;
	dydt(3)= +kb*r*mp+b*zmp-ku*rmp-gamma/nx*rmp-f*rmp-lam*rmp;
	dydt(4)= +kb*r*mq+b*zmq-ku*rmq-gamma/nx*rmq-f*rmq-lam*rmq;
	dydt(5)= +kb*r*mt+b*zmt-ku*rmt-gamma/nx*rmt-f*rmt-lam*rmt;
	dydt(6)= +gamma/nx*rmt-lam*et;
	dydt(7)= +kb*r*mm+b*zmm-ku*rmm-gamma/nx*rmm-f*rmm-lam*rmm;
	dydt(8)= +f*rmm-b*zmm-lam*zmm;
	dydt(9)= +f*rmr-b*zmr-lam*zmr;
	dydt(10)= +f*rmp-b*zmp-lam*zmp;
	dydt(11)= +f*rmq-b*zmq-lam*zmq;
	dydt(12)= +f*rmt-b*zmt-lam*zmt;
	dydt(13)= +(we*a/(thetax + a))+ku*rmt+gamma/nx*rmt-kb*r*mt-dm*mt-lam*mt;
	dydt(14)= +(we*a/(thetax + a))+ku*rmm+gamma/nx*rmm-kb*r*mm-dm*mm-lam*mm;
	dydt(15)= +gamma/nx*rmq-lam*q;
	dydt(16)= +gamma/nx*rmp-lam*p;
	dydt(17)= +(et*vt*s0/(Kt + s0))-nucat-lam*si;
	dydt(18)= +(wq*a/(thetax + a)/(1 + (q/Kq)^nq))+ku*rmq+gamma/nx*rmq-kb*r*mq-dm*mq-lam*mq;
	dydt(19)= +(wp*a/(thetax + a))+ku*rmp+gamma/nx*rmp-kb*r*mp-dm*mp-lam*mp;
	dydt(20)= +(wr*a/(thetar + a))+ku*rmr+gamma/nr*rmr-kb*r*mr-dm*mr-lam*mr;
	dydt(21)= +ku*rmr+ku*rmt+ku*rmm+ku*rmp+ku*rmq+gamma/nr*rmr+gamma/nr*rmr+gamma/nx*rmt+gamma/nx*rmm+gamma/nx*rmp+gamma/nx*rmq-kb*r*mr-kb*r*mt-kb*r*mm-kb*r*mp-kb*r*mq-lam*r;
	dydt(22)= +ns*nucat-ttrate-lam*a;
