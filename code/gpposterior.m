function [pmean,pcov] = gpposterior(x, y, xt, my, ell, sigma, omega)
% non-stationary (scalar) gaussian kernel
		
	Ktt = gausskernel(x,x,ell,sigma, omega);
	Kts = gausskernel(x,xt,ell,sigma, 0);
	Kss = gausskernel(xt,xt,ell,sigma, 0);
	Kst = Kts';
	
	A = Kst / Ktt;
	
	pmean = A*(y-my) + my;
	pcov = Kss - A*Kts;
end



