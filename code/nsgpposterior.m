function [pmean,pstd,lt,st,ot,pmeanderiv,pstdderiv] = nsgpposterior(pars, xt)
% non-stationary (scalar) gaussian kernel
	
	% extrapolate latent functions
	[ell,sigma,omega] = latentchols(pars);
	ot = gpposterior(pars.xtr, omega, xt, pars.muomega, pars.betaomega, pars.alphaomega, pars.tol);
	lt = gpposterior(pars.xtr, ell,   xt, pars.muell,   pars.betaell,   pars.alphaell,   pars.tol);
	st = gpposterior(pars.xtr, sigma, xt, pars.musigma, pars.betasigma, pars.alphasigma, pars.tol);
	
	% extrapolate unknown function
	Ktt = nsgausskernel(pars.xtr, pars.xtr, ell, ell,   sigma, sigma, omega);
	Kts = nsgausskernel(pars.xtr, xt,       ell, lt,    sigma, st,    log(0));
	Kss = nsgausskernel(xt,       xt,       lt,  lt,    st,    st,    log(0));
	Kst = Kts';
	
	A = Kst / Ktt;
	
	pmean = A*(pars.ytr-mean(pars.ytr));
	pcov = Kss - A*Kts;
	pstd = sqrt(diag(pcov));
	
	% compute also the derivative GP
	if nargout > 5
		Kdst = nsgausskernelderiv(xt,  pars.xtr, lt, ell, st, sigma);
		Kdss = nsgausskernelderiv2(xt, xt,       lt, lt,  st, st);
		pmeanderiv = Kdst / Ktt * pars.ytr;
		pstdderiv = sqrt(diag(Kdss - Kdst / Ktt * Kdst'));
	end
	
	ot = exp(ot);
	lt = exp(lt);
	st = exp(st);	
end

