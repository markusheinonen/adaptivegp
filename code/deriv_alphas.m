function [dl,ds,do] = deriv_alphas(pars)
% derivative of the alphas over MLL


	Kl = gausskernel(pars.xtr,pars.xtr,pars.betaell,   pars.alphaell, pars.tol);
	Ks = gausskernel(pars.xtr,pars.xtr,pars.betasigma, pars.alphasigma, pars.tol);
	Ko = gausskernel(pars.xtr,pars.xtr,pars.betaomega, pars.alphaomega, pars.tol);
	
	[ell,sigma,omega] = latentchols(pars);

	al = Kl\(ell   - pars.muell);
	as = Ks\(sigma - pars.musigma);
	ao = Ko\(omega - pars.muomega);

	dKl = 2 * pars.alphaell^(-1)   * Kl;
	dKs = 2 * pars.alphasigma^(-1) * Ks;
	dKo = 2 * pars.alphaomega^(-1) * Ko;
	
	dl = 0.5*sum(diag((al*al' - inv(Kl))*dKl));
	ds = 0.5*sum(diag((as*as' - inv(Ks))*dKs));
	do = 0.5*sum(diag((ao*ao' - inv(Ko))*dKo));
end


