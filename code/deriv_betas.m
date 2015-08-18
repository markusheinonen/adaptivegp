function [dl,ds,do] = deriv_betas(pars)
% derivative of the betas over MLL

	Kl = gausskernel(pars.xtr,pars.xtr,pars.betaell,   pars.alphaell, pars.tol);
	Ks = gausskernel(pars.xtr,pars.xtr,pars.betasigma, pars.alphasigma, pars.tol);
	Ko = gausskernel(pars.xtr,pars.xtr,pars.betaomega, pars.alphaomega, pars.tol);
	
	[ell,sigma,omega] = latentchols(pars);

	al = Kl\(ell   - pars.muell);
	as = Ks\(sigma - pars.musigma);
	ao = Ko\(omega - pars.muomega);

	dKl = pars.betaell^(-3)   * pars.D .* Kl;
	dKs = pars.betasigma^(-3) * pars.D .* Ks;
	dKo = pars.betaomega^(-3) * pars.D .* Ko;
	
	dl = 0.5*sum(diag((al*al' - inv(Kl))*dKl));% + 0.5*pars.betaell^(-1) - 2; % gamma prior
	ds = 0.5*sum(diag((as*as' - inv(Ks))*dKs));% + 0.5*pars.betasigma^(-1) - 2;
	do = 0.5*sum(diag((ao*ao' - inv(Ko))*dKo));% + 0.5*pars.betaomega^(-1) - 2;	
end


