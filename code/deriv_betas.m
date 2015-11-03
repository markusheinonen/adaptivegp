function [dl,ds,do] = deriv_betas(gp)
% derivative of the betas over MLL

	Kl = gausskernel(gp.xtr,gp.xtr,gp.betaell,   gp.alphaell, gp.tol);
	Ks = gausskernel(gp.xtr,gp.xtr,gp.betasigma, gp.alphasigma, gp.tol);
	Ko = gausskernel(gp.xtr,gp.xtr,gp.betaomega, gp.alphaomega, gp.tol);
	
	al = Kl\(gp.l_ell   - gp.l_muell);
	as = Ks\(gp.l_sigma - gp.l_musigma);
	ao = Ko\(gp.l_omega - gp.l_muomega);

	dKl = gp.betaell^(-3)   * gp.D .* Kl;
	dKs = gp.betasigma^(-3) * gp.D .* Ks;
	dKo = gp.betaomega^(-3) * gp.D .* Ko;
	
	dl = 0.5*sum(diag((al*al' - inv(Kl))*dKl));% + 0.5*pars.betaell^(-1) - 2; % gamma prior
	ds = 0.5*sum(diag((as*as' - inv(Ks))*dKs));% + 0.5*pars.betasigma^(-1) - 2;
	do = 0.5*sum(diag((ao*ao' - inv(Ko))*dKo));% + 0.5*pars.betaomega^(-1) - 2;	
end


