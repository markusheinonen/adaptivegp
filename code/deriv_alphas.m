function [dl,ds,do] = deriv_alphas(gp)
% derivative of the alphas over MLL


	Kl = gausskernel(gp.xtr,gp.xtr,gp.betaell,   gp.alphaell, gp.tol);
	Ks = gausskernel(gp.xtr,gp.xtr,gp.betasigma, gp.alphasigma, gp.tol);
	Ko = gausskernel(gp.xtr,gp.xtr,gp.betaomega, gp.alphaomega, gp.tol);
		
	al = Kl\(gp.l_ell   - gp.l_muell);
	as = Ks\(gp.l_sigma - gp.l_musigma);
	ao = Ko\(gp.l_omega - gp.l_muomega);

	dKl = 2 * gp.alphaell^(-1)   * Kl;
	dKs = 2 * gp.alphasigma^(-1) * Ks;
	dKo = 2 * gp.alphaomega^(-1) * Ko;
	
	dl = 0.5*sum(diag((al*al' - inv(Kl))*dKl));
	ds = 0.5*sum(diag((as*as' - inv(Ks))*dKs));
	do = 0.5*sum(diag((ao*ao' - inv(Ko))*dKo));
end


