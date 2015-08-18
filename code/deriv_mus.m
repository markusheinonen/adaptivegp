function [dl, ds, do] = deriv_mus(pars)
% derivative of the mean parameters against MLL

	% update if necessary
	if sum(ismember('ab', pars.nsfuncs))
		pars.Kl = gausskernel(pars.xtr,pars.xtr,pars.betaell,   pars.alphaell,   pars.tol);
		pars.Ks = gausskernel(pars.xtr,pars.xtr,pars.betasigma, pars.alphasigma, pars.tol);
		pars.Ko = gausskernel(pars.xtr,pars.xtr,pars.betaomega, pars.alphaomega, pars.tol);
	end
	
	[ell,sigma,omega] = latentchols(pars);
	
	dl = sum(pars.Kl\(ell - pars.muell));
	ds = sum(pars.Ks\(sigma - pars.musigma));
	do = sum(pars.Ko\(omega - pars.muomega));	
end


