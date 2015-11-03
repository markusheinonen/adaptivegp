function [val,valy,vall,vals,vale] = nsgpmll(pars)
	
	n = length(pars.xtr);
	if sum(ismember('ab', pars.nsfuncs))
		pars.Kl = gausskernel(pars.xtr, pars.xtr, pars.betaell,   pars.alphaell,   pars.tol);
		pars.Ks = gausskernel(pars.xtr, pars.xtr, pars.betasigma, pars.alphasigma, pars.tol);
		pars.Ko = gausskernel(pars.xtr, pars.xtr, pars.betaomega, pars.alphaomega, pars.tol);
	end
	Ky = nsgausskernel(pars.xtr, pars.xtr, pars.l_ell, pars.l_ell, pars.l_sigma, pars.l_sigma, pars.l_omega);
	zs = zeros(pars.n,pars.p);
	
	% check if non-sdp or low condition number
	[~,p] = chol(Ky);
	rc = rcond(Ky);
	if p > 0 || rc < 1e-15
		val = -inf;
		return;
	end
	
	% assuming exp-transformation here
	valy = diag(logmvnpdf(pars.ytr, zs, Ky));
	
	if length(pars.ell) == 1
		vall = logmvnpdf(pars.ell*ones(n,1), muell, pars.Kl);
	else
		vall = logmvnpdf(pars.l_ell, pars.l_muell, pars.Kl);
	end
	
	if length(pars.sigma) == 1
		vals = logmvnpdf(pars.sigma*ones(n,1), musigma, pars.Ks);
	else
		vals = logmvnpdf(pars.l_sigma, pars.l_musigma, pars.Ks);
	end
	
	if length(pars.omega) == 1
		vale = logmvnpdf(pars.omega*ones(n,1), muomega, pars.Ko);
	else
		vale = logmvnpdf(pars.l_omega, pars.l_muomega, pars.Ko);
	end
	
	val = valy + vall + vals + vale;
end



