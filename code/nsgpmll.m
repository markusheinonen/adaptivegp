function [val,valy,vall,vals,vale,vall1,vall2] = nsgpmll(pars)

	[ell,sigma,omega] = latentchols(pars);
	
	n = length(pars.xtr);
	if sum(ismember('ab', pars.nsfuncs))
		pars.Kl = gausskernel(pars.xtr, pars.xtr, pars.betaell,   pars.alphaell,   pars.tol);
		pars.Ks = gausskernel(pars.xtr, pars.xtr, pars.betasigma, pars.alphasigma, pars.tol);
		pars.Ko = gausskernel(pars.xtr, pars.xtr, pars.betaomega, pars.alphaomega, pars.tol);
	end
	Ky = nsgausskernel(pars.xtr, pars.xtr, ell, ell, sigma, sigma, omega);
	zs = zeros(n,1);
	
	% check if non-sdp or low condition number
	[~,p] = chol(Ky);
	rc = rcond(Ky);
	if p > 0 || rc < 1e-15
		val = -inf;
		return;
	end
	
	% assuming exp-transformation here
	valy = logmvnpdf(pars.ytr, zs, Ky);
	
	if length(pars.ell) == 1
		vall = logmvnpdf(pars.ell*ones(n,1), pars.muell, pars.Kl);
	else
		[vall,vall1,vall2] = logmvnpdf(ell, pars.muell, pars.Kl);
	end
	
	if length(pars.sigma) == 1
		vals = logmvnpdf(pars.sigma*ones(n,1), pars.musigma, pars.Ks);
	else
		vals = logmvnpdf(sigma, pars.musigma, pars.Ks);
	end
	
	if length(pars.omega) == 1
		vale = logmvnpdf(pars.omega*ones(n,1), pars.muomega, pars.Ko);
	else
		vale = logmvnpdf(omega, pars.muomega, pars.Ko);
	end
	
	val = valy + vall + vals + vale;
end



