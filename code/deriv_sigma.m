function [ds] = deriv_sigma(pars, scalar)
% derivative of the sigma latent function wrt MLL
	
	if ~exist('scalar','var')
		scalar = 0;
	end

	[ell,sigma,omega] = latentchols(pars);
	
	n = length(pars.xtr);
	
	if sum(ismember('ab', pars.nsfuncs))
		pars.Ks = gausskernel(pars.xtr,pars.xtr,pars.betasigma, pars.alphasigma, pars.tol);
	end
	
	Ky = nsgausskernel(pars.xtr, pars.xtr, ell, ell, sigma, sigma, omega);
	Kf = nsgausskernel(pars.xtr, pars.xtr, ell, ell, sigma, sigma, log(0));
	
	a = Ky\pars.ytr;
	A = a*a' - inv(Ky);

	ds = 2*diag( A * Kf ) - pars.Ks\(sigma-pars.musigma);
	
	if scalar
		ds = ones(n,1) * sum(ds);
	end
	
	if isfield(pars, 'Ls')
		ds = pars.Ls'*ds;
	end
end


