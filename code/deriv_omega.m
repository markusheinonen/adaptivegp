function [do] = deriv_omega(pars, scalar)
% derivative of the noise latent function over the MLL
%

	if ~exist('scalar','var')
		scalar = 0;
	end
	
	n = length(pars.xtr);
	
	[ell,sigma,omega] = latentchols(pars);	
	
	if sum(ismember('ab', pars.nsfuncs))
		pars.Ko = gausskernel(pars.xtr,pars.xtr,pars.betaomega, pars.alphaomega, pars.tol);
	end
		
	Ky = nsgausskernel(pars.xtr,pars.xtr,ell, ell, sigma, sigma, omega);
	a = Ky\pars.ytr;
	A = a*a' - inv(Ky);

	dK = diag( 2.*exp(2*omega) );
	do = 0.5*diag(A*dK) - pars.Ko\(omega - pars.muomega);

	if scalar
		do = ones(n,1) * sum(do);
	end
	
	if isfield(pars, 'Lo')
		do = pars.Lo'*do;
	end
end


