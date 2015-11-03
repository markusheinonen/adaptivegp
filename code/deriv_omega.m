function dwl_o = deriv_omega(pars, scalar)
% derivative of the noise latent function over the MLL
%

	if ~exist('scalar','var')
		scalar = 0;
	end
	
	n = length(pars.xtr);
		
	if sum(ismember('ab', pars.nsfuncs))
		pars.Ko = gausskernel(pars.xtr,pars.xtr,pars.betaomega, pars.alphaomega, pars.tol);
	end
		
	Ky = nsgausskernel(pars.xtr,pars.xtr,pars.l_ell, pars.l_ell, pars.l_sigma, pars.l_sigma, pars.l_omega);
	a = Ky\pars.ytr;
	A = a*a' - inv(Ky);

	dK = diag( 2.*exp(2*pars.l_omega) );
	dl_o = 0.5*diag(A*dK) - pars.Ko\(pars.l_omega - pars.l_muomega);

	if scalar
		dl_o = ones(n,1) * sum(dl_o);
	end
	
%	if isfield(pars, 'Lo')
	if ismember('o', pars.nsfuncs)
		dwl_o = pars.Lo'*dl_o;
	else
		dwl_o = pars.Lo\dl_o;
	end
end


