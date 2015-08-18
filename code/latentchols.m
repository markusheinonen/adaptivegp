function [ell,sigma,omega] = latentchols(pars)
	
	ell = pars.ell;
	sigma = pars.sigma;
	omega = pars.omega;

	if isfield(pars, 'Ll')
		ell = pars.Ll * pars.ell;
	end
	if isfield(pars, 'Ls')
		sigma = pars.Ls * pars.sigma;
	end
	if isfield(pars, 'Lo')
		omega = pars.Lo * pars.omega;
	end
end