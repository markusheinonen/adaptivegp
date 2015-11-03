function [fmean,fstd,lt,st,ot] = nsgpposterior(gp, xt)
% Model posteriors over target points 'xt'
% returns:
%  - fmean      : f posterior mean
%  - fstd       : f posterior std
%  - lt         :   ell posterior mean [lengthscale]
%  - st         : sigma posterior mean [signal variance]
%  - ot         : omega posterior mean [noise variance]
%  - fderivmean : f derivative posterior mean
%  - fderivstd  : f derivative posterior std

	if ~exist('xt','var')
		xt = gp.xtr;
	end
	
	ot = gpposterior(gp.xtr, gp.l_omega, xt, gp.l_muomega, gp.betaomega, gp.alphaomega, gp.tol);
	lt = gpposterior(gp.xtr, gp.l_ell,   xt, gp.l_muell,   gp.betaell,   gp.alphaell,   gp.tol);
	st = gpposterior(gp.xtr, gp.l_sigma, xt, gp.l_musigma, gp.betasigma, gp.alphasigma, gp.tol);

	% extrapolate unknown function
	Ktt = nsgausskernel(gp.xtr, gp.xtr, gp.l_ell, gp.l_ell,   gp.l_sigma, gp.l_sigma, gp.l_omega);
	Kts = nsgausskernel(gp.xtr, xt,     gp.l_ell, lt,         gp.l_sigma, st,         log(0));
	Kss = nsgausskernel(xt,     xt,     lt,       lt,         st,         st,         log(0));
	Kst = Kts';
	
	A = Kst / Ktt;
	
	fmean = A*(gp.ytr - gp.muf) + gp.muf;
	fcov = Kss - A*Kts;
	fstd = sqrt(diag(fcov));
	
	% compute also the derivative GP
%	if nargout > 5
%		Kdst = nsgausskernelderiv(xt,  gp.xtr, lt, gp.l_ell, st, gp.l_sigma);
%		Kdss = nsgausskernelderiv2(xt, xt,       lt, lt,  st, st);
%		fderivmean = Kdst / Ktt * gp.ytr;
%		fderivstd = sqrt(diag(Kdss - Kdst / Ktt * Kdst'));
%		
%		[~,~,fmean,fstd,lt,st,ot,fderivmean,fderivstd] = denormalise(gp, fmean, fstd, exp(lt), exp(st), exp(ot), fderivmean, fderivstd);
%	end

	lt = exp(lt);
	st = exp(st);
	ot = exp(ot);

	[~,~,fmean,fstd,lt,st,ot] = denormalise(gp, fmean, fstd, lt, st, ot);
end

