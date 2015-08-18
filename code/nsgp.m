function [gpmodel,samples,mll,mse,nmse,nlpd] = nsgp(x,y,nsfuncs,optim, gpmodel)
% main function for nonstationary GP modeling
% 
% INPUTS
% nsfuncs : string of nonstationary components, 'lsemab', options:
%           l : lengthscale
%           s : signal variance
%           o : noise std
%           m : prior means
%           a : prior variances
%           b : prior lenghtscales
% optim   : either 'hmc' or 'grad'
%
% OUTPUTS
% pars    : parameter values learned [grad] or initialised [hmc]
% samples : hmc samples of the posterior? [hmc] or empty [grad]
%


	if ~exist('nsfuncs','var')
		nsfuncs = 'lso';
	end
	if ~exist('optim','var')
		optim = 'grad';
	end
	
	% filter NaN's
	nonnuls = ~isnan(y) & prod(~isnan(x),2);
	x = x(nonnuls,:);
	y = y(nonnuls);

	if ~exist('pars','var')
		gpmodel = defaultpars(x);
	end
	
	% normalize
	n = length(y);
	gpmodel.xbias = min(x);
	gpmodel.xscale = max(x)-min(x);
	gpmodel.yscale = max(abs(y-mean(y)));
	gpmodel.ybias = mean(y);
	gpmodel.xtr = (x - repmat(gpmodel.xbias,n,1)) ./ repmat(gpmodel.xscale,n,1);
	gpmodel.ytr = (y - gpmodel.ybias) / gpmodel.yscale;

%	pars.xtr = x;
%	pars.ytr = y;
	gpmodel.xts = gpmodel.xtr;
	gpmodel.yts = gpmodel.ytr;
		
	gpmodel.nsfuncs = nsfuncs;
	gpmodel.optim = optim;
	gpmodel.n = length(gpmodel.xtr);
	
	if length(gpmodel.ell) == 1
		gpmodel.ell = gpmodel.ell * ones(gpmodel.n,1);
	end
	if length(gpmodel.sigma) == 1
		gpmodel.sigma = gpmodel.sigma * ones(gpmodel.n,1);
	end
	if length(gpmodel.omega) == 1
		gpmodel.omega = gpmodel.omega * ones(gpmodel.n,1);
	end
	
	% precompute these when we are not optimising alphas/betas
	gpmodel.Kl = gausskernel(gpmodel.xtr, gpmodel.xtr, gpmodel.betaell,   gpmodel.alphaell,   gpmodel.tol);
	gpmodel.Ks = gausskernel(gpmodel.xtr, gpmodel.xtr, gpmodel.betasigma, gpmodel.alphasigma, gpmodel.tol);
	gpmodel.Ko = gausskernel(gpmodel.xtr, gpmodel.xtr, gpmodel.betaomega, gpmodel.alphaomega, gpmodel.tol);
	gpmodel.D = pdist2(gpmodel.xtr,gpmodel.xtr).^2;	
	
	if gpmodel.white
		% transform ell's into cholesky domain
		if ismember('l', nsfuncs)
			gpmodel.Ll = chol(gpmodel.Kl)';
			gpmodel.ell = gpmodel.Ll \ gpmodel.ell;
		end
		if ismember('s', nsfuncs)
			gpmodel.Ls = chol(gpmodel.Ks)';
			gpmodel.sigma = gpmodel.Ls \ gpmodel.sigma;
		end
		if ismember('o', nsfuncs)
			gpmodel.Lo = chol(gpmodel.Ko)';
			gpmodel.omega = gpmodel.Lo \ gpmodel.omega;
		end
	end
	
	
	samples = [];
	switch optim
		% HMC sampling, do 1 chain (call this function repeatedly for more chains)
		case 'hmc'
			[samples,mll] = nsgpnuts(gpmodel);
		
		% gradient descent (several restarts included)
		case 'grad'
			[gpmodel,mll] = nsgpgrad(gpmodel);
	end
	
	[mse,nmse,nlpd] = testerrors(gpmodel, samples);
end





