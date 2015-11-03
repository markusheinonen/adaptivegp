% only class in the package: a gp-model with all parameters
classdef gpmodel < handle
	properties
		% data
		n,d,p
		xtr,ytr
		xbias
		xscale
		ybias
		yscale
		D
		
		% kernel parameters
		init_ell = 0.05;
		init_sigma = 0.30;
		init_omega = 0.05;
		betaell  = 0.20;
		betasigma = 0.20;
		betaomega = 0.30;
		alphaell = 1;
		alphasigma = 1;
		alphaomega = 1;
		muf = 0;
		tol = 1e-3;
		wl_ell
		wl_sigma
		wl_omega % white-log vectors, keep these in memory

		% kernel matrices and choleskies
		Kl,Ks,Ko
		Ll,Ls,Lo
		
		% remaining parameters
		nsfuncs = 'lso';
		optim = 'grad';
		warmups = 300;
		iters = 250;
		delta = 0.65;
		maxdepth = 9;
		leapfrog = 0.02;
		restarts = 3;
		graditers = 5000;
		plotiters = false;
		verbose = true;
	end
	% properties to be computed on the fly
	properties (Dependent)
		l_ell
		l_sigma
		l_omega
		ell
		sigma
		omega
		l_muell
		l_musigma
		l_muomega
	end
	methods
		% constructor
		function gp = gpmodel(x,y)
			% dimensions
			gp.n = size(x,1);
			gp.d = size(x,2);
			gp.p = size(y,2);
			
			% normalise inputs to [0,1] and output to max(abs(y)) = 1
			gp.xbias = min(x);
			gp.xscale = max(x)-min(x);
			gp.yscale = max(abs(y-mean(y)));
			gp.ybias = mean(y);
			x = (x - repmat(gp.xbias,gp.n,1)) ./ repmat(gp.xscale,gp.n,1);
			y = (y - gp.ybias) / gp.yscale;
			gp.muf = mean(y);
			
			gp.xtr = x;
			gp.ytr = y;
			gp.D = pdist2(gp.xtr,gp.xtr).^2;
			gp.n = size(gp.xtr,1);
			
			gp.init();
		end
		
		% init function
		function gp = init(gp)
			% prior covariance matrices
			gp.Kl = gausskernel(gp.xtr, gp.xtr, gp.betaell,   gp.alphaell,   gp.tol);
			gp.Ks = gausskernel(gp.xtr, gp.xtr, gp.betasigma, gp.alphasigma, gp.tol);
			gp.Ko = gausskernel(gp.xtr, gp.xtr, gp.betaomega, gp.alphaomega, gp.tol);
			
			% cholesky decompositions of these
			gp.Ll = chol(gp.Kl)';
			gp.Ls = chol(gp.Ks)';
			gp.Lo = chol(gp.Ko)';
			
			% set initial parameters in white-log domain
			gp.wl_ell = gp.Ll \ (log(gp.init_ell) * ones(gp.n,1));
			gp.wl_sigma = gp.Ls \ (log(gp.init_sigma) * ones(gp.n,1));
			gp.wl_omega = gp.Lo \ (log(gp.init_omega) * ones(gp.n,1));
		end
				
		% all getters that do something
		function l = get.ell(gp)
			l = exp(gp.Ll*gp.wl_ell);
		end
		function s = get.sigma(gp)
			s = exp(gp.Ls*gp.wl_sigma);
		end
		function o = get.omega(gp)
			o = exp(gp.Lo*gp.wl_omega);
		end
		function ll = get.l_ell(gp)
			ll = gp.Ll*gp.wl_ell;
		end
		function ls = get.l_sigma(gp)
			ls = gp.Ls*gp.wl_sigma;
		end
		function lo = get.l_omega(gp)
			lo = gp.Lo*gp.wl_omega;
		end
		% mu is always fixed to mean of that parameters' values
		function l_muell = get.l_muell(gp)
			l_muell = mean(gp.l_ell);
		end
		function l_musigma = get.l_musigma(gp)
			l_musigma = mean(gp.l_sigma);
		end
		function l_muomega = get.l_muomega(gp)
			l_muomega = mean(gp.l_omega);
		end
	end
end




