function [samples,mlls] = nsgpnuts(pars)
% HMC-NUTS-DA sampling with Hoffman/Gelman code
% hmc parameters:
%  pars.nsfuncs is a string of derivatives to include, defaults to 'sle'
%    where 'l' = ell   [lengthscale]
%          's' = sigma [signal]
%          'o' = omega [error]
%          'm' = prior means
%          'a' = prior alphas
%          'b' = prior betas
%  pars.warmups  : default 1000
%  pars.iters    : default 1000
%  pars.maxdepth : default 9
%  pars.delta    : default 0.5
%

	function [logp,grad] = gradients(theta)
		% put sample into pars
		pars = parsesample(theta,pars);
		
		% compute logp and grads
		logp = nsgpmll(pars);
		grad = computegrads(pars);
	
		if pars.plotiters && iter > 20 && mod(iter,100)==0
			plotnsgp(pars,1,1);
			drawnow;
		end
		
		% iter counter
		iter = iter + 1;
	end

	%% defaults if not specified before
	if ~isfield(pars, 'iters')
		pars.iters = 1000;
	end
	if ~isfield(pars, 'warmups')
		pars.warmups = 1000;
	end
	if ~isfield(pars, 'delta')
		pars.delta = 0.5;
	end
	if ~isfield(pars, 'nsfuncs')
		pars.nsfuncs = 'slo';
	end

	f = @gradients;
	iter = 1;
	
%	% take MAP optimal as starting point
%	parsstat = gradient(xtr, ytr, pars, '');
%	parsmap = gradient(xtr, ytr, parsstat);

	% make initial sample
	theta0 = makesample(pars);
	
	% compute HMC-NUTS
	if ismember('l', pars.nsfuncs)
		samples = nuts(pars.leapfrog, f, pars.iters, theta0);
	else
		samples = nuts_da(f, pars.iters, pars.warmups, theta0, pars.delta, pars.maxdepth, pars.verbose);
	end
	
	% compute MLLS
	mlls = zeros(pars.iters,1);
	for i=1:pars.iters
		pars = parsesample(samples(i,:), pars);
		mlls(i) = nsgpmll(pars);
	end
end





