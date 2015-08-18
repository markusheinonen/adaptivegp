function pars = gradient(pars, nsfuncs)
% Compute gradient descent over the marginal log likelihood (MLL) of the
% adaptiveGP against the 9 parameters and 3 latent functions

	if ~exist('nsfuncs','var')
		nsfuncs = pars.nsfuncs;
	end
		
	% initialise
	step = 0.001;

	% MLL of the initial values
	mlls = zeros(pars.graditers,1);
	mlls(1) = nsgpmll(pars);

	if pars.verbose
%		display(sprintf('  iter     1 step %6.2f mll %8.2f', log10(step), mlls(1)));
%		display(sprintf('%5sgp %6.d %8.2f %9.2f', nsfuncs, 1, log10(step), mlls(1)));
	end
	if pars.plotiters
		plotnsgp(pars,1,1);
		drawnow;
	end
	
	
	nsl = 0;
	nss = 0;
	nso = 0;
	if ismember('l', nsfuncs)
		nsl = 1;
	end
	if ismember('s', nsfuncs)
		nss = 1;
	end
	if ismember('o', nsfuncs)
		nso = 1;
	end
	
	dml = 0; dms = 0; dmo = 0;
	dbl = 0; dbs = 0; dbo = 0;
	dal = 0; das = 0; dao = 0;

	% gradient steps over all parameters
	for iter=2:pars.graditers
		
%		[ell,sigma,omega] = latentchols(pars);
%		subplot(1,4,1); scatter(pars.xtr(:,1), pars.xtr(:,2), 30, pars.ytr, 'filled'); colorbar;
%		subplot(1,4,2); scatter(pars.xtr(:,1), pars.xtr(:,2), 30, exp(ell), 'filled'); colorbar;
%		subplot(1,4,3); scatter(pars.xtr(:,1), pars.xtr(:,2), 30, exp(sigma), 'filled'); colorbar;
%		subplot(1,4,4); scatter(pars.xtr(:,1), pars.xtr(:,2), 30, exp(omega), 'filled'); colorbar;
		drawnow;
		
		% derivatives of latent functions and means
		dl = deriv_ell(pars,~nsl);
		ds = deriv_sigma(pars,~nss);
		do = deriv_omega(pars,~nso);
		
		
		% deriv means
%		if ismember('m', pars.nsfuncs)
			[dml,dms,dmo] = deriv_mus(pars);
%		end
		if ismember('b', pars.nsfuncs)
			[dbl,dbs,dbo] = deriv_betas(pars);
		end
		if ismember('a', pars.nsfuncs)
			[dal,das,dao] = deriv_alphas(pars);
		end
		
		% save old functions
		p = pars;
%		l = pars.ell;
%		s = pars.sigma;
%		o = pars.omega;
		
		% gradient steps
		pars.ell = pars.ell + step*dl;
		pars.omega = pars.omega + step*do;
		pars.sigma = pars.sigma + step*ds;
		
		% gradients of mu's
		pars.muell = pars.muell + step*dml;
		pars.musigma = pars.musigma + step*dms;
		pars.muomega = pars.muomega + step*dmo;
%		pars.muell = mean(pars.ell);
%		pars.musigma = mean(pars.sigma); 
%		pars.muomega = mean(pars.omega);
		
		% gradients of alpha
		pars.alphaell   = pars.alphaell   + step*dal;
		pars.alphasigma = pars.alphasigma + step*das;
		pars.alphaomega = pars.alphaomega + step*dao;
		
		% gradients of beta
		pars.betaell   = pars.betaell   + step*dbl;
		pars.betasigma = pars.betasigma + step*dbs;
		pars.betaomega = pars.betaomega + step*dbo;
		
		% compute MLL
		mlls(iter) = nsgpmll(pars);
		
		% update step
		if mlls(iter) < mlls(iter-1)   % if overshooting, go back and decrease step size
%			pars.ell = l;
%			pars.sigma = s;
%			pars.omega = o;
			pars = p;
			mlls(iter) = mlls(iter-1);
			
			step = 0.70 * step; % drop 25% if failing
		else
			step = 1.10 * step; % increase 5% if going nicely
		end
		
		if pars.verbose && mod(iter,100) == 0
%			display(sprintf('  iter %5.d step %6.2f mll %8.2f', iter, log10(step), mlls(iter)));
			display(sprintf('%5sgp %6.d %8.2f %9.2f', nsfuncs, iter, log10(step), mlls(iter)));
		end
		
		% TODO fails between negative and positive
		if (log10(step) < -7) || (iter > 50 && (mlls(iter)-mlls(iter-30)) < 0.1) 
			break;
		end
		
		if pars.plotiters && mod(iter,10) == 0
			plotnsgp(pars,1,1);
			drawnow;
		end
	end
end



