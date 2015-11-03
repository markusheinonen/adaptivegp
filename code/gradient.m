function gp = gradient(gp, nsfuncs)
% Compute gradient descent over the marginal log likelihood (MLL) of the
% adaptiveGP against the 9 parameters and 3 latent functions

	if ~exist('nsfuncs','var')
		nsfuncs = gp.nsfuncs;
	end
		
	% initial step size
	step = 1e-5;

	% MLL of the initial values
	mlls = zeros(gp.graditers,gp.p);
	mlls(1,:) = nsgpmll(gp);

	if gp.verbose
		display(sprintf('%5sgp %6.d %8.2f %9.2f', nsfuncs, 1, log10(step), mean(mlls(1,:))));
	end
	
	if gp.plotiters
		plotnsgp(gp,1,1);
		drawnow;
	end
	
%	dml = 0; dms = 0; dmo = 0;
	dbl = 0; dbs = 0; dbo = 0;
	dal = 0; das = 0; dao = 0;

	% gradient steps over all parameters
	for iter=2:gp.graditers
		
		% derivatives of latent functions and means
		dwl_l =   deriv_ell(gp, ~ismember('l',nsfuncs));
		dwl_s = deriv_sigma(gp, ~ismember('s',nsfuncs));
		dwl_o = deriv_omega(gp, ~ismember('o',nsfuncs));
		
		% deriv means
%		if ismember('m', pars.nsfuncs)
%			[dml,dms,dmo] = deriv_mus(pars);
%		end
		if ismember('b', gp.nsfuncs)
			[dbl,dbs,dbo] = deriv_betas(gp);
		end
		if ismember('a', gp.nsfuncs)
			[dal,das,dao] = deriv_alphas(gp);
		end
		
		% save old parameters
		l_cp = gp.wl_ell;
		s_cp = gp.wl_sigma;
		o_cp = gp.wl_omega;
		
		% gradient steps
		gp.wl_ell   = gp.wl_ell   + step*dwl_l;
		gp.wl_omega = gp.wl_omega + step*dwl_o;
		gp.wl_sigma = gp.wl_sigma + step*dwl_s;
		
		% gradients of alpha
		gp.alphaell   = gp.alphaell   + step*dal;
		gp.alphasigma = gp.alphasigma + step*das;
		gp.alphaomega = gp.alphaomega + step*dao;
		
		% gradients of beta
		gp.betaell   = gp.betaell   + step*dbl;
		gp.betasigma = gp.betasigma + step*dbs;
		gp.betaomega = gp.betaomega + step*dbo;
		
		% compute MLL
		mlls(iter,:) = nsgpmll(gp);
		
		% update step
		if all(mlls(iter,:) < mlls(iter-1,:))   % if overshooting, go back and decrease step size
			gp.wl_ell = l_cp;
			gp.wl_sigma = s_cp;
			gp.wl_omega = o_cp;
			mlls(iter,:) = mlls(iter-1,:);
			
			step = 0.70 * step; % drop 25% if failing
		else
			step = 1.10 * step; % increase 5% if going nicely
		end
		
		if gp.verbose && mod(iter,100) == 0
			display(sprintf('%5sgp %6.d %8.2f %9.2f', nsfuncs, iter, log10(step), mean(mlls(iter,:))));
		end
		
		if (log10(step) < -7) || (iter > 50 && (mean(mlls(iter,:)-mlls(iter-30,:))) < 0.1) 
			display(sprintf('%5sgp %6.d %8.2f %9.2f', nsfuncs, iter, log10(step), mean(mlls(iter,:))));
			break;
		end
		
		if gp.plotiters && mod(iter,10) == 0
			plotnsgp(gp,1,1);
			drawnow;
		end
	end
end



