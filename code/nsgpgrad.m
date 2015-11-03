function [gp,mll] = nsgpgrad(gp)
% learn a nonstationary GP 
	
	display(sprintf('Optimizing for %d restarts ...', gp.restarts));
	display('  model   iter stepsize       mll');
	
	% perform gradient search using 'gp.restarts' initial conditions
	% the first initial condition is always fixed to the 'gp.init_ell', etc
	% remaining are randomised
	
	% store models
	gps = cell(gp.restarts,1);
	mlls = zeros(gp.restarts,1);
	
	gps{1} = gradient(gp);
	mlls(1) = nsgpmll(gp);
	
	for iter=2:gp.restarts
		% set random initial values
		gp.init_ell = unifrnd(0.03, 0.3);
		gp.init_sigma = unifrnd(0.1, 0.5);
		gp.init_omega = unifrnd(0.01, 0.10);
		gp.init();  % do choleskies, etc.
		
		gps{iter} = gradient(gp);
		mlls(iter) = nsgpmll(gp);
	end
	
	[mll,i] = max(mlls);
	gp = gps{i};
	
	display(sprintf('Best model mll=%.2f',mll));

	
% 	
% 	gp = gradient(gp);
% 	mll = nsgpmll(gp);
% 	return;
% 
% 	bestgp = [];
% 	bestmll = -inf;
% 
% 	if strcmp(gp.nsfuncs,'lso')
% 		
% 		% directly go to full model
% 		gp = gradient(gp);
% 		mll = nsgpmll(gp);
% 		if mll > bestmll
% 			bestmll = mll;
% 			bestgp = gp;
% 		end
% 
% 		% iteratively add more stuff
% 		ords = perms(gp.nsfuncs);
% 		for j=1:size(ords,1)
% 			gp = gradient(gp, ''); % stationary
% 			for i=1:length(gp.nsfuncs)
% 				gp = gradient(gp, ords(j,1:i));
% 			end
% 			mll = nsgpmll(gp);
% 			
% 			if mll > bestmll
% 				bestmll = mll;
% 				bestgp = gp;
% 			end
% 		end
% 				
% 		% go from stationary to full
% 		gp = gradient(gp, '');
% 		gp = gradient(gp);
% 		mll = nsgpmll(gp);
% 		if mll > bestmll
% 			bestmll = mll;
% 			bestgp = gp;
% 		end
% 		
% 		gp = bestgp;
% 		mll = bestmll;
% 		
% 		return;
% 	end
% 
% 	
% 	
% 	%% find best standard GP solution from multiple restarts
% 	% generate initial values
% 	ls = log([0.10 unifrnd(0.02, 0.20, gp.restarts-1,1)']);
% 	ss = log([0.5*max(gp.ytr) unifrnd(0.20, 0.80, gp.restarts-1,1)']);
% 	os = log([mean(abs(diff(gp.ytr))) unifrnd(0.01, 0.20, gp.restarts-1,1)']);
% 	
% 	n = gp.n;
% 	
% 	bestmll = -inf;
% 	bestpars = gp;
% 	for r=1:gp.restarts
% %		display(sprintf(' [%d/%d] Gradients of stationary parameters', r, pars.restarts));
% 
% %		pars.muell = ls(r);
% %		pars.musigma = ss(r);
% %		pars.muomega = os(r);
% 		
% 		gp.ell = ls(r)*ones(n,1);
% 		gp.sigma = ss(r)*ones(n,1);
% 		gp.omega = os(r)*ones(n,1);
% 		
% 		if isfield(gp, 'Ll')
% 			gp.ell = gp.Ll\(ls(r)*ones(n,1));
% 		end
% 		if isfield(gp, 'Ls')
% 			gp.sigma = gp.Ls\(ss(r)*ones(n,1));
% 		end
% 		if isfield(gp, 'Lo')
% 			gp.omega = gp.Lo\(os(r)*ones(n,1));
% 		end
% 		
% 		gp = gradient(gp, '');
% 		
% 		mll = nsgpmll(gp);
% 	
% 		if mll > bestmll
% 			bestmll = mll;
% 			bestpars = gp;
% 		end
% 	end
% 	
% 	gp = bestpars;
% 	pars1 = gp;
% 	pars2 = gp;
% 	
% 	%% (1) add nonstationary functions one-by-one
% 	for k=1:pars1.maxrounds
% 		% optimise functions one-by-one in order (several rounds)
% 		for i=1:length(pars1.nsfuncs)
% %			display(sprintf(' [%d/%d] Gradients of non-stationary function "%s"', k, pars1.maxrounds, pars1.nsfuncs(i)));
% 			pars1 = gradient(pars1, pars1.nsfuncs(i));
% 		end
% 		
% 		% always do stationary learning at end
% %		display(' Gradients of stationary function');
% 		pars1 = gradient(pars1, '');
% 	end
% 	
% 	% full learning at end
% %	display(sprintf(' [1/1] Gradients of non-stationary function "%s"', pars1.nsfuncs));
% 	pars1 = gradient(pars1);
% 	mll1 = nsgpmll(pars1);
% 	
% 	%% (2) do all nonstationary functions simultaneously
% 	for k=1:pars2.maxrounds
% %		display(sprintf(' [%d/%d] Gradients of non-stationary function "%s"', k, pars2.maxrounds, pars2.nsfuncs));
% 		pars2 = gradient(pars2);
% 		
% 		% always do stationary learning at end
% %		display(' Gradients of stationary function');
% 		pars2 = gradient(pars2, '');
% 	end
% 	mll2 = nsgpmll(pars2);
% 	
% 	%% choose better
% 	if mll1 > mll2
% 		gp = pars1;
% 		mll = mll1;
% 	else
% 		gp = pars2;
% 		mll = mll2;
% 	end
% 	
% 	display(sprintf('Best model mll=%.2f',mll));
% 	
end




