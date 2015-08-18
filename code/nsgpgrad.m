function [pars,mll] = nsgpgrad(pars)
% learn a nonstationary GP by incrementing nonstationarity
% i.e. start from stationary solution, add nonstationarities, etc.
	
	display('Optimizing...');
	display('  model   iter stepsize       mll');


	if strcmp(pars.nsfuncs,'lso')
		pars1 = gradient(pars); % full
		mll1 = nsgpmll(pars1);
	
		pars2 = gradient(pars, ''); % stationary
		pars2 = gradient(pars2);  % full
		mll2 = nsgpmll(pars2);
	
		mll3 = -inf;
		pars3 = [];
		if pars.nsfuncs
			pars3 = gradient(pars, ''); % stationary
			for k=1:2
				for i=1:length(pars.nsfuncs) % singletons
					pars3 = gradient(pars3, pars.nsfuncs(i));
				end
			end
			pars3 = gradient(pars3);  % full
			mll3 = nsgpmll(pars3);
		end
		
		if mll1 >= max(mll2,mll3)
			pars = pars1;
			mll = mll1;
		end
		if mll2 >= max(mll1,mll3)
			pars = pars2;
			mll = mll2;
		end
		if mll3 >= max(mll1,mll2)
			pars = pars3;
			mll = mll3;
		end	
		return;
	end

	
	
	%% find best standard GP solution from multiple restarts
	% generate initial values
	ls = log([0.10 unifrnd(0.02, 0.20, pars.restarts-1,1)']);
	ss = log([0.5*max(pars.ytr) unifrnd(0.20, 0.80, pars.restarts-1,1)']);
	os = log([mean(abs(diff(pars.ytr))) unifrnd(0.01, 0.20, pars.restarts-1,1)']);
	
	n = pars.n;
	
	bestmll = -inf;
	bestpars = pars;
	for r=1:pars.restarts
%		display(sprintf(' [%d/%d] Gradients of stationary parameters', r, pars.restarts));

		pars.muell = ls(r);
		pars.musigma = ss(r);
		pars.muomega = os(r);
		
		pars.ell = ls(r)*ones(n,1);
		pars.sigma = ss(r)*ones(n,1);
		pars.omega = os(r)*ones(n,1);
		
		if isfield(pars, 'Ll')
			pars.ell = pars.Ll\(ls(r)*ones(n,1));
		end
		if isfield(pars, 'Ls')
			pars.sigma = pars.Ls\(ss(r)*ones(n,1));
		end
		if isfield(pars, 'Lo')
			pars.omega = pars.Lo\(os(r)*ones(n,1));
		end
		
		pars = gradient(pars, '');
		
		mll = nsgpmll(pars);
	
		if mll > bestmll
			bestmll = mll;
			bestpars = pars;
		end
	end
	
	pars = bestpars;
	pars1 = pars;
	pars2 = pars;
	
	%% (1) add nonstationary functions one-by-one
	for k=1:pars1.maxrounds
		% optimise functions one-by-one in order (several rounds)
		for i=1:length(pars1.nsfuncs)
%			display(sprintf(' [%d/%d] Gradients of non-stationary function "%s"', k, pars1.maxrounds, pars1.nsfuncs(i)));
			pars1 = gradient(pars1, pars1.nsfuncs(i));
		end
		
		% always do stationary learning at end
%		display(' Gradients of stationary function');
		pars1 = gradient(pars1, '');
	end
	
	% full learning at end
%	display(sprintf(' [1/1] Gradients of non-stationary function "%s"', pars1.nsfuncs));
	pars1 = gradient(pars1);
	mll1 = nsgpmll(pars1);
	
	%% (2) do all nonstationary functions simultaneously
	for k=1:pars2.maxrounds
%		display(sprintf(' [%d/%d] Gradients of non-stationary function "%s"', k, pars2.maxrounds, pars2.nsfuncs));
		pars2 = gradient(pars2);
		
		% always do stationary learning at end
%		display(' Gradients of stationary function');
		pars2 = gradient(pars2, '');
	end
	mll2 = nsgpmll(pars2);
	
	%% choose better
	if mll1 > mll2
		pars = pars1;
		mll = mll1;
	else
		pars = pars2;
		mll = mll2;
	end
	
	display(sprintf('Best model mll=%.2f',mll));
	
end




