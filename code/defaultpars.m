function pars = defaultpars(xtr)


	pars = struct;
	
	n = size(xtr,1);
	pars.n = size(xtr,1);
	pars.d = size(xtr,2);
	
	% initial latent vectors
	pars.ell   = log( 0.1) * ones(n,1);
	pars.sigma = log( 0.3) * ones(n,1);
	pars.omega = log(0.05) * ones(n,1);
	
	% hyperparameters
	pars.betaell = 0.15;
	pars.betasigma = 0.2;
	pars.betaomega = 0.3;
	pars.alphaomega = 1;
	pars.alphasigma = 1;
	pars.alphaell = 1;
	pars.muell = log(0.1);
	pars.musigma = log(0.3);
	pars.muomega = log(0.05);
	pars.tol = 1e-3;
	pars.white = 1;
	
	% hmc parameters
	pars.warmups = 1000; % warmups are skipped with ns-lengthscale
	pars.iters = 200;
	pars.delta = 0.65;
	pars.maxdepth = 9;
	pars.leapfrog = 0.01;
	
	% gradient parameters
	pars.maxrounds = 3;
	pars.restarts = 10;
	pars.graditers = 5000;
	
	% verbosity parameters
	pars.plotiters = false;
	pars.verbose = true;
end


