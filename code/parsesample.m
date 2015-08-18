function pars = parsesample(theta, pars)

	% retrieve results from 'theta' into 'pars'
	
	n = length(pars.ell);
	j = 0;
	if ismember('l', pars.nsfuncs)
		pars.ell   = theta(j+1:j+n)'; 
		j=j+n;
	else
		pars.ell   = theta(j+1) * ones(n,1); 
		j=j+1;
	end
	if ismember('s', pars.nsfuncs)
		pars.sigma = theta(j+1:j+n)';
		j=j+n;
	else
		pars.sigma = theta(j+1) * ones(n,1); 
		j=j+1;
	end
	if ismember('o', pars.nsfuncs)
		pars.omega = theta(j+1:j+n)';
		j=j+n;
	else
		pars.omega = theta(j+1) * ones(n,1); 
		j=j+1;
	end
	if ismember('m', pars.nsfuncs)
		pars.muell   = theta(j+1);
		pars.musigma = theta(j+2);
		pars.muomega = theta(j+3);
		j = j + 3;
	end
	if ismember('a', pars.nsfuncs)
		pars.alphaell   = theta(j+1);
		pars.alphasigma = theta(j+2);
		pars.alphaomega = theta(j+3);
		j = j + 3;
	end
	if ismember('b', pars.nsfuncs) 
		pars.betaell   = theta(j+1);
		pars.betasigma = theta(j+2);
		pars.betaomega = theta(j+3);
	end
end

