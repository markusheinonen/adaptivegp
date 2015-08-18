function theta = makesample(pars)

	theta = [];
	if ismember('l', pars.nsfuncs); 
		theta = [theta pars.ell'];
	else
		theta = [theta pars.ell(1)];
	end
	if ismember('s', pars.nsfuncs); 
		theta = [theta pars.sigma'];
	else
		theta = [theta pars.sigma(1)];
	end
	if ismember('o', pars.nsfuncs); 
		theta = [theta pars.omega'];
	else
		theta = [theta pars.omega(1)];
%		theta = [theta pars.omega'];
	end
	if ismember('m', pars.nsfuncs); theta = [theta pars.muell pars.musigma pars.muomega]; end
	if ismember('a', pars.nsfuncs); theta = [theta pars.alphaell pars.alphasigma pars.alphaomega]; end
	if ismember('b', pars.nsfuncs); theta = [theta pars.betaell pars.betasigma pars.betaomega]; end

end




