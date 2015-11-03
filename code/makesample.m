function theta = makesample(gp)

	theta = [];
	if ismember('l', gp.nsfuncs); 
		theta = [theta gp.wl_ell'];
	else
		theta = [theta gp.wl_ell(1)];
	end
	if ismember('s', gp.nsfuncs); 
		theta = [theta gp.wl_sigma'];
	else
		theta = [theta gp.wl_sigma(1)];
	end
	if ismember('o', gp.nsfuncs); 
		theta = [theta gp.wl_omega'];
	else
		theta = [theta gp.wl_omega(1)];
	end
	
	if ismember('m', gp.nsfuncs); 
		theta = [theta gp.l_muell gp.l_musigma gp.l_muomega]; 
	end
	if ismember('a', gp.nsfuncs); theta = [theta gp.alphaell gp.alphasigma gp.alphaomega]; end
	if ismember('b', gp.nsfuncs); theta = [theta gp.betaell gp.betasigma gp.betaomega]; end

end




