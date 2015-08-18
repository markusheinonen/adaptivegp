function grad = computegrads(pars)

	n = length(pars.ell);

%	nd = ismember('l',pars.nsfuncs)*n + ...
%		 ismember('s',pars.nsfuncs)*n + ...
%		 ismember('o',pars.nsfuncs)*n + ...
%		 ismember('m',pars.nsfuncs)*3 + ...
%		 ismember('a',pars.nsfuncs)*3 + ...
%		 ismember('b',pars.nsfuncs)*3;
	

	nd = 1;

	% compute derivs
	j = 0;
	grad = zeros(1,nd);
	if ismember('l', pars.nsfuncs); 
		grad(j+1:j+n) = deriv_ell(pars); j=j+n; 
	else
		grad(j+1) = mean(deriv_ell(pars,1)); j=j+1;
	end
	if ismember('s', pars.nsfuncs); 
		grad(j+1:j+n) = deriv_sigma(pars); j=j+n; 
	else
		grad(j+1) = mean(deriv_sigma(pars,1)); j=j+1;
	end
	if ismember('o', pars.nsfuncs); 
		grad(j+1:j+n) = deriv_omega(pars); j=j+n; 
	else
		grad(j+1) = mean(deriv_omega(pars,1)); j=j+1;
	end
	if ismember('m', pars.nsfuncs); grad(j+1:j+3) = deriv_mus(pars); j=j+3; end
	if ismember('a', pars.nsfuncs); grad(j+1:j+3) = deriv_alphas(pars); j=j+3; end
	if ismember('b', pars.nsfuncs); grad(j+1:j+3) = deriv_betas(pars); end
end



