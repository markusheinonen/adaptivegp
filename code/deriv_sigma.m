function dwl_s = deriv_sigma(gp, scalar)
% derivative of the sigma latent function wrt MLL
	
	if ~exist('scalar','var')
		scalar = 0;
	end

	n = length(gp.xtr);
	
	if sum(ismember('ab', gp.nsfuncs))
		gp.Ks = gausskernel(gp.xtr,gp.xtr,gp.betasigma, gp.alphasigma, gp.tol);
	end
	
	Ky = nsgausskernel(gp.xtr, gp.xtr, gp.l_ell, gp.l_ell, gp.l_sigma, gp.l_sigma, gp.l_omega);
	Kf = nsgausskernel(gp.xtr, gp.xtr, gp.l_ell, gp.l_ell, gp.l_sigma, gp.l_sigma, log(0));
	
	a = Ky\gp.ytr;
	A = a*a' - inv(Ky);

	dl_s = 2*diag( A * Kf ) - gp.Ks\(gp.l_sigma - gp.l_musigma);
	
	if scalar
		dl_s = ones(n,1) * sum(dl_s);
	end
	
	if ismember('s', gp.nsfuncs)
		dwl_s = gp.Ls'*dl_s;
	else
		dwl_s = gp.Ls\dl_s;
	end
end


