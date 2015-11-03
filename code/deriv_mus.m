function [dl, ds, do] = deriv_mus(gp)
% derivative of the mean parameters against MLL

	% update if necessary
	if sum(ismember('ab', gp.nsfuncs))
		gp.Kl = gausskernel(gp.xtr,gp.xtr,gp.betaell,   gp.alphaell,   gp.tol);
		gp.Ks = gausskernel(gp.xtr,gp.xtr,gp.betasigma, gp.alphasigma, gp.tol);
		gp.Ko = gausskernel(gp.xtr,gp.xtr,gp.betaomega, gp.alphaomega, gp.tol);
	end
		
	dl = sum(gp.Kl\(gp.l_ell - gp.l_muell));
	ds = sum(gp.Ks\(gp.l_sigma - gp.l_musigma));
	do = sum(gp.Ko\(gp.l_omega - gp.l_muomega));	
end


