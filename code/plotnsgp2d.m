function [] = plotnsgp2d(pars)
	% grid
	ng = 70;
	nl = 10;
	[X1,X2] = meshgrid(linspace(0,1,ng)',linspace(0,1,ng)');
	Xt = [X1(:) X2(:)];

	[ell,sigma,omega] = latentchols(pars);

	% model at grid
	[yt,ytcov] = nsgpposterior(pars, Xt);
	ot = exp(gpposterior(pars.xtr, omega, Xt, pars.muomega, pars.betaomega, pars.alphaomega, pars.tol));
	lt = exp(gpposterior(pars.xtr, ell,   Xt, pars.muell,   pars.betaell,   pars.alphaell,   pars.tol));
	st = exp(gpposterior(pars.xtr, sigma, Xt, pars.musigma, pars.betasigma, pars.alphasigma, pars.tol));

	subplot(2,3,1); contour(X1,X2, reshape(yt,ng,ng),nl); title('E[f]'); colorbar;
	hold on; 
		scatter(pars.xtr(:,1), pars.xtr(:,2), 20, pars.ytr, 'filled'); 
	hold off;
	subplot(2,3,2); contour(X1,X2, reshape(diag(ytcov),ng,ng),nl); title('Var[f]'); colorbar;
	subplot(2,3,4); contour(X1,X2, reshape(lt,ng,ng),nl); title('ell'); colorbar;
	subplot(2,3,5); contour(X1,X2, reshape(st,ng,ng),nl); title('sigma'); colorbar;
	subplot(2,3,6); contour(X1,X2, reshape(ot,ng,ng),nl); title('omega'); colorbar;
end
		