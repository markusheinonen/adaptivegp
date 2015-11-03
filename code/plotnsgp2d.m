function [] = plotnsgp2d(gp)
	% grid
	ng = 70;
	nl = 10;
	[X1,X2] = meshgrid(linspace(0,1,ng)',linspace(0,1,ng)');
	Xt = [X1(:) X2(:)];

% 	% grid
% 	mg = 80;
% 	xrng = max(gp.xtr) - min(gp.xtr);
% 	x1 = linspace(min(gp.xtr(:,1))-0.05, max(gp.xtr(:,1))+0.05, mg)';
% 	x2 = linspace(min(gp.xtr(:,2))-0.05, max(gp.xtr(:,2))+0.05, mg)';
% 	[X1,X2] = meshgrid(x1,x2);
% 	Xt = [X1(:) X2(:)];

%		[xtr,ytr,ft,ftstd,lt,st,ot,ftderiv,ftderivstd] = denormalise(gp, ft, ftstd, lt, st, ot);
%		[ell,sigma,omega] = latentchols(gp);
%		ell = exp(ell);
%		sigma = exp(sigma);
%		omega = exp(omega);
%		xtr = gp.xtr;
%		ytr = gp.ytr;

%		ell = exp(ell) * geomean(gp.xscale);
%		sigma = exp(sigma) * gp.yscale;
%		omega = exp(omega) * gp.yscale;
%		[xtr,ytr] = denormalise(gp);


	[ft,ftstd,lt,st,ot] = nsgpposterior(gp, Xt);

	% model at grid
%	ot = exp(gpposterior(pars.xtr, omega, Xt, pars.muomega, pars.betaomega, pars.alphaomega, pars.tol));
%	lt = exp(gpposterior(pars.xtr, ell,   Xt, pars.muell,   pars.betaell,   pars.alphaell,   pars.tol));
%	st = exp(gpposterior(pars.xtr, sigma, Xt, pars.musigma, pars.betasigma, pars.alphasigma, pars.tol));

	subplot(2,3,1); contour(X1,X2, reshape(yt,ng,ng),nl); title('E[f]'); colorbar;
	hold on; 
		scatter(gp.xtr(:,1), gp.xtr(:,2), 20, gp.ytr, 'filled'); 
	hold off;
	subplot(2,3,2); contour(X1,X2, reshape(diag(ytcov),ng,ng),nl); title('Var[f]'); colorbar;
	subplot(2,3,4); contour(X1,X2, reshape(lt,ng,ng),nl); title('ell'); colorbar;
	subplot(2,3,5); contour(X1,X2, reshape(st,ng,ng),nl); title('sigma'); colorbar;
	subplot(2,3,6); contour(X1,X2, reshape(ot,ng,ng),nl); title('omega'); colorbar;
end

