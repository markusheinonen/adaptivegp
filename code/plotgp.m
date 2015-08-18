function [] = plotgp(x,y,xt,pars)

	% MLL of the initial values
	mll = gpmll(x,y,pars);
	
	xt = linspace(0, 25, 200)';
	nt = length(xt);
	[yt,ytcov] = gpposterior(x, y, xt, 0, pars.ell, pars.sigma, pars.eps, pars.tol);
	ytvar = diag(ytcov);
	ytstd = sqrt(ytvar);
	
	subplot(2,1,1);
	h = fill([xt; flip(xt)], [yt+2*ytstd; flip(yt-2*ytstd)], 'red','facealpha', 0.20);
	set(h, 'EdgeColor','none');
	hold on;
		plot(x,y,'ko');
		plot(xt,yt, 'k-', 'linewidth',1.5);
		plot(xt, yt + 2*sqrt(ytstd.^2 + pars.eps^2), 'k--');
		plot(xt, yt - 2*sqrt(ytstd.^2 + pars.eps^2), 'k--');		
		

%		xt2 = linspace(0, 25, 50)';
%		[yt2,ytcov2] = gpposterior(x, y, xt2, 0, pars.ell, pars.sigma, pars.eps, pars.tol);
%		Z = mvnrnd(yt2, ytcov2 + 0.0001*eye(50),5)';
%		plot(xt2,Z, '-');
	hold off;
	title(sprintf('Stationary GP posterior [MLL  %.3f]', mll));
	xlabel('time');
	ylabel('value');
	legend('95% posterior','training data','mean posterior','95% noisy posterior');	
	
	subplot(2,1,2);
	K = gausskernel(xt,xt,pars.ell,pars.sigma,0);
	Z = mvnrnd(zeros(nt,1), K,5)';
	plot(xt,Z, '-x');
	title('5 samples from the prior N(0,K)');	
end




