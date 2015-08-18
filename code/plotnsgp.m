function [] = plotnsgp(pars, plotlatent, plotderivs, plotkernel, ts, truemodel)

	if ~exist('truemodel','var')
		truemodel = [];
	end
	if ~exist('plotkernel','var')
		plotkernel = false;
	end
	if ~exist('plotderivs','var')
		plotderivs = false;
	end
	if ~exist('plotlatent','var')
		plotlatent = false;
	end
	if ~exist('ts','var')
		ts = false;
	end
	
	%% 2D plot
	if size(pars.xtr,2) == 2
		% grid
		mg = 70;
		x1 = linspace(min(pars.xtr(:,1))*1.05,max(pars.xtr(:,1))*1.05,mg)';
		x2 = linspace(min(pars.xtr(:,2))*1.05,max(pars.xtr(:,2))*1.05,mg)';
		[X1,X2] = meshgrid(x1,x2);
		Xt = [X1(:) X2(:)];
		
		yt = nsgpposterior(pars, pars.xtr);
		[ft,ftstd,lt,st,ot] = nsgpposterior(pars, Xt);
				
		subplot(2,3,1); contour(X1,X2, reshape(ft,mg,mg),10); title('E[f]'); colorbar;
		hold on; 
			scatter(pars.xtr(:,1), pars.xtr(:,2), 20, pars.ytr, 'filled'); 
		hold off;
		
		subplot(2,3,2); scatter(pars.xtr(:,1), pars.xtr(:,2), 20, abs(pars.ytr-yt), 'filled'); colorbar;
%		subplot(2,3,2); contour(X1,X2, reshape(ftstd.^2,mg,mg),10); title('Var[f]'); colorbar;
		subplot(2,3,4); contour3(X1,X2, reshape(lt,mg,mg),10); title('ell'); colorbar;
%		hold on;
%			scatter(pars.xtr(:,1), pars.xtr(:,2), 20, pars.ytr, 'filled'); 
%		hold off;
		subplot(2,3,5); contour3(X1,X2, reshape(st,mg,mg),10); title('sigma'); colorbar;
%		hold on; 
%			scatter(pars.xtr(:,1), pars.xtr(:,2), 20, pars.ytr, 'filled'); 
%		hold off;
		subplot(2,3,6); contourf(X1,X2, reshape(ot,mg,mg),10); title('omega'); colorbar;
%		hold on; 
%			scatter(pars.xtr(:,1), pars.xtr(:,2), 20, pars.ytr, 'filled'); 
%		hold off;
		return;
	end
	

	
	squares = 1 + plotkernel + plotlatent + plotderivs;
	cols = [248 118 109; 0 186 56; 97 156 255] / 255; % ggplot2 colors
	xt = linspace(0,1,200)';
	
	% MLL of the initial values
	nsmll = nsgpmll(pars);
	
	% do inverse cholesky transformation
	
	[ft,ftstd,lt,st,ot] = nsgpposterior(pars,xt);

	if plotderivs
		[~,~,~,~,~,ftderiv,ftstdderiv] = nsgpposterior(pars,xt);
		ftderiv = ftderiv * pars.yscale;
		ftderivstd = ftstdderiv * pars.yscale;
	end
	
	% normalisation
	xtr = pars.xtr * pars.xscale + pars.xbias;
	ytr = pars.ytr * pars.yscale + pars.ybias;
	xt = xt * pars.xscale + pars.xbias;
	ft = ft * pars.yscale + pars.ybias;
	ftstd = ftstd * pars.yscale;
	lt = lt * pars.xscale;
	ot = ot * pars.yscale;
	st = st * pars.yscale;
	
	% test predictions
	if ts
		[mse,~,nlpd] = testerrors(pars);
	end

	%% function posterior
	if squares > 1
		subplot(squares,1,1);
	end
	
	% plot 1D curve
	h = fill([xt; flip(xt)], [ft+2*ftstd; flip(ft-2*ftstd)], 'black','facealpha', 0.25);
	set(h, 'EdgeColor','none');
	hold on;
		h1 = plot(xtr,ytr,'ko');
		h2 = plot(xt, ft, 'k-', 'linewidth',1.5);
		h3 = plot(xt, ft + 2*sqrt(ftstd.^2 + ot.^2), 'k--');

		% plot test data
		if ts
			plot(pars.xts, pars.yts, 'ro');
		end

		if ~isempty(truemodel)
			h4 = plot(truemodel.x, truemodel.f, 'r-', 'linewidth',1.0);
%			plot(truemodel.x, truemodel.f + 2*truemodel.sd, 'r--', 'linewidth',1.5);
%			plot(truemodel.x, truemodel.f - 2*truemodel.sd, 'r--', 'linewidth',1.5);
		end
		plot(xt, ft - 2*sqrt(ftstd.^2 + ot.^2), 'k--');
	hold off;
	if ts
		title(sprintf('Nonstationary function 95%% posterior [MLL=%.3f, MSE=%.3f, NLPD=%.3f]', nsmll,mse,nlpd));
	else
%		title(sprintf('Nonstationary function 95%% posterior [MLL  %.3f]', nsmll));
		if strcmp(pars.nsfuncs,'')
			title('Stationary GP function 95% posterior');
		else
			title(sprintf('Nonstationary %s-GP function 95%% posterior', upper(pars.nsfuncs)));
		end
	end
	xlabel('time');
	ylabel('value');
	legend('95% posterior','training data','mean posterior','95% noisy posterior');
	if ~isempty(truemodel)
		legend([h1 h2 h h3 h4], {'Data','Posterior mean','Posterior','Noisy posterior', 'True function'});
	end
	set(gca,'fontsize',12)

	%% latent functions
	if plotlatent		
		% latent functions
		subplot(squares,1,2);
		plot(xt, lt, 'color',cols(1,:));
		hold on;
			plot(xt, st, 'color',cols(2,:));
			plot(xt, ot, 'color',cols(3,:));
			
			if ~isempty(truemodel)
				plot(truemodel.x, truemodel.l, '--', 'color',cols(1,:));
				plot(truemodel.x, truemodel.s, '--', 'color',cols(2,:));
				plot(truemodel.x, truemodel.o, '--', 'color',cols(3,:));
			end

%			plot(xtr, exp(ell),   '.', 'color',cols(1,:));
%			plot(xtr, exp(sigma), '.', 'color',cols(2,:));
%			plot(xtr, exp(omega), '.', 'color',cols(3,:));
		hold off;
		xlabel('time');
		ylabel('value');
		ylim([0,max([max(st),max(lt),max(ot)])*1.1]);

		legend('Lengthscale','Signal variance','Noise variance');
		
%		legend(sprintf('lengthscale (ell)  [\\mu=%.2f \\alpha=%.2f \\beta=%.2f]',       exp(pars.muell),   pars.alphaell,   pars.betaell),...
%			   sprintf('signal variance (sigma)  [\\mu=%.2f \\alpha=%.2f \\beta=%.2f]', exp(pars.musigma), pars.alphasigma, pars.betasigma),...
%			   sprintf('noise std (omega)  [\\mu=%.2f \\alpha=%.2f \\beta=%.2f]',       exp(pars.muomega), pars.alphaomega, pars.betaomega));
		if ~isempty(truemodel)
			legend('Lengthscale','Signal variance','Noise variance','Generating lengthscale','Generating signal variance','Generating noise variance');
		end
		title('Latent functions');
		
		set(gca,'fontsize',12)
	end
	
	
	if plotderivs
		subplot(squares, 1, 3);
		
		fill([xt; flip(xt)], [ftderiv + 2*ftstdderiv; flip(ftderiv - 2*ftstdderiv) ], 'black', 'facealpha', 0.25, 'edgecolor', 'none');
		hold on;
			plot(xt, ftderiv, 'k-');
			plot([min(xt),max(xt)], [0 0], 'k--');
		hold off;
		
		title('Derivative function 95% posterior');
		xlabel('time');
		ylabel('value');
	end
	
	%% derivatives
% 	if plotderivs
% 		subplot(squares,1,3);
% 		
% 		dl = deriv_ell(pars,   ~ismember('l', pars.nsfuncs));
% 		ds = deriv_sigma(pars, ~ismember('s', pars.nsfuncs));
% 		do = deriv_omega(pars, ~ismember('o', pars.nsfuncs));
% 		
% 		if isfield(pars, 'Ll')
% 			dl = pars.Ll*dl;
% 		end
% 		if isfield(pars, 'Ls')
% 			ds = pars.Ls*ds;
% 		end
% 		if isfield(pars, 'Lo')
% 			do = pars.Lo*do;
% 		end
% 		
% 		stem(xtr, do, 'color',cols(3,:));
% 		hold on;
% 			stem(xtr, ds, 'color',cols(2,:));
% 			stem(xtr, dl, 'color',cols(1,:));
% 		hold off;
% %		ylim('auto');
% %		legend('lengthscale (ell)','signal variance (sigma)','noise std (omega)');
% 		title('Latent derivatives');
% 		set(gca,'fontsize',12)
% 	end
	
	%% kernel
	if plotkernel
		subplot(squares,1,squares);
		K = nsgausskernel(xt,xt,log(lt/pars.yscale),log(lt/pars.yscale),log(st/pars.yscale),log(st/pars.yscale),log(ot/pars.yscale));
		imagesc(real(K));
		title('Kernel matrix');
		set(gca,'fontsize',12)
	end
end




