function [] = plotnsgp(gp, plotlatent, plotderivs, plotkernel, ts, truemodel)
% Plots the GP model
% parameters:
%  - plotlatent : plot latent functions (ell,sigma,omega)
%  - plotderivs : plot derivative GP
%  - plotkernel : plot covariance matrix
%  - ts         : test data to include
%  - truemodel  : highlight true function
%

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
	if size(gp.xtr,2) == 2
		plotnsgp2d(gp);
		return;
	end
	
	%%
	
	squares = 1 + plotkernel + plotlatent + plotderivs;
	cols = [248 118 109; 0 186 56; 97 156 255] / 255; % ggplot2 colors
	xt = linspace(0,1,250)';
	
	% MLL of the initial values
	nsmll = nsgpmll(gp);
	
	[ft,ftstd,lt,st,ot] = nsgpposterior(gp,xt);
	xt = linspace(0,1,250)' * gp.xscale + gp.xbias;

	[xtr,ytr] = denormalise(gp);
%	xtr = gp.xtr;
%	ytr = gp.ytr;
	
	% test predictions
	if ts
		[mse,~,nlpd] = testerrors(gp);
	end

	%% function posterior
	if squares > 1
		subplot(squares,1,1);
	end
	
	% plot 1D curve	
	for i=1:gp.p
		h = fill([xt; flip(xt)], [ft(:,i)+2*ftstd; flip(ft(:,i)-2*ftstd)], 'black','facealpha', 0.45, 'edgecolor','none');
		hold on;
	end
	hold on;
		for i=1:gp.p
			h3 = fill([xt; flip(xt)], [ft(:,i)+2*ftstd; flip(ft(:,i)+2*sqrt(ftstd.^2+ot.^2))], 'black','facealpha', 0.2, 'edgecolor','none');
			fill([xt; flip(xt)], [ft(:,i)-2*ftstd; flip(ft(:,i)-2*sqrt(ftstd.^2+ot.^2))], 'black','facealpha', 0.2, 'edgecolor','none');
		end
		
		h1 = plot(xtr, ytr, 'o', 'color', 'black', 'MarkerSize', 6);
		plot(xtr, ytr, '.', 'color', 'white', 'MarkerSize', 12);

		if ~isempty(truemodel)
			h4 = plot(truemodel.x, truemodel.f, 'r-', 'linewidth',1.0);
		end
	hold off;
	if ts
		title(sprintf('Nonstationary function 95%% posterior [MLL=%.3f, MSE=%.3f, NLPD=%.3f]', nsmll,mse,nlpd));
	else
		if strcmp(gp.nsfuncs,'')
			title('Stationary GP function 95% posterior');
		else
			title(sprintf('Nonstationary %s-GP function 95%% posterior', upper(gp.nsfuncs)));
		end
	end
	xlabel('time');
	ylabel('value');
%	legend([h1 h2 h h3], {'Data','Posterior mean','95% posterior','95% noisy posterior'}, 'Location','NorthWest');
%	if ~isempty(truemodel)
%		legend([h1 h2 h h3 h4], {'Data','Posterior mean','Posterior','Noisy posterior', 'True function'}, 'Location','NorthWest');
%	end
	set(gca,'fontsize',12)

	%% latent functions
	if plotlatent		
		% latent functions
		subplot(squares,1,2);
		plot(xt, lt, 'color',cols(1,:));
		hold on;
			plot(xt, st, 'color',cols(2,:));
			plot(xt, ot, 'color',cols(3,:));
			
%			plot(xtr, exp(ell),   '.', 'color',cols(1,:));
%			plot(xtr, exp(sigma), '.', 'color',cols(2,:));
%			plot(xtr, exp(omega), '.', 'color',cols(3,:));
		hold off;
		xlabel('time');
		ylabel('value');
		ylim([0,max([max(st),max(lt),max(ot)])*1.1]);
		l = legend('$\ell$ lengthscale','$\sigma^2$ signal variance','$\omega^2$ noise variance', 'Location','NorthWest');
		set(l,'Interpreter','latex')
		
%		legend(sprintf('lengthscale (ell)  [\\mu=%.2f \\alpha=%.2f \\beta=%.2f]',       exp(gp.muell),   gp.alphaell,   gp.betaell),...
%			   sprintf('signal variance (sigma)  [\\mu=%.2f \\alpha=%.2f \\beta=%.2f]', exp(gp.musigma), gp.alphasigma, gp.betasigma),...
%			   sprintf('noise std (omega)  [\\mu=%.2f \\alpha=%.2f \\beta=%.2f]',       exp(gp.muomega), gp.alphaomega, gp.betaomega));
%		if ~isempty(truemodel)
%			legend('Lengthscale','Signal variance','Noise variance','Generating lengthscale','Generating signal variance','Generating noise variance');
%		end
		title('Parameters');
		
		set(gca,'fontsize',12)
	end
	
	%% derivatives
	if plotderivs
		subplot(squares,1,3);
		
		dwl_l = deriv_ell(gp,   ~ismember('l', gp.nsfuncs));
		dwl_s = deriv_sigma(gp, ~ismember('s', gp.nsfuncs));
		dwl_o = deriv_omega(gp, ~ismember('o', gp.nsfuncs));
		
		dl_l = gp.Ll*dwl_l;
		dl_s = gp.Ls*dwl_s;
		dl_o = gp.Lo*dwl_o;
		
		stem(xtr, dl_o, 'color',cols(3,:));
		hold on;
			stem(xtr, dl_s, 'color',cols(2,:));
			stem(xtr, dl_l, 'color',cols(1,:));
		hold off;
		title('Latent derivatives');
		set(gca,'fontsize',12)
	end
	
	%% kernel
	if plotkernel
		subplot(squares,1,squares);
		K = nsgausskernel(xt,xt,log(lt/gp.yscale),log(lt/gp.yscale),log(st/gp.yscale),log(st/gp.yscale),log(ot/gp.yscale));
		imagesc(real(K));
		title('Kernel matrix');
		set(gca,'fontsize',12)
	end
end


