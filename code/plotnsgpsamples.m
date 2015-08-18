function [] = plotnsgpsamples(pars, samples, plotlatent, n, mapgp, truemodel)

	if ~exist('plotlatent','var')
		plotlatent = false;
	end
	if ~exist('n','var')
		n = 100;
	end
	if ~exist('truemodel','var')
		truemodel = [];
	end
	if ~exist('mapgp','var')
		mapgp = [];
	end

	squares = 2 + plotlatent*3;
	cols = [248 118 109; 0 186 56; 97 156 255] / 255;
	
	% MLL of the initial values
	xt = linspace(0,1,200)';
	nt = length(xt);
	m = size(samples,1);

	%% first check mll's
	mlls = zeros(n,1);
	for i=1:m
		pars = parsesample(samples(i,:), pars);
		mlls(i) = nsgpmll(pars);
	end
	
	% threshold worst 10% away
	thr = quantile(mlls,0.10);
	samples = samples(mlls>thr,:);
	m = size(samples,1);

	%% choose indices from remaining
	if n < m
		I = randperm(m,n);
	else
		I = 1:m;
	end

	%% precompute all functions and latents over all samples
	ots = zeros(n,nt);
	lts = zeros(n,nt);
	sts = zeros(n,nt);
	stds = zeros(n,nt);
	means = zeros(n,nt);
	for i=1:n
		pars = parsesample(samples(I(i),:), pars);
		[means(i,:),stds(i,:),lts(i,:),sts(i,:),ots(i,:)] = nsgpposterior(pars, xt);
	end
		
%	xtr = pars.xtr;
%	ytr = pars.ytr;

	%% map GP
	if ~isempty(mapgp)
		[ftmap,~,ltmap,stmap,otmap] = nsgpposterior(mapgp, xt);
		ftmap = ftmap * mapgp.yscale + mapgp.ybias;
		ltmap = ltmap * mapgp.yscale;
		stmap = stmap * mapgp.yscale;
	end

	% scale
	xt = xt * pars.xscale + pars.xbias;
	means = means * pars.yscale + pars.ybias;
	stds = stds * pars.yscale;
	ots = ots * pars.yscale;
	sts = sts * pars.yscale;
	xtr = pars.xtr * pars.xscale + pars.xbias;
	ytr = pars.ytr * pars.yscale + pars.ybias;
	
	
	%% plot function
	falpha = 0.65;
	
	if squares > 1
		subplot(squares, 1, [1 2]);
	end
	
	% plot posterior samples
	hold on;
	for i=1:n
		yt1 = means(i,:)'+2*sqrt(stds(i,:)'.^2 + ots(i,:)'.^2);
		yt2 = means(i,:)'-2*sqrt(stds(i,:)'.^2 + ots(i,:)'.^2);
		l(1) = fill([xt; flip(xt)], [yt1; flip(yt2)], 'black', 'facealpha', falpha/n, 'edgecolor', 'none');
	end
	
	lid = 2;
	if ~isempty(mapgp)
		lid = lid + 1;
		l(lid) = plot(xt, ftmap, 'k--', 'linewidth', 1.5);
	end
	if ~isempty(truemodel)
		lid = lid + 1;
		l(lid) = plot(truemodel.x, truemodel.f, 'r-');
	end
	
	l(2) = plot(xtr, ytr, 'o', 'color', 'black', 'MarkerSize', 6);
	plot(xtr, ytr, '.', 'color', 'white', 'MarkerSize', 12);

	if ~isempty(mapgp) && ~isempty(truemodel)
		legend(l, {'Posterior samples','Data','MAP posterior','True function'});
	elseif ~isempty(mapgp)
		legend(l, {'Posterior samples','Data','MAP posterior'});
	elseif ~isempty(truemodel)
		legend(l, {'Posterior samples','Data','True function'});
	else
		legend(l, {'Posterior samples','Data'});
	end
	
	ylim( [quantile(min(means - 2*sqrt(stds.^2 + ots.^2),[],2), 0.05), quantile(max(means + 2*sqrt(stds.^2 + ots.^2),[],2), 0.95)] );
	ylabel('value');
	
	if strcmp(pars.nsfuncs,'')
		title('Stationary GP function 95% posterior');
	else
		title(sprintf('Nonstationary %s-GP function 95%% posterior', upper(pars.nsfuncs)));
	end

	
	%% plot latent
	lw_s = 0.020;
	lw_l = 0.010;
	lw_o = 0.005;
	alpha_s = 5;
	alpha_l = 4;
	alpha_o = 2;

	if plotlatent		
		subplot(squares,1,3); hold on; lid = 1;
		for i=1:n
			h(lid) = fill([xt; flip(xt)], [lts(i,:)' + lw_l; flip(lts(i,:)' - lw_l)], cols(1,:), 'facealpha', alpha_l/n, 'edgecolor','none');
		end
		if ~isempty(mapgp)
			lid = lid + 1;
			h(lid) = plot(xt, ltmap, '--', 'color', cols(1,:), 'linewidth',2);
		end
		if ~isempty(truemodel)
			lid = lid + 1;
			h(lid) = plot(truemodel.x, truemodel.l, 'k-');
		end

		if ~isempty(mapgp) && ~isempty(truemodel)
			legend(h, {'Samples','MAP lengthscale', 'True lengthscale'});
		elseif ~isempty(truemodel)
			legend(h, {'Samples','True lengthscale'});
		elseif ~isempty(mapgp)
			legend(h, {'Samples','MAP lengthscale'});
		else
			legend(h, {'Samples'});
		end

		hold off;
		ylabel('value');
		ylim([0 quantile(max(lts'),0.99)]);
		title('Lengthscale posterior');
		
		subplot(squares,1,4); hold on; lid = 1;
		for i=1:n
			h(lid) = fill([xt; flip(xt)], [ots(i,:)' + lw_o; flip(ots(i,:)' - lw_o)], cols(3,:), 'facealpha', alpha_o/n, 'edgecolor','none');
		end
		if ~isempty(mapgp)
			lid = lid + 1;
			h(lid) = plot(xt, otmap, '--', 'color', cols(3,:), 'linewidth',2);
		end
		if ~isempty(truemodel)
			lid = lid + 1;
			h(lid) = plot(truemodel.x, truemodel.o, 'k-');
		end

		if ~isempty(mapgp) && ~isempty(truemodel)
			legend(h, {'Samples','MAP noise variance', 'True noise variance'});
		elseif ~isempty(truemodel)
			legend(h, {'Samples','True noise variance'});
		elseif ~isempty(mapgp)
			legend(h, {'Samples','MAP noise variance'});
		else
			legend(h, {'Samples'});
		end
		hold off;
%		xlabel('time');
		ylabel('value');
		ylim([0 quantile(max(ots'),0.99)]);
		title('Noise variance posterior');
		
		subplot(squares,1,5); hold on; lid = 1;
		for i=1:n
			h(lid) = fill([xt; flip(xt)], [sts(i,:)' + lw_s; flip(sts(i,:)' - lw_s)], cols(2,:), 'facealpha', alpha_s/n, 'edgecolor','none');
		end
		if ~isempty(mapgp)
			lid = lid + 1;
			h(lid) = plot(xt, stmap, '--', 'color', cols(2,:), 'linewidth',2);
		end
		if ~isempty(truemodel)
			lid = lid + 1;
			h(lid) = plot(truemodel.x, truemodel.s, 'k-');
		end

		if ~isempty(mapgp) && ~isempty(truemodel)
			legend(h, {'Samples','MAP signal variance', 'True variance'});
		elseif ~isempty(truemodel)
			legend(h, {'Samples','True signal variance'});
		elseif ~isempty(mapgp)
			legend(h, {'Samples','MAP signal variance'});
		else
			legend(h, {'Samples'});
		end
		
		hold off;
		ylabel('value');
		xlabel('time');
		ylim([0 quantile(max(sts'),0.99)]);
		title('Signal variance posterior');
	end
end




