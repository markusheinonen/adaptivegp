function [mse,nmse,nlpd] = testerrors(pars, samples)

	if ~exist('samples', 'var')
		samples = [];
	end

	m = size(samples,1);
	yts = pars.yts * pars.yscale + pars.ybias;
	ytr = pars.ytr * pars.yscale + pars.ybias;
	
	if samples
		mse = zeros(m,1);
		nmse = zeros(m,1);
		nlpd = zeros(m,1);

		for i=1:m
			pars = parsesample(samples(i,:), pars);
			[mts,stds,~,~,ots] = nsgpposterior(pars, pars.xts);
			
			% compute over samples
			mse(i) = computemse(mts, yts);
			nmse(i) = computenmse(mts, yts, ytr);
			nlpd(i) = computenlpd(mts, yts, stds.^2 + ots.^2);
		end
	else
		% compute posterior for test points
		[ft,stds,~,~,ot] = nsgpposterior(pars, pars.xts);
		
		% compute statistics: MLL, NMSE, NLPD
		mse = computemse(ft, yts);
		nmse = computenmse(ft, yts, ytr);
		nlpd = computenlpd(ft, yts, stds.^2 + ot.^2);
	end
end



