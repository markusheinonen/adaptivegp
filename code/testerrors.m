function [mse,nmse,nlpd] = testerrors(gp, samples)

	if ~exist('samples', 'var')
		samples = [];
	end

	m = size(samples,1);
%	yts = pars.yts * pars.yscale + pars.ybias;
%	ytr = pars.ytr * pars.yscale + pars.ybias;
	
	nmse = 0;

	if samples
		mse = zeros(m,1);
%		nmse = zeros(m,1);
		nlpd = zeros(m,1);

		for i=1:m
			gp = parsesample(samples(i,:), gp);
			[mts,stds,~,~,ots] = nsgpposterior(gp, gp.xts);
			
			% compute over samples
			mse(i) = computemse(mts, yts);
%			nmse(i) = computenmse(mts, yts, ytr);
			nlpd(i) = computenlpd(mts, yts, stds.^2 + ots.^2);
		end
	else
		% compute posterior for test points
		[ft,stds,~,~,ot] = nsgpposterior(gp, gp.xts);
				
		% compute statistics: MLL, NMSE, NLPD
		mse = computemse(ft, gp.yts);
%		nmse = computenmse(ft, gp.yts, gp.ytr);
		nlpd = computenlpd(ft, gp.yts, stds.^2 + ot.^2);
	end
end



