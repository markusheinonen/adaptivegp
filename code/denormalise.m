function [xtr,ytr,ft,ftstd,lt,st,ot] = denormalise(pars, ft, ftstd, lt, st, ot)
	
%	if ~exist('ftderiv', 'var')
%		ftderiv = 0;
%	end
%	if ~exist('ftderivstd', 'var')
%		ftderivstd = 0;
%	end
	if ~exist('ot', 'var')
		ot = 0;
	end
	if ~exist('lt', 'var')
		lt = 0;
	end
	if ~exist('st', 'var')
		st = 0;
	end
	if ~exist('ftstd', 'var')
		ftstd = 0;
	end
	if ~exist('ft', 'var')
		ft = 0;
	end
	
	
	% denormalisation
	xtr = pars.xtr .* repmat(pars.xscale,pars.n,1) + repmat(pars.xbias,pars.n,1);
	ytr = pars.ytr * pars.yscale + pars.ybias;
	ft = ft * pars.yscale + pars.ybias;
	ftstd = ftstd * pars.yscale;
	lt = lt * geomean(pars.xscale);
	ot = ot * pars.yscale;
	st = st * pars.yscale;
	
%	ftderiv = ftderiv * pars.yscale + pars.ybias;
%	ftderivstd = ftderivstd * pars.yscale;
end
