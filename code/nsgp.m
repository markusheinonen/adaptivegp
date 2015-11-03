function [gp,samples,mll,mse,nmse,nlpd] = nsgp(x,y,nsfuncs,optim,varargin)
% Main function to learn the adaptive GP
% 
% INPUTS
%  - x       : input data vector
%  - y       : output data vector
%  - nsfuncs : string of components set as nonstationary [default 'lso']
%      l : lengthscale
%      s : signal variance
%      o : noise variance
%  - optim : either 'grad' [default] or 'hmc'
%  - varargin : define parameters
%               give params-struct or define pairs of values
%
% OUTPUTS
%  - params  : parameter values learned [grad] or initialised [hmc]
%  - samples : hmc samples of the posterior? [hmc] or empty [grad]
%

	if ~exist('nsfuncs','var')
		nsfuncs = 'lso';
	end
	if ~exist('optim','var')
		optim = 'grad';
	end
	
	% no params given
	if isempty(varargin)
		gp = gpmodel(x,y);
	end
	
	% single param or odd number given: assume first param is a struct
	if mod(length(varargin),2)==1
		gp = varargin{1};
		varargin = varargin(2:end);
	end
	
	% set all remaining arguments
	for i=1:length(varargin)/2
		gp.(varargin{i*2-1}) = varargin{i*2};
	end
	
	gp.nsfuncs = nsfuncs;
	gp.optim = optim;
	gp.init();
	
	samples = [];
	switch optim
		% HMC sampling, do 1 chain (call this function repeatedly for more chains)
		case 'hmc'
			[samples,mll] = nsgpnuts(gp);
		
		% gradient descent (several restarts included)
		case 'grad'
			[gp,mll] = nsgpgrad(gp);
	end
	
%	mse = nan;
%	nmse = nan;
%	nlpd = nan;
%	if cv < 1
%		[mse,nmse,nlpd] = testerrors(gp, samples);
%		display(sprintf('mse=%.4f nmse=%.4f mlpd=%.4f', mse, nmse, nlpd));
%	end
end





