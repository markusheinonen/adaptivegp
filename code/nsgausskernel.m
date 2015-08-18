function K = nsgausskernel(x1, x2, l1, l2, s1, s2, o)
% non-stationary (scalar) gaussian kernel
	
	l1 = exp(l1);
	l2 = exp(l2);
	s1 = exp(s1);
	s2 = exp(s2);
	o = exp(o);
	
	sumcross = @(v1,v2) repmat(v1,1,size(v2,1)) + repmat(v2',size(v1,1),1);
	
	% loop through all inputs if variables are matrices
	K = s1*s2' .* sqrt((2*l1*l2')./sumcross(l1.^2,l2.^2)) .* exp(- (pdist2(x1,x2).^2 ./ sumcross(l1.^2,l2.^2)));
	K = K + diag(o.^2);
end

