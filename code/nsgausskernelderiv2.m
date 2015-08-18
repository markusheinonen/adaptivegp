function Kdd = nsgausskernelderiv2(x1, x2, l1, l2, s1, s2)
% ns-gauss kernel's second derivative
	
	sumcross = @(v1,v2) repmat(v1,1,length(v2)) + repmat(v2',length(v1),1);
	Kf = nsgausskernel(x1,x2,l1,l2,s1,s2,log(0));
	
	D = pdist2(x1,x2).^2;
	L = sumcross(exp(l1).^2,exp(l2).^2);
	
	Kdd = (2*L.^(-1) - 4*(D.*L.^(-2))) .* Kf;
end
