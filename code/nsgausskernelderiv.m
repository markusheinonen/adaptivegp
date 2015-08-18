function Kd = nsgausskernelderiv(x1, x2, l1, l2, s1, s2)
% ns-gauss kernel's first derivative
		
	sumcross = @(v1,v2) repmat(v1,1,length(v2)) + repmat(v2',length(v1),1);
	Kf = nsgausskernel(x1,x2,l1,l2,s1,s2,log(0));
	
	D = repmat(x1,1,length(x2)) - repmat(x2,1,length(x1))';
	L = sumcross(exp(l1).^2,exp(l2).^2);
	
	Kd = - 2 * D .* L.^(-1) .* Kf;
end

