function dl = deriv_ell(pars, scalar)
% derivative of the ell latent function wrt MLL

	if ~exist('scalar','var')
		scalar = 0;
	end
	
	[ell,sigma,omega] = latentchols(pars);
	
	if sum(ismember('ab', pars.nsfuncs))
		pars.Kl = gausskernel(pars.xtr, pars.xtr, pars.betaell, pars.alphaell, pars.tol);
	end
	
	n = length(pars.xtr);
	Ky = nsgausskernel(pars.xtr,pars.xtr,ell, ell, sigma, sigma, omega);
	a = Ky\pars.ytr;
	A = a*a' - inv(Ky);
		
	% correct in the log-transform?
	if length(pars.ell) == 1 || scalar
		Kf = nsgausskernel(pars.xtr, pars.xtr, ell, ell, sigma, sigma, log(0));
		dKl = exp(mean(ell))^(-2) * (pars.D .* Kf);

		dl = 0.5*(diag( A*dKl )) - (pars.Kl\(ell - pars.muell));
		
		dl = ones(n,1)*sum(dl);
		
		if isfield(pars, 'Ll')
			dl = pars.Ll\dl;
		end
		
	else
		ell = exp(ell);
		sigma = exp(sigma);

		% compute dK matrix first, then cross-slice from it
		L = repmat(ell.^2,1,n) + repmat(ell.^2,1,n)';
		E = exp(-pars.D./L);
		R = sqrt( 2*(ell*ell') ./ L );
		dK = (ell*ell') .* (sigma*sigma') .* E .* (R.^(-1)) .* (L.^(-3)) .* (4 * pars.D .* repmat(ell.^2,1,n) - repmat(ell.^4,1,n) + repmat(ell'.^4,n,1));

%		dK = zeros(n,n);		
%		
%		for i=1:n
%			li = ell(i);
%			for j=1:n
%				lj = ell(j);
%				lij = li^2 + lj^2;
%				eij = exp(-pars.D(i,j)/lij);
%				rij = sqrt( (2*li*lj)/(li^2 + lj^2) );
%				sij = sigma(i)*sigma(j);
%				dK(i,j) = li*lj * sij*eij*rij^(-1)*lij^(-3) * (4*pars.D(i,j)*li^2 - li^4 + lj^4);
%			end
%		end
		
		dl = zeros(n,1);
		for i=1:n
%			ei = zeros(n,1);
%			ei(i) = 1;
%			Mi = bsxfun(@or, ei, ei');
			
			% original code:
%			dl(i) = 0.5*sum(diag(A * (Mi .* dK) ));
			% optimised:
			% i) replace diag() with sum(,2) since we only need diagonal results
			% ii) make Mi.*dK sparse
%			dl(i) = 0.5*sum( sum(A .* sparse(Mi .* dK)',2) );
			% iii) construct sparse matrix directly from the col/row
			dl(i) = 0.5*sum( sum(A .* sparse([1:n i*ones(1,n)], [i*ones(1,n) 1:n], [dK(:,i) dK(i,:)'])',2));
		end
		
		
		dl = dl - pars.Kl\(log(ell)-pars.muell);
%		dl = dl - Kl\(log(ell)-mean(log(ell)));
		
		% transform
		if isfield(pars, 'Ll')
			dl = pars.Ll'*dl;
		end
	end
end


