function val = gpmll(x,y,pars)
	val = logmvnpdf(y, zeros(length(x),1), gausskernel(x,x,mean(pars.ell), mean(pars.sigma), mean(pars.eps)));
end



