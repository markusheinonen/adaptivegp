function [] = plotkernel(pars,xt)
	
	[~,~,lt,st,ot] = nsgpposterior(pars,xt);

	K = nsgausskernel(xt, xt, log(lt/pars.yscale), log(lt/pars.yscale), log(st/pars.yscale), log(st/pars.yscale), log(ot/pars.yscale));
	imagesc(real(K));
	title('Kernel matrix');
	set(gca,'fontsize',12)
end
