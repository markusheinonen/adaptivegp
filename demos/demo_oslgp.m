
addpath('../data');
addpath('../code');

load datasets;

% fit fully ns-gp
nsgp = nsgp(Dlso.x, Dlso.y, 'lso', 'grad');
figure;
plotnsgp(nsgp,true);



% fit with HMC and plot?(takes a while)
[lgp,samples] = nsgp(Dl.x, Dl.y, 'l', 'hmc');
figure;
plotnsgpsamples(lgp,samples,true);
  
