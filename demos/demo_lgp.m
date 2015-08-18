
addpath('../data');
addpath('../code');

load datasets;

% fit stationary GP
gp = nsgp(Dl.x, Dl.y, '', 'grad');
figure; 
plotnsgp(gp,true);

% fit L-GP
lgp = nsgp(Dl.x, Dl.y, 'l', 'grad');
figure;
plotnsgp(lgp,true);

% fit fully ns-gp
nsgp = nsgp(Dl.x, Dl.y, 'lso', 'grad');
figure; 
plotnsgp(nsgp,true);



% fit with HMC and plot?(takes a while)
[lgp,samples] = nsgp(Dl.x, Dl.y, 'l', 'hmc');
figure;
plotnsgpsamples(lgp,samples,true);
  
