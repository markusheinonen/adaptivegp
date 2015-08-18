function nmse = computenmse(ypred, yts, ytr)
	nmse = mean((yts - ypred).^2) / mean((yts - mean(ytr)).^2);
end



