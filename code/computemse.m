function mse = computemse(ypred, yts)
	mse = mean((yts - ypred).^2);
end



