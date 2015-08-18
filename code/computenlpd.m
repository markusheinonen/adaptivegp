function nlpd = computenlpd(ypred, yts, var)

	nlpd = 0;
	for i=1:length(yts)
		nlpd = nlpd + logmvnpdf(yts(i), ypred(i), var(i));
	end
	nlpd = - nlpd / length(yts);
end

