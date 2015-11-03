function nlpd = computenlpd(ypred, yts, var)

	nlpd = zeros(length(yts),1);
	for i=1:length(yts)
		nlpd(i) = logmvnpdf(yts(i), ypred(i), var(i));
	end
	nlpd = mean(nlpd);
%	nlpd = - nlpd / length(yts);
end

