function mean_diam_g = weightedGeoMean(vals, weights)
	%for lognormal distributions the geometric mean diameter equals 
	%the count median diameter
    mean_diam_g = exp( ...
                        sum(log(vals).*weights)./sum(weights) ...
                      );
end

