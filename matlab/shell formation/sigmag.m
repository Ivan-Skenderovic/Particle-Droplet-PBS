function vals = sigmag(vals, weights)
vol_mean = sum(log(vals).*weights)./sum(weights);
vals = exp(sqrt(sum(weights.*((vol_mean - log(vals)).^2))./sum(weights)));
end

