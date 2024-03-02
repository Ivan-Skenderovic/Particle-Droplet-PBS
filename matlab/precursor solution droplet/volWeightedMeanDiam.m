function D43 = volWeightedMeanDiam(diams, weights)

    D43 = sum(weights.*diams.^4)/sum(weights.*diams.^3);

end

