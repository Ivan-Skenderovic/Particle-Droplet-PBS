function [binMidPoints, binWeights, binWidth] = normalize_psd( vals, weights, ...
    noBins, lowerBinBoundary, upperBinBoundary )

BinWeight = zeros(1, noBins);
binBoundaries = logspace( lowerBinBoundary, upperBinBoundary, noBins+1);    

  for binIdx = 1:noBins
      BinWeight(binIdx) = sum( ...
          weights( ...
                  vals > binBoundaries(binIdx) & ...
                  vals <= binBoundaries(binIdx+1) ...
                 ) ...
          );
  end
    
binMidPoints = (binBoundaries(2:end) + binBoundaries(1:end-1))/2;
binWidth = log(binBoundaries(2:end)) - log(binBoundaries(1:end-1));
binWeights = BinWeight./binWidth;
 
end

