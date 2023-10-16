function [ prepared_vals, prepared_weights] = prepare_psd( vals, weights, noBins, lowerBinBoundary, upperBinBoundary )
[sorted_vals, order_vals] = sort(vals);
sorted_weights = weights(order_vals);

BinWeight=zeros(1,noBins);
binBoundaries = logspace( lowerBinBoundary, upperBinBoundary, noBins+1);    

  for binIdx=1:noBins
      BinWeight(binIdx)=sum(sorted_weights(sorted_vals>binBoundaries(binIdx)&sorted_vals<=binBoundaries(binIdx+1))) + 1e-1;
      Indices=(sorted_vals>binBoundaries(binIdx)&sorted_vals<=binBoundaries(binIdx+1));
  end
    
BinVals=(binBoundaries(2:end) + binBoundaries(1:end-1))/2;
BinNorm=(log(binBoundaries(2:end)) - log(binBoundaries(1:end-1)));

prepared_vals = BinVals;
prepared_weights = BinWeight./BinNorm;
 
end

