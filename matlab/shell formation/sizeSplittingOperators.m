function splittingOperators = sizeSplittingOperators(gridVolume, noNodes)

splittingOperators = zeros(noNodes, noNodes, noNodes);

   for i = 3 : noNodes - 2
       for j = 3 : noNodes - 2
           for k = 3 : noNodes - 2
                
                combinedVolume = gridVolume(i) + gridVolume(j);
                
                if (gridVolume(k) <= combinedVolume && combinedVolume <= gridVolume(k+1))
                    
                    splittingOperators(i,j,k) = (gridVolume(k+1) - combinedVolume )/...
                        (gridVolume(k+1) - gridVolume(k));
                    
                else if (gridVolume(k-1) <= combinedVolume && combinedVolume <= gridVolume(k))
                
                    splittingOperators(i,j,k) = ( combinedVolume - gridVolume(k-1))/...
                        (gridVolume(k) - gridVolume(k-1));
                    
                else
                    
                    splittingOperators(i,j,k) = 0.0;
                    
                end
                
                
            end
       end
   end

end

