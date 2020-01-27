%% Interpolate model domain to match domain size
function [domain2D] = interpDomain(bath,config)

%     bath = bath.bath;
    domain2D = ones(config.Nx,config.Ny);
    r = config.Nx / size(bath,1);
    for ii = 1:size(bath,1)
        for jj = 1:size(bath,2)
            gridcell = ones(r,r);
%             if ii > 1 && ii < size(bath,1) && jj > 1 && jj < size(bath,2)
%                 if isnan(bath(ii-1,jj)) && isnan(bath(ii,jj-1)) && ...
%                         ~isnan(bath(ii+1,jj)) && ~isnan(bath(ii,jj+1))
%                     for col = 1:r
%                         gridcell(1:r-col+1,col) = NaN;
%                     end
%                 elseif ~isnan(bath(ii-1,jj)) && ~isnan(bath(ii,jj-1)) && ...
%                         isnan(bath(ii+1,jj)) && isnan(bath(ii,jj+1))
%                     for col = 1:r
%                         gridcell(r-col+1:end,col) = NaN;
%                     end
%                 elseif isnan(bath(ii-1,jj)) && ~isnan(bath(ii,jj-1)) && ...
%                         ~isnan(bath(ii+1,jj)) && isnan(bath(ii,jj+1))
%                     for col = 1:r
%                         gridcell(1:col,col) = NaN;
%                     end
%                 elseif ~isnan(bath(ii-1,jj)) && isnan(bath(ii,jj-1)) && ...
%                         isnan(bath(ii+1,jj)) && ~isnan(bath(ii,jj+1))
%                     for col = 1:r
%                         gridcell(col:end,col) = NaN;
%                     end
%                 else
%                     gridcell(:,:) = bath(ii,jj);
%                 end
%             else
%                 gridcell(:,:) = bath(ii,jj);
%             end
            gridcell(:,:) = bath(ii,jj);
            domain2D((ii-1)*r+1:ii*r,(jj-1)*r+1:jj*r) = gridcell;
        end
    end



end