function B = binMatrixSum(A, binSize, dim)
    if nargin < 3
        dim = 1; % Standard: entlang erster Dimension
    end

    sz = size(A);
    nBins = floor(sz(dim) / binSize);

    szB = sz;
    szB(dim) = nBins;
    B = zeros(szB);

    idx = reshape(1:binSize*nBins, binSize, nBins);

    if dim == 1
        for k = 1:nBins
            B(k,:) = sum(A(idx(:,k),:),1);
        end
    elseif dim == 2
        for k = 1:nBins
            B(:,k) = sum(A(:,idx(:,k)),2);
        end
    else
        error('Nur Dimension 1 oder 2 werden unterstützt.');
    end
end
