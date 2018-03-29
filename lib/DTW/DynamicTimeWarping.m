function [p,q,D] = DynamicTimeWarping(A, B)%picos_source, picos_target
% [p,q] = DynamicTimeWarping(A, B)
% A: d*N matrix, d is the dimensionality of a feature vector
% B: d*M matrix
%    Use dynamic programming to find a min-cost path through matrix M.
%    Return state sequence in p,q

% EA = sqrt(sum(A.^2));
% EB = sqrt(sum(B.^2));

% this is 10x faster
% M = (A'*B)./(EA'*EB);
% M = A' * B;


M = pdist2(A', B');
% M = M.^2;

% [rows,R]=size(A);
% [rows,C]=size(B);
% for r=1:R
%     for c=1:C
%         M(r,c) = sum((A(:,r)-B(:,c)).^2);
%     end
% end
%AAA = repmat(A(:), 1, c);
%BBB = repmat(B(:)', r, 1);

%M=(AAA - BBB).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[r,c] = size(M);

% costs
D = zeros(r+1, c+1);
D(1,:) = NaN;
D(:,1) = NaN;
D(1,1) = 0;
D(2:(r+1), 2:(c+1)) = M;

% traceback
phi = zeros(r,c);

for i = 1:r; 
  for j = 1:c;
    [dmin, tb] = min([D(i, j), D(i, j+1), D(i+1, j)]);
    D(i+1,j+1) = D(i+1,j+1)+dmin;
    phi(i,j) = tb;
  end
end

% Traceback from top left
i = r; 
j = c;
p = i;
q = j;
while i > 1 & j > 1
  tb = phi(i,j);
  if (tb == 1)
    i = i-1;
    j = j-1;
  elseif (tb == 2)
    i = i-1;
  elseif (tb == 3)
    j = j-1;
  else    
    error;
  end
  p = [i,p];
  q = [j,q];
end

% Strip off the edges of the D matrix before returning
D = D(2:(r+1),2:(c+1));
