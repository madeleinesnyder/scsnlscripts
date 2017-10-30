function outmatrix = tril2full(invector, nrows)

outmatrix = NaN(nrows, nrows);
index = 0;
for j = 1:nrows
  for i = (j+1):nrows
    index = index + 1;
    outmatrix(i,j) = invector(index);
    outmatrix(j,i) = invector(index);
  end
end

end