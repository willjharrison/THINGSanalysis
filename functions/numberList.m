function output = numberList(N,repeats)

output = reshape(repmat(1:N,repeats,1),repeats*N,1);