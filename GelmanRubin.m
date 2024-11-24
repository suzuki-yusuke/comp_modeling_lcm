function R_hat = GelmanRubin(L,m,Z)

num_parameters = size(Z,2);
k = ceil(L/m);
x = mat2cell(Z,ones(1,m).*k,ones(1,num_parameters));

phi = cellfun(@mean,x);
R_hat = zeros(1,num_parameters);
for h = 1:num_parameters
    B = zeros(1,m);
    s = zeros(1,m);
    for j = 1:m
        B(j) = (phi(j,h)-mean(phi(:,h)))^2;
        s(j) = (1/(k-1))*sum((x{j,h}-phi(j,h)).^2);
    end
    B = (k/(m-1))*sum(B);
    W = mean(s);
    R_hat(h) = sqrt(((k-1)/k)+(1/k)*(B/W));
end



