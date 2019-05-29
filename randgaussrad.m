function randrad = randgaussrad(mu,sigma)

% Generate normally distributed random numbers
randrad = randn(1);

% Shift to match given mean and std
randrad = mu + randrad * sigma;

end
