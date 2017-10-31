function x = GP_rnd(GP)
x = mvnrnd(GP.mu, GP.S,GP.M)';
end