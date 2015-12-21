function [ alpha_0 ] = initial_guess(K, y, delta, tau)
%Subroutine to calculate an initial guess for initialization of the Newton
%root-finding algorithm.
%   A random number between [1, norm(y)/delta] is generated which is used
%   in the subsequent lemma to calculate the starting point.

a = 1;

if delta == 0
    b = norm(y)/(delta + eps); 
    %To prevent numerical error due to division by zero
else
    b = norm(y)/delta;
end

c = a + (b - a)*rand(1);   
alpha_tilde = (norm(K)^2*c*delta)/(norm(y) - c*delta);
alpha_0 = norm(K)*(tau + (sqrt(tau^2 + (alpha_tilde + tau^2)*...
    (c^2 - 1)) + tau)/(c^2 - 1));
end