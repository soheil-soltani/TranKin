%This function implements the Newton root-finding routine to calculate the
%   regularization parameter.
%   The user is first asked to determine the under-relaxation factor for
%   the iterative root-finding algorithm. For a very noisy signal, more 
%   under-relaxation may be required; therefore, a value of about 0.2 can
%   be recommended.
%   A starting value for the regularization parameter is then calculated
%   such that a positive value for discrepancy is obtained.
%   Using this starting value, Newton's method is initialized to calculate
%   the final value of the regularization parameter for either the special
%   case of tau = 0 or a general case of tau ~= 0.

% Copyright © 2014 Soheil Soltani

% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files 
% (the "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function [alpha,A,B,index] = frootf(alpha_0, G, Sigma, T, y, t, delta,...
    tau, K)

%1-Initialization
UR_input = ['Please enter the under-relaxation factor '...
    '(from the interval (0, 1]) for the iterative root-finding '...
    'routine... '];
omega = input(UR_input);
if (omega > 1 || omega <= 0)
    error('Incorrect input. Please restart the program.')
end
tic

%2- Starting value of the Newton root-finding method
discrepancy = -1;       %Dummy value to initialize the while-loop
while discrepancy < 0
    alpha = alpha_0;
    solution = T*((Sigma./(Sigma.^2 + alpha)).*(G'*y));
    if tau == 0
        discrepancy = norm(K*solution - y).^2 - delta^2;
    else
        diff_x = zeros(size(t));
        diff_x(2:end,1) = diff(solution)./diff(t);
        x_norm = sqrt(trapz(t, (solution.^2 + diff_x.^2)));
        discrepancy = norm(K*solution - y).^2 - (delta + tau*x_norm)^2;
    end
    if discrepancy > 0
        fprintf(2,'Correct starting value was computed.');
        break   %This way, the already calculated initial guess 
                    %is preserved.
    else
        alpha_0 = alpha_0*2
    end
end

%3a- regularization parameter for the special case of tau = 0
if tau == 0             %Exact kernel    
    alpha = alpha_0;
    threshold = eps;   
    j = 1;              %Counter
    discrepancy = 1;    %Dummy value to initialize the while-loop
    while abs(discrepancy) > threshold
        A(j) = alpha;
        solution = T*((Sigma./(Sigma.^2 + alpha)).*(G'*y));
        discrepancy = norm(K*solution - y).^2 - delta^2;
        B(j) = discrepancy;        
        derivative = alpha*(((1./(Sigma.^2 + alpha)).*solution)'*solution); 
        if (discrepancy <= threshold)
            break
        end
        alpha = alpha - omega*(discrepancy/derivative);
        j = j + 1;     
    end
index = 1:1:j;

exec_time = toc;

output_msg = ['Convergence criterion was reached after %g iteration(s)'...
    '.\nCalculations took %g seconds.\nThe initial guess for the'...
    ' regularization parameter was: %g.\nThe calculated regularization'...
    ' parameter is: %g.\n '];
msgbox(sprintf(output_msg,j,exec_time,alpha_0,alpha),...
    'Calculations Completed');

%3b- regularization parameter for a general case (tau ~= 0)
else                    %Approximate kernel, i.e. tau ~= 0
    alpha = alpha_0;
    threshold = eps;   
    j = 1;              %Counter
    discrepancy = 1;    %Dummy value to initialize the while-loop
    while abs(discrepancy) > threshold
        A(j) = alpha;
        solution = T*((Sigma./(Sigma.^2 + alpha)).*(G'*y));
        diff_x = zeros(size(t));
        diff_x(2:end,1) = diff(solution)./diff(t);
        x_norm = sqrt(trapz(t, (solution.^2 + diff_x.^2)));
        discrepancy = norm(K*solution - y).^2 - (delta + tau*x_norm)^2;
        B(j) = discrepancy;
        derivative = (alpha + tau*delta/x_norm + tau^2)*...
            (((1./(Sigma.^2 + alpha)).*solution)'*solution);
        if (discrepancy <= threshold)
            break
        end
        alpha = alpha - omega*(discrepancy/derivative);      
        j = j + 1;
    end
index = 1:1:j;   

exec_time = toc;

output_msg = ['Convergence criterion was reached after %g iteration(s)'...
    '.\nCalculations took %g seconds.\nThe initial guess for the'...
    ' regularization parameter was: %g.\nThe calculated regularization'...
    ' parameter is: %g.\n '];
msgbox(sprintf(output_msg,j,exec_time,alpha_0,alpha),...
    'Calculations Completed');
end