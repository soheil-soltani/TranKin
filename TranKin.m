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


clear
pack
clc

%1-Initialization
F_input_msg = ['Please enter the name and extension (e.g. data.txt)\n'...
    'of the .txt file containing the variables t, y, g, respectively: '];
file_name = input(F_input_msg, 's');
fileID = fopen(file_name,'r');
data = fscanf(fileID, '%g', [3 Inf]);
data = data';

t = data(:,1);      %Time series at which data has been acquired
y = data(:,2);      %Specifies the right-hand side
g = data(:,3);      %Specifies the kernel

l = length(t);

%2-Singular Value Decomposition
K = tril(toeplitz(g));    %Reforming the kernel into a convolution operator
[G,Sigma,T] = svd(K);
Sigma = diag(Sigma);

%3-Error
%delta and tau are experimental uncertainties in specifying the 
%   right-hand side and the kernel, respectively.
delta_msg = ['Please enter the error in specifying the righ-hand side.'...
    ' delta = '];
delta = input(delta_msg);

tau_msg = ['Please enter the error in specifying the kernel. tau = '];
tau = input(tau_msg);

%4-Regularization parameter
alpha_input_msg = ['Please enter the regularization parameter (from the '...
    'interval [0, 1]) or press <Enter> to let the '...
    'software \ncalculate it: '];
getkey = input(alpha_input_msg);
if getkey ~= 13
	if (getkey<=1 && getkey >=0)
        fprintf('alpha has been entered successfully.')
        alpha = getkey;
        plot_flag = 0;
	else
		errormsg = ['Incorrect input. Note that 0 <= alpha <= 1. '...
		'Please restart the program.'];
        error(errormsg)
	end
else
    plot_flag = 1;
	inputmsg = ['Please enter the initial guess for the Newton '...
        'root-finding routine (from the interval [0, 1]) or press '...
        '<Enter> to \nlet the software calculate it: '];
	getkey = input(inputmsg);
	if getkey ~= 13
		if (getkey<=1 && getkey >=0)
			fprintf('The initial guess has been entered successfully.')
        	alpha_0 = getkey;
		else
			errormsg = ['Incorrect input. Note that 0 <= alpha <= 1. '...
                'Please restart the program.'];
			error(errormsg)
		end
	else
		alpha_0 = initial_guess(K, y, delta, tau);
	end
    [alpha,A,B,index] = frootf(alpha_0, G, Sigma, T, y, t, delta, tau, K);
end

%5-Computing the regularized solution
solution = T*((Sigma./(Sigma.^2 + alpha)).*(G'*y));

%6-Post-processing
fprintf(2,'\nPlot the results?[Y/N]... ');
plot_confirm = input('','s');
if plot_confirm == 'Y'
    figure(1)
    plot(t,solution)
    title('Regularized solution', 'FontSize', 12)
    xlabel('Time, [arb. u.]','FontSize',12)
    ylabel('[arb. u.]','FontSize',12)
    
    if plot_flag == 1
        figure(2)
        subplot(2,1,1)
        plot(index,A)
        title('Convergence monitor','FontSize',12)
        xlabel('Iteration','FontSize',12)
        ylabel('\alpha','FontSize',12)
        subplot(2,1,2)
        plot(index,B,'r')
        xlabel('Iteration','FontSize',12)
        ylabel('\rho_\eta(\alpha)','FontSize',12)
    end
end