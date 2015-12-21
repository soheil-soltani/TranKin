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
<<<<<<< HEAD
F_input_msg = ['\nPlease enter the name and extension (e.g. data.txt)\n'...
=======
F_input_msg = ['Please enter the name and extension (e.g. data.txt)\n'...
>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab
    'of the .txt file containing the variables t, y, g, respectively: '];
file_name = input(F_input_msg, 's');
fileID = fopen(file_name,'r');
data = fscanf(fileID, '%g', [3 Inf]);
data = data';

t = data(:,1);      %Time series at which data has been acquired
y = data(:,2);      %Specifies the right-hand side
g = data(:,3);      %Specifies the kernel

l = length(t);
<<<<<<< HEAD
% 
=======

>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab
%2-Singular Value Decomposition
K = tril(toeplitz(g));    %Reforming the kernel into a convolution operator
[G,Sigma,T] = svd(K);
Sigma = diag(Sigma);

%3-Error
%delta and tau are experimental uncertainties in specifying the 
%   right-hand side and the kernel, respectively.
<<<<<<< HEAD
delta_msg = ['\nPlease enter the error in specifying the righ-hand '...
    'side. delta = '];
delta = input(delta_msg);
if (delta < 0)
    errormsg = ['A positive value for error must be entered. Please'...
        ' restart the program.'];
    error(errormsg)
end

tau_msg = ['\nPlease enter the error in specifying the kernel. tau = '];
tau = input(tau_msg);
if (tau < 0)
    errormsg = ['A positive value for error must be entered. Please'...
        ' restart the program.'];
    error(errormsg)
end

%4-Regularization parameter
alpha_input_msg = ['\nPlease enter the regularization parameter (from '...
    ' the interval (0, Inf) ) or \npress <Enter> to let the software'...
    ' calculate it: '];
user_inp_alpha = input(alpha_input_msg);
if (isempty(user_inp_alpha) == false)
    %There is no iterative procedure in this case.
	if (user_inp_alpha > 0)
        fprintf('\nalpha has been entered successfully.')
        alpha = user_inp_alpha;
=======
delta_msg = ['Please enter the error in specifying the righ-hand side.'...
    ' delta = '];
delta = input(delta_msg);

tau_msg = ['Please enter the error in specifying the kernel. tau = '];
tau = input(tau_msg);

%4-Regularization parameter
alpha_input_msg = ['Please enter the regularization parameter (from the '...
    'interval (0, Inf) ) or press <Enter> to let the '...
    'software \ncalculate it: '];
getkey = input(alpha_input_msg);
if getkey ~= 13			%There is no iterative procedure in this case.
	if (getkey > 0)
        fprintf('alpha has been entered successfully.')
        alpha = getkey;
>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab
        plot_flag = 0;	%This prevents plotting the convergence history, 
                        %which is relevant to the iterative routine only.
	else
		errormsg = ['Incorrect input. Note that alpha > 0. '...
		'Please restart the program.'];
        error(errormsg)
	end
else
    plot_flag = 1; 		%This ensures plotting the convergence history of 
                        %the following iterative routine.
<<<<<<< HEAD
                        
%4-a Choosing the iterative root-finding scheme: Newton vs. Bisection                        
    method = input(['\nPlease choose the root-finding scheme: '...
        '(Newton or Bisection)[N/B]...'],'s');
%4-b Initial guess for the Newton method

    if (method == 'N')          %Newton method has been chosen   
        inputmsg = ['\nPlease enter the initial guess for the Newton '...
            'root-finding routine \n(from the interval (0, Inf) ) or '...
            'press <Enter> \nto let the software calculate it: '];   
        user_initial = input(inputmsg);
        if (isempty(user_initial) == false)
                                %User has entered the initial guess
            if (user_initial > 0)
                fprintf('\nThe initial guess has been entered correctly.')
                alpha_0 = user_initial;
                [alpha,A,B,index] = frootf(alpha_0, G, Sigma, T, y, t,...
                delta, tau, K);
            else
                errormsg = ['Incorrect input. Note that alpha > 0. '...
                    'Please restart the program.'];
                error(errormsg)
            end
        else                    %The software calculated the initial guess
            alpha_0 = initial_guess(K, y, delta, tau);
            [alpha,A,B,index] = frootf(alpha_0, G, Sigma, T, y, t,...
                delta, tau, K);
        end
    elseif (method == 'B')      %Bisection method has been chosen
        
%4-c The initial interval for starting the bisection method

        alpha_0 = initial_guess(K, y, delta, tau);
        inputmsg = ['\nPlease enter alpha_1 and alpha_2 enclosed in '...
            'brackets (without comma) \nto set the starting interval '...
            '[alpha_1 alpha_2] for initializing the bisection method '...
            '\nor press <Enter> to let the software calculate it: '];
        interval = input(inputmsg);
        if ((size(interval,1) == 1 &&  size(interval,2) == 2))
            [alpha,A,B,index] = frootf(alpha_0, G, Sigma, T, y, t,...
                delta, tau, K, interval);                 
        elseif (isempty(interval) == true)
            [alpha,A,B,index] = frootf(alpha_0, G, Sigma, T, y, t,...
                delta, tau, K, interval);
        else
            error('Incorrect input. Please restart the program.');
        end
    else
        error('Incorrect choice. Please restart the program.');
    end
=======
	inputmsg = ['Please enter the initial guess for the Newton '...
        'root-finding routine (from the interval (0, Inf) ) or press '...
        '<Enter> to \nlet the software calculate it: '];
	getkey = input(inputmsg);
	if getkey ~= 13
		if (getkey > 0)
			fprintf('The initial guess has been entered successfully.')
        	alpha_0 = getkey;
		else
			errormsg = ['Incorrect input. Note that alpha > 0. '...
                'Please restart the program.'];
			error(errormsg)
		end
    else
		alpha_0 = initial_guess(K, y, delta, tau);
	end
    [alpha,A,B,index] = frootf(alpha_0, G, Sigma, T, y, t, delta, tau, K);
>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab
end

%5-Computing the regularized solution
solution = T*((Sigma./(Sigma.^2 + alpha)).*(G'*y));

%6-Post-processing
<<<<<<< HEAD
fprintf('\nPlot the results?[Y/N]... ');
=======
fprintf(2,'\nPlot the results?[Y/N]... ');
>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab
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
<<<<<<< HEAD
        if (method == 'N')
            title('Convergence monitor (Newton method)','FontSize',12)
        else
            title('Convergence monitor (Bisection method)','FontSize',12)
        end
=======
        title('Convergence monitor','FontSize',12)
        xlabel('Iteration','FontSize',12)
>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab
        ylabel('\alpha','FontSize',12)
        subplot(2,1,2)
        plot(index,B,'r')
        xlabel('Iteration','FontSize',12)
<<<<<<< HEAD
        ylabel('\rho_\eta(\alpha)','FontSize',12)        
=======
        ylabel('\rho_\eta(\alpha)','FontSize',12)
>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab
    end
end