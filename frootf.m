<<<<<<< HEAD
%   This function implements two root-finding schemes to calculate the
%   regularization parameter; the Newton method and the bisection method.
%   A starting value for the regularization parameter is first calculated
=======
%This function implements the Newton root-finding routine to calculate the
%   regularization parameter.
%   The user is first asked to determine the under-relaxation factor for
%   the iterative root-finding algorithm. For a very noisy signal, more 
%   under-relaxation may be required; therefore, a value of about 0.2 can
%   be recommended.
%   A starting value for the regularization parameter is then calculated
>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab
%   such that a positive value for discrepancy is obtained.
%   Using this starting value, Newton's method is initialized to calculate
%   the final value of the regularization parameter for either the special
%   case of tau = 0 or a general case of tau ~= 0.
<<<<<<< HEAD
%   If the Newton method was chosen, the user is first asked to determine
%   the under-relaxation factor for the iterative root-finding algorithm.
%   For a very noisy signal, more under-relaxation may be required; 
%   therefore, a value of about 0.2 can be recommended.
%   If the bisection method was chosen, the algortihm is expected to be
%   much more stable. The routine controls that the starting interval for
%   the iterative bisection method is chosen such that the sign of the
%   discrepancy functional is different at the two lower- and upper bounds.
=======
>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab

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
<<<<<<< HEAD
    tau, K, interval)

%1-Initialization
threshold = eps;   
j = 1;              %Counter

%2- Specifying alpha_0 such that discrepancy > 0
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
            fprintf('\nCorrect starting value was computed ');
            fprintf('(discrepancy > 0).\n ');
            break   %This way, the already calculated initial guess 
                        %is preserved.
        else
            alpha_0 = alpha_0*2;
        end
    end
 if (nargin == 9)                            %i.e. Newton method was chosen
    Var_Exist_Chk = evalin('base','exist(''omega_from_GUI'', ''var'');');
    if Var_Exist_Chk == 0    
        UR_input = ['\nPlease enter the under-relaxation factor '...
            '(from the interval (0, 1]) for the \niterative root-finding '...
            'routine... '];
        omega = input(UR_input);
        if (omega > 1 || omega <= 0)
            error('\nIncorrect input. Please restart the program.')
        end
    else
        omega = evalin('base', 'omega_from_GUI');
        if  isnan(omega)
            UR_input = ['\nPlease enter the under-relaxation factor '...
            '(from the interval (0, 1]) for the \niterative root-finding '...
            'routine... '];
            omega = input(UR_input);
        end
            
        if (omega > 1 || omega <= 0)
            UR_errormsg = ['\nIncorrect under-relaxation input. Please '...
                'restart the program.'];
            error(UR_errormsg)
        end
    end
    tic

%3-a regularization parameter for the special case of tau = 0
    if tau == 0             %Exact kernel    
        alpha = alpha_0;
        discrepancy = 1;    %Dummy value to initialize the while-loop
        while abs(discrepancy) > threshold
            A(j) = alpha;
            solution = T*((Sigma./(Sigma.^2 + alpha)).*(G'*y));
            discrepancy = norm(K*solution - y).^2 - delta^2;        
            B(j) = discrepancy;        
            derivative = ...
                alpha*(((1./(Sigma.^2 + alpha)).*solution)'*solution); 
            if (discrepancy <= threshold)
                break
            end
            alpha = alpha - omega*(discrepancy/derivative);
            j = j + 1;     
        end
    index = 1:1:j;
    exec_time = toc;
    output_msg = ['\nConvergence criterion was reached after'...
        ' %g iteration(s).\nCalculations took %g seconds.\nThe'...
        ' initial guess for the regularization parameter was:'...
        ' %g.\nThe calculated regularization parameter is: %g.\n '];
    msgbox(sprintf(output_msg,j,exec_time,alpha_0,alpha),...
        'Calculations Completed');

%3-b regularization parameter for a general case (tau ~= 0)
    else                     %Approximate kernel, i.e. tau ~= 0
        alpha = alpha_0;        
        discrepancy = 1;     %Dummy value to initialize the while-loop
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
    output_msg = ['\nConvergence criterion was reached after'...
        ' %g iteration(s).\nCalculations took %g seconds.\nThe'...
        ' initial guess for the regularization parameter was:'...
        ' %g.\nThe calculated regularization parameter is: %g.\n '];
    msgbox(sprintf(output_msg,j,exec_time,alpha_0,alpha),...
        'Calculations Completed');
    end
else                                     %i.e. Bisection method was chosen
%4-a    
    tic
    if (size(interval,1) == 1 &&  size(interval,2) == 2)
        a = interval(1); b = interval(2);
        
% The algorithm now controls to ensure that the sign of the discrepancy 
% functional changes over the specified interval
% -------------------------------------------------------------------------
%4-b
        if tau == 0         %Exact kernel
            solution_a = T*((Sigma./(Sigma.^2 + a)).*(G'*y));
            solution_b = T*((Sigma./(Sigma.^2 + b)).*(G'*y));
            f_a = norm(K*solution_a - y).^2 - delta^2;
            f_b = norm(K*solution_b - y).^2 - delta^2;
            if f_a*f_b > 0
                error_bisec = ['\nIncorrect interval. f(a)f(b) > 0.'...
                    ' Please restart the program.'];
                error(error_bisec)
            end
%4-c            
        else                % Approximate kernel, i.e. tau ~= 0
            solution_a = T*((Sigma./(Sigma.^2 + a)).*(G'*y));
            solution_b = T*((Sigma./(Sigma.^2 + b)).*(G'*y));
            
            diff_x_a = zeros(size(t));
            diff_x_a(2:end,1) = diff(solution_a)./diff(t);
            
            diff_x_b = zeros(size(t));
            diff_x_b(2:end,1) = diff(solution_b)./diff(t);
            
            a_norm = sqrt(trapz(t, (solution_a.^2 + diff_x_a.^2)));
            b_norm = sqrt(trapz(t, (solution_b.^2 + diff_x_b.^2)));
            
            f_a = norm(K*solution_a - y).^2 - (delta + tau*a_norm)^2;
            f_b = norm(K*solution_b - y).^2 - (delta + tau*b_norm)^2;
            if f_a*f_b > 0
                error_bisec = ['\nIncorrect interval. f(a)f(b) > 0.'...
                    ' Please restart the program.'];
                error(error_bisec)
            end
        end
%--------------------------------------------------------------------------         
    else                    % The interval is calculated by the algorithm
%4-d        
        
% The algorithm starts from the calcualted initial guess for which it was
% already ensured in part 2 of the current routine that a positive value 
% for discrepancy is calculated. Next, the algorithm keeps halving the 
% intial guess until a negative value for discrepancy is encountered. This 
% value is then going to be the lower bound of the interval. The upper
% bound is obtained by multiplying the lower bound by 2. Thereby, the root
% is trapped in a close neighbourhood.

        if tau == 0         %Exact kernel
            err = 1;        %Dummy value to initialize the while-loop
            a = alpha_0;
            while err > 0
                solution = T*((Sigma./(Sigma.^2 + a)).*(G'*y));
                err = norm(K*solution - y).^2 - delta^2;
                if err < 0
                    b = a*2;
                    break
                end
                a = a/2;
            end
% The algorithm now controls to ensure that the sign of the discrepancy 
% functional changes over the specified interval
% -------------------------------------------------------------------------
            solution_a = T*((Sigma./(Sigma.^2 + a)).*(G'*y));
            solution_b = T*((Sigma./(Sigma.^2 + b)).*(G'*y));
            f_a = norm(K*solution_a - y).^2 - delta^2;
            f_b = norm(K*solution_b - y).^2 - delta^2;
            if f_a*f_b > 0
                error('\nIncorrect interval. Software failure.')
            end 
% -------------------------------------------------------------------------

        else                %Approximate kernel, i.e. if tau ~= 0
            err = 1;
            a = alpha_0;
            while err > 0
                solution = T*((Sigma./(Sigma.^2 + a)).*(G'*y));
                diff_solution = zeros(size(t));
                diff_solution(2:end,1) = diff(solution)./diff(t);
                solution_norm = ...
                    sqrt(trapz(t, (solution.^2 + diff_solution.^2)));
                err = norm(K*solution - y).^2 - ...
                    (delta + tau*solution_norm)^2;                
                if err < 0
                    b = a*2;
                    break
                end
                a = a/2;
            end
% The algorithm now controls to ensure that the sign of the discrepancy 
% functional changes over the specified interval
% -------------------------------------------------------------------------            
            solution_a = T*((Sigma./(Sigma.^2 + a)).*(G'*y));
            solution_b = T*((Sigma./(Sigma.^2 + b)).*(G'*y));
            
            diff_x_a = zeros(size(t));
            diff_x_a(2:end,1) = diff(solution_a)./diff(t);
            
            diff_x_b = zeros(size(t));
            diff_x_b(2:end,1) = diff(solution_b)./diff(t);
            
            a_norm = sqrt(trapz(t, (solution_a.^2 + diff_x_a.^2)));
            b_norm = sqrt(trapz(t, (solution_b.^2 + diff_x_b.^2)));
            
            f_a = norm(K*solution_a - y).^2 - (delta + tau*a_norm)^2;
            f_b = norm(K*solution_b - y).^2 - (delta + tau*b_norm)^2;
            if f_a*f_b > 0
                error('\nIncorrect interval. Software failure.')
            end         
% -------------------------------------------------------------------------                        
        end
    end
    
%4-e Now the interval for starting the bisection method is known.    
            
    p = (a + b)/2;
    alpha = alpha_0;
    discrepancy = 1;    %Dummy value to initialize the while-loop
    if tau == 0         %Exact kernel
        while abs(discrepancy) > threshold        
            A(j) = alpha;
            solution_a = T*((Sigma./(Sigma.^2 + a)).*(G'*y));
            solution_b = T*((Sigma./(Sigma.^2 + b)).*(G'*y));
            solution_p = T*((Sigma./(Sigma.^2 + p)).*(G'*y));

            f_a = norm(K*solution_a - y).^2 - delta^2;
            f_b = norm(K*solution_b - y).^2 - delta^2;
            f_p = norm(K*solution_p - y).^2 - delta^2;
            
            if f_a*f_p < 0
                b = p;
            else
                a = p;
            end
            p = (a + b)/2;
            alpha = p;
            solution_p = T*((Sigma./(Sigma.^2 + p)).*(G'*y));
            discrepancy = norm(K*solution_p - y).^2 - delta^2;
            B(j) = discrepancy;
            if abs(discrepancy) <= threshold
                break
            end
            j = j + 1;
        end
        index = 1:1:j;
        exec_time = toc;
        output_msg = ['\nConvergence criterion was reached after'...            
        '%g iteration(s). \nCalculations took %g seconds.\nThe initial '...
        ' guess for the regularization parameter was: %g.\nThe '...
        ' calculated regularization parameter is: %g.\n '];
        msgbox(sprintf(output_msg,j,exec_time,alpha_0,alpha),...
            'Calculations Completed');
    else
        while abs(discrepancy) > threshold        
            A(j) = alpha;
            solution_a = T*((Sigma./(Sigma.^2 + a)).*(G'*y));
            solution_b = T*((Sigma./(Sigma.^2 + b)).*(G'*y));
            solution_p = T*((Sigma./(Sigma.^2 + p)).*(G'*y));
            
            diff_x_a = zeros(size(t));
            diff_x_a(2:end,1) = diff(solution_a)./diff(t);
            
            diff_x_b = zeros(size(t));
            diff_x_b(2:end,1) = diff(solution_b)./diff(t);
            
            diff_x_p = zeros(size(t));
            diff_x_p(2:end,1) = diff(solution_p)./diff(t);
            
            a_norm = sqrt(trapz(t, (solution_a.^2 + diff_x_a.^2)));
            b_norm = sqrt(trapz(t, (solution_b.^2 + diff_x_b.^2)));
            p_norm = sqrt(trapz(t, (solution_p.^2 + diff_x_p.^2)));
            
            f_a = norm(K*solution_a - y).^2 - (delta + tau*a_norm)^2;
            f_b = norm(K*solution_b - y).^2 - (delta + tau*b_norm)^2;           
            f_p = norm(K*solution_p - y).^2 - (delta + tau*p_norm)^2;
            
            if f_a*f_p < 0
                b = p;
            else
                a = p;
            end
            p = (a + b)/2;
            alpha = p;
            solution_p = T*((Sigma./(Sigma.^2 + p)).*(G'*y));
            discrepancy = norm(K*solution_p - y).^2 -...
                (delta + tau*p_norm)^2;
            B(j) = discrepancy;
            if abs(discrepancy) <= threshold
                break
            end
            j = j + 1;
        end
        index = 1:1:j;
        exec_time = toc;
        output_msg = ['\nConvergence criterion was reached after '...
        '%g iteration(s). \nCalculations took %g seconds.\nThe initial '...
        ' guess for the regularization parameter was: %g.\nThe '...
        ' calculated regularization parameter is: %g.\n '];
        msgbox(sprintf(output_msg,j,exec_time,alpha_0,alpha),...
            'Calculations Completed');
    end
 end
=======
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
>>>>>>> 37411c16ab68b52ea1b64f09192ec4d8d43284ab
