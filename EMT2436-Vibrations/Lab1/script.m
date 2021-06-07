% Aurthor:          Jilks Smith
% Field of Study:   Bsc. Mechatronic Enhineering
% Course:           Vibrations
% Year of Study:    4th


% Beginning of script

disp("Linear Constant Coefficient ODE of the form:   a y'' + b y' + c y = 0 ");



roots_type = 1;
% 1 for real and distinct roots
% 2 for real and repeated roots
% 3 for complex roots

% Prompt user input, a,b,c and intial conditions

a = double(input("Enter the value of a: "));
b = double(input("Enter the value of b: "));
c = double(input("Enter the value of c: "));

first_initial_condition = double(input("Enter the initial condition y(0): "));
second_initial_condition = double(input("Enter the initial condition y'(0): "));

roots_type = determineTypeOfRoots(a, b, c);

if(roots_type == 1)
    sol = solveEquationWithRealDistinctRoots(a,b,c,first_initial_condition,second_initial_condition);
    disp("Final solution for real and distinct roots problem: ");
    disp(sol);
end

if(roots_type == 2)
    sol = solveEquationWithRepeatedRealRoots(a,b,c,first_initial_condition,second_initial_condition);
    disp("Final solution for real and repeated roots problem: ");
    disp(sol);
end

if(roots_type == 3)
    sol = solveEquationWithComplexRoots(a,b,c,first_initial_condition,second_initial_condition);
    disp("Final solution for complex roots problem: ");
    disp(sol);
end

dSolveSolution = solveEquationWithDSolve(a,b,c, first_initial_condition, second_initial_condition);
disp("DSolve Solution");
disp(dSolveSolution);




% Funtion that returns integer for check for nature of roots b^2 - 4ac
function root_type_int = determineTypeOfRoots(a, b, c)
if(b*b > 4*a*c)
    root_type_int = 1; % Real and distinct roots
end

if(b*b == 4*a*c)
    root_type_int = 2; % Real and repeated roots
end

if(b*b < 4*a*c)
    root_type_int = 3; % Complex roots
end
    
end 




% Function to solve for solution of real and distinct roots
function problemSolution = solveEquationWithRealDistinctRoots(a,b,c,y1,y2)
    syms c1 c2 x

    first_root = ( -b + sqrt(b*b - 4*a*c) )/ (2*a);
    second_root = ( -b - sqrt(b*b - 4*a*c) )/ (2*a);

    generalSolution = c1 * exp(first_root*x) +  c2 * exp(second_root*x);

    eqn1 = c1 * exp(first_root * 0) +  c2 * exp(second_root * 0) == y1; % = y1 when x=0

    diff_eqn1 = first_root * c1 * exp(first_root * 0) +  second_root * c2 * exp(second_root * 0) == y2; % = y2 when x=0

    sol = solve([eqn1, diff_eqn1], [c1, c2]);
    c1_solution = sol.c1;
    c2_solution = sol.c2;


    %The final solution
    problemSolution = c1_solution * exp(first_root * x) +  c2_solution * exp(second_root * x);
    


end



% Function to solve for solution of real and repeated roots 
function problemSolution = solveEquationWithRepeatedRealRoots(a,b,c,y1,y2)
    syms c1 c2 x

    % first_root = second_root = root
    root = ( -b + sqrt(b*b - 4*a*c) )/ (2*a);

    generalSolution = ((c1 * x) + c2) * exp(root * x);

    eqn1 = ((c1 * 0) + c2) * exp(root * 0) == y1; % = y1,  x=0
    diff_eqn1 = c1 * exp(root * 0) + root * (c1 * 0 + c2) * exp(root * 0) == y2; % = y2,  x=0

    sol = solve([eqn1, diff_eqn1], [c1, c2]);
    c1_solution = sol.c1;
    c2_solution = sol.c2;

    %The final solution
    problemSolution = ((c1_solution * x) + c2_solution) * exp(root * x);


end

% Function to solve for solution of complex roots roots

function problemSolution = solveEquationWithComplexRoots(a,b,c,y1,y2)
    syms c1 c2 x

    first_root = ( -b + sqrt(b*b - 4*a*c) )/ (2*a);
    second_root = ( -b - sqrt(b*b - 4*a*c) )/ (2*a);

    alpha_part = real(first_root);
    beta_part= imag(first_root);

    generalSolution = exp(alpha_part * x) * (c1 * cos(beta_part * x) + c2 * sin(beta_part * x));

    eqn1 = exp(alpha_part * 0) * (c1 * cos(beta_part * 0) + c2 * sin(beta_part * 0)) == y1; % = y1 when x=0
    diff_eqn1 = ( (alpha_part * exp(alpha_part * 0) ) * (c1 * cos(beta_part * 0) + c2 * sin(beta_part * 0)) ) + ( exp(alpha_part * 0) * ( c1 * beta_part * -sin(beta_part * 0) + c2 * beta_part * cos(beta_part * 0)) ) == y2 ; % = y2 when x=0
   

    sol = solve([eqn1, diff_eqn1], [c1, c2]);
    c1_solution = sol.c1;
    c2_solution = sol.c2;

    %The final solution
    problemSolution = exp(alpha_part * x) * (c1_solution * cos(beta_part * x) + c2_solution * sin(beta_part * x));
    

end

function dSolveSolution = solveEquationWithDSolve(a, b, c, y1, y2)
    
    syms y(x)
    Dy = diff(y);

    ode_eqn = a * diff(y,x,2) + b * diff(y,x) + c * y == 0;

    cond1 = y(0) == y1;
    cond2 = Dy(0) == y2;

    dSolveSolution = dsolve(ode_eqn, [cond1, cond2]);


end
