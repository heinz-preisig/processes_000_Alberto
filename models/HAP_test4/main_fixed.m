function output_vars = main_fixed
  % Index sets
  N = [1,2];
  N_lbl = "N";
  A = [1];
  A_lbl = "A";
  I = [];
  I_lbl = "I";
  S = [1,2];
  S_lbl = "S";
  indexOrder = cell(1, 4);
  indexOrder(1)  = N_lbl;
  indexOrder(2)  = A_lbl;
  indexOrder(3)  = I_lbl;
  indexOrder(4)  = S_lbl;

  % Variables
  % molecular masses
  Mm = MultiDimVar({S_lbl}, [2], indexOrder);
  Mm(2) = 1.5;
  Mm(1) = 1.5;


  % molar convective flow in x-direction
  fnc_x = MultiDimVar({A_lbl, S_lbl}, [1,2], indexOrder);
  fnc_x(1,2) = 10;
  fnc_x(1,1) = 10;


  % incidence matrix
  F = MultiDimVar({N_lbl, A_lbl}, [2,1], indexOrder);
  F(1,1) = -1;
  F(2,1) = 1;


  % initial molar mass
  no = MultiDimVar({N_lbl, S_lbl}, [2,2], indexOrder);
  no(2,2) = 0;
  no(2,1) = 0;


  % diameter_x
  d_x = MultiDimVar({N_lbl}, [2], indexOrder);
  d_x(2) = 1;


  % accumulation due to diffusion in x-direction
  and_x = MultiDimVar({N_lbl, S_lbl}, [2,2], indexOrder);
  and_x(2,2) = 0;
  and_x(2,1) = 0;


  % density
  rho = MultiDimVar({N_lbl}, [2], indexOrder);
  rho(2) = 1000;


  % liquid level cylinder
  h_L = MultiDimVar({N_lbl}, [2], indexOrder);


  %
  pi = MultiDimVar({}, [1], indexOrder);
  pi(1) = 1.3;


  % starting time
  to = MultiDimVar({}, [1], indexOrder);
  to(1) = 0;


  % numerical value one half
  oneHalf = MultiDimVar({}, [1], indexOrder);
  oneHalf(1) = 0.5;


  % cross sectional area xy
  Axy = MultiDimVar({N_lbl}, [2], indexOrder);


  % mass
  m = MultiDimVar({N_lbl}, [2], indexOrder);


  % volume
  V = MultiDimVar({N_lbl}, [2], indexOrder);


  % end time
  te = MultiDimVar({}, [1], indexOrder);
  te(1) = 10;


  % differential mass balance without reaction
  an = MultiDimVar({N_lbl, S_lbl}, [2,2], indexOrder);


  % time
  t = MultiDimVar({}, [1], indexOrder);


  % fundamental state -- molar mass
  n = MultiDimVar({N_lbl, S_lbl}, [2,2], indexOrder);


  % accumulation of molar mass due to convection
  anc = MultiDimVar({N_lbl, S_lbl}, [2,2], indexOrder);



  % Integrators
  % Index sets
  N_E_93 = [2];
  S_E_93 = [1,2];
  % Initial conditions
  %phi_0(1:2) = reshape((no(N_E_93, S_E_93)).value, [], 1);
  phi_0 = reshape(no.value(N_E_93, S_E_93), [], 1);

  % Integration interval
  integration_interval = [to.value, te.value];
  dt = 5;
  options = odeset('InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-3, 'AbsTol', 1e-3);
  % Integrator routine
  [t, phi] = ode45(@f, integration_interval, phi_0, options);
##  [t, phi] = ode45(@f, integration_interval, phi_0);

  % Printing output
  % HEADERS
  fprintf("%s", "time");
  fprintf("\n");

  % DATA
  time_steps = length(t);
  fileID = fopen('fixed_results.txt','w')
  for i = 1:time_steps
    [integrand, other_vars] = f(t(i), phi(i,:));
    fprintf("%d\t", t(i));
    fprintf("%d\t", other_vars);
    fprintf("\n");
    fprintf(fileID,"%d\t", t(i));
    fprintf(fileID,"%d\t", other_vars);
    fprintf(fileID,"\n");
  endfor
  fclose(fileID)

  % Integrands
  function [dphidt, extra_output] = f(t, phi)
    n(N_E_93, S_E_93) = reshape(phi(1:2), 1,2);

    A_E_87 = [1];
    S_E_87 = [1,2];
    N_E_87 = [2];
    anc(N_E_87, S_E_87) = E_87(fnc_x, F, N_E_87, A_E_87, S_E_87);

    N_E_30 = [2];
    S_E_30 = [1,2];
    m(N_E_30) = E_30(n, Mm, N_E_30, S_E_30);

    N_E_92 = [2];
    S_E_92 = [1,2];
    an(N_E_92, S_E_92) = E_92(and_x, anc, N_E_92, S_E_92);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sys1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sol = sys1(V, rho, m);

    N_E_97 = [2];
    S_E_97 = [1,2];
    V(N_E_97) = sol(1:1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N_E_95 = [2];
    S_E_95 = [1,2];
    Axy(N_E_95) = E_95(oneHalf, pi, d_x, N_E_95, S_E_95);

    N_E_96 = [2];
    S_E_96 = [1,2];
    h_L(N_E_96) = E_96(V, Axy, N_E_96, S_E_96);

    extra_output = [ V(N_E_97).value, h_L(N_E_96).value
    ]; % HAP fix

    dphidt(1:2) = reshape(an.value(N_E_93, S_E_93), [], 1);
  endfunction
endfunction

% Functions for the equations
function sol = E_87(fnc_x, F, N, A, S)
  sol = einsum(F(N, A), fnc_x(A, S), {"A"});
endfunction

function sol = E_30(n, Mm, N, S)
  sol = einsum(Mm(S), n(N, S), {"S"});
endfunction

function sol = E_92(and_x, anc, N, S)
  sol = anc(N, S) + and_x(N, S);
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sys1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol = sys1(V, rho, m)
  % Load initial conditions for the solvers
  run("init_fixed.m");   % HAP fix
  initial_guess = sys1_init;
  sol = fsolve(@sys, initial_guess);

  function F = sys(x)
    N_E_97 = [2];
    S_E_97 = [1,2];
    V(N_E_97) = x;

    F(1:2) = res_E_97(V, rho, m, N_E_97, S_E_97);
  endfunction
endfunction

% Residue equations
function sol = res_E_97(V, rho, m, N, S)
  sol = rho(N) - (einsum(1 ./ V(N), m(N)));

  sol = sol.value;
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = E_95(oneHalf, pi, d_x, N, S)
  sol = einsum(einsum(pi, (einsum(oneHalf, d_x(N)))), (einsum(oneHalf, d_x(N))));
endfunction

function sol = E_96(V, Axy, N, S)
  sol = einsum(V(N), 1 ./ Axy(N));
endfunction


% Auxiliar functions
function result = indexunion(varargin)
  result = varargin{1};
  for i = 2:nargin
    result = union(result, varargin{i});
  endfor
endfunction
