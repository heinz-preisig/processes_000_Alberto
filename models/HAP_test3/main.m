function output_vars = main
  % Index sets
  N = [1,2];
  N_lbl = "N";
  A = [1];
  A_lbl = "A";
  I = [];
  I_lbl = "I";
  S = [1];
  S_lbl = "S";
  indexOrder = cell(1, 4);
  indexOrder(1)  = N_lbl;
  indexOrder(2)  = A_lbl;
  indexOrder(3)  = I_lbl;
  indexOrder(4)  = S_lbl;

  % Variables
  % accumulation of molar mass due to convection
  anc = MultiDimVar({N_lbl, S_lbl}, [2,1], indexOrder);


  % density
  rho = MultiDimVar({N_lbl}, [2], indexOrder);
  rho(2) = 15;


  % end time
  te = MultiDimVar({}, [1], indexOrder);
  te(1) = 10;


  % molar convective flow in x-direction
  fnc_x = MultiDimVar({A_lbl, S_lbl}, [1,1], indexOrder);
  fnc_x(1,1) = 5;


  % initial molar mass
  no = MultiDimVar({N_lbl, S_lbl}, [2,1], indexOrder);
  no(2,1) = 0;


  % volume
  V = MultiDimVar({N_lbl}, [2], indexOrder);


  % starting time
  to = MultiDimVar({}, [1], indexOrder);
  to(1) = 0;


  % molecular masses
  Mm = MultiDimVar({S_lbl}, [1], indexOrder);
  Mm(1) = 10;


  % time
  t = MultiDimVar({}, [1], indexOrder);


  % cross sectional area xy
  Axy = MultiDimVar({N_lbl}, [2], indexOrder);


  % accumulation due to diffusion in x-direction
  and_x = MultiDimVar({N_lbl, S_lbl}, [2,1], indexOrder);
  and_x(2,1) = 0;


  % mass
  m = MultiDimVar({N_lbl}, [2], indexOrder);


  %
  pi = MultiDimVar({}, [1], indexOrder);
  pi(1) = 3.14;


  % fundamental state -- molar mass
  n = MultiDimVar({N_lbl, S_lbl}, [2,1], indexOrder);


  % numerical value one half
  oneHalf = MultiDimVar({}, [1], indexOrder);
  oneHalf(1) = 0.5;


  % liquid level cylinder
  h_L = MultiDimVar({N_lbl}, [2], indexOrder);


  % incidence matrix
  F = MultiDimVar({N_lbl, A_lbl}, [2,1], indexOrder);
  F(1,1) = -1;
  F(2,1) = 1;


  % diameter_x
  d_x = MultiDimVar({N_lbl}, [2], indexOrder);
  d_x(2) = 2;


  % differential mass balance without reaction
  an = MultiDimVar({N_lbl, S_lbl}, [2,1], indexOrder);



  % Integrators
  % Index sets
  N_E_93 = [2];
  S_E_93 = [1];
  % Initial conditions
  phi_0(1:1) = reshape((no(N_E_93, S_E_93)).value, [], 1);

  % Integration interval
  integration_interval = [to.value, te.value];
  dt = 1;
  options = odeset('InitialStep', dt, 'MaxStep', dt, 'RelTol', 1e-3, 'AbsTol', 1e-3);
  % Integrator routine
  [t, phi] = ode45(@f, integration_interval, phi_0, options);

  % Printing output
  % HEADERS
  fprintf("%s", "time");
  fprintf("\n");

  % DATA
  time_steps = length(t);
  for i = 1:time_steps
    [integrand, other_vars] = f(t(i), phi(i,:));
    fprintf("%d\t", t(i));
    fprintf("%d\t", other_vars);
    fprintf("\n");
  endfor

  % Integrands
  function [dphidt, extra_output] = f(t, phi)
    n(N_E_93, S_E_93) = reshape(phi(1:1), 1,1);

    A_E_87 = [1];
    S_E_87 = [1];
    N_E_87 = [2];
    anc(N_E_87, S_E_87) = E_87(F, fnc_x, N_E_87, A_E_87, S_E_87);

    N_E_30 = [2];
    S_E_30 = [1];
    m(N_E_30) = E_30(n, Mm, N_E_30, S_E_30);

    N_E_92 = [2];
    S_E_92 = [1];
    an(N_E_92, S_E_92) = E_92(and_x, anc, N_E_92, S_E_92);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sys1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sol = sys1(rho, V, m);

    N_E_97 = [2];
    S_E_97 = [1];
    V(N_E_97) = sol(1:1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N_E_95 = [2];
    S_E_95 = [1];
    Axy(N_E_95) = E_95(oneHalf, d_x, pi, N_E_95, S_E_95);

    N_E_96 = [2];
    S_E_96 = [1];
    h_L(N_E_96) = E_96(V, Axy, N_E_96, S_E_96);

    extra_output = [
    ];

    dphidt(1:1) = reshape((an(N_E_93, S_E_93)).value, [], 1);
  endfunction
endfunction

% Functions for the equations
function sol = E_87(F, fnc_x, N, A, S)
  sol = einsum(F(N, A), fnc_x(A, S), {"A"});
endfunction

function sol = E_30(n, Mm, N, S)
  sol = einsum(Mm(S), n(N, S), {"S"});
endfunction

function sol = E_92(and_x, anc, N, S)
  sol = anc(N, S) + and_x(N, S);
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sys1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol = sys1(rho, V, m)
  % Load initial conditions for the solvers
  run("init.m");
  initial_guess = sys1_init;
  sol = fsolve(@sys, initial_guess);

  function F = sys(x)
    N_E_97 = [2];
    S_E_97 = [1];
    V(N_E_97) = x(1:1);

    F(1:1) = res_E_97(rho, V, m, N_E_97, S_E_97);
  endfunction
endfunction

% Residue equations
function sol = res_E_97(rho, V, m, N, S)
  sol = rho(N) - (einsum(1 ./ V(N), m(N)));

  sol = sol.value;
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sol = E_95(oneHalf, d_x, pi, N, S)
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
