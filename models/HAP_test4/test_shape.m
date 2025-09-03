

%% checking shapes

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

  % initial molar mass
  no = MultiDimVar({N_lbl, S_lbl}, [2,2], indexOrder);
  no(2,2) = 0;
  no(2,1) = 0;


  N_E_93 = [2];
  S_E_93 = [1,2];
  % Initial conditions
  %phi_0(1:2) = reshape((no(N_E_93, S_E_93)).value, [], 1);
  phi_0(1:2) = reshape(no.value(N_E_93, S_E_93), [], 1);

  phi_0


