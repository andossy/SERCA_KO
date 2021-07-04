function [p, x0, varargout]  = KOInit(varargin) %loginfo, names,
                                                  %currents
% Initial conditions and parameters 
% to Mouse SERCA-KO model from Li et al. 2012
% Implemented by Andy Edwards April 5 2013

  if nargout < 2
    error('To few output arguments.')
  end
  nxargout = max(nargout-2,0);
  
%-------------------------------------------------------------------------------
% Constants
%-------------------------------------------------------------------------------

p.K_L = 0.6;   % micromolar (in L_type_calcium_current)
p.LCC_C1 = 33.0;   % millivolt (in L_type_calcium_current)
p.LCC_C2 = 8.2;   % millivolt (in L_type_calcium_current)
p.LCC_C3 = 0.1;   % dimensionless (in L_type_calcium_current)
p.LCC_C4 = 40.0;   % millivolt (in L_type_calcium_current)
p.LCC_C5 = 6.0;   % millivolt (in L_type_calcium_current)
p.LCC_C6 = 10.0;   % millisecond (in L_type_calcium_current)
p.LCC_C7 = 315.0;   % millisecond (in L_type_calcium_current)
p.LCC_C8 = 26.0;   % millivolt (in L_type_calcium_current)
p.LCC_C9 = 4.5;   % millivolt (in L_type_calcium_current)
p.P_CaL = 2.0e-7;   % microlitre_per_millisecond (in L_type_calcium_current)
p.V_L = -5.0;   % millivolt (in L_type_calcium_current)
p.a = 0.3;   % dimensionless (in L_type_calcium_current)
p.b = 0.4;   % dimensionless (in L_type_calcium_current)
p.delta_V_L = 8.0;   % millivolt (in L_type_calcium_current)
p.i_CaL_max = 20.0;   % picoA_per_picoF (in L_type_calcium_current)
p.phi_L = 2.5;   % dimensionless (in L_type_calcium_current)
p.t_L = 9.0;   % millisecond (in L_type_calcium_current)
p.tau_L = 400.0;   % millisecond (in L_type_calcium_current)
p.anion_o = 0.0;   % micromolar (in anion_flux)
p.p_anion = 1.0e-7;   % cm_per_second (in anion_flux)
p.E_Cl = -40.0;   % millivolt (in calcium_activated_chloride_current)
p.Km_Cl = 10.0;   % micromolar (in calcium_activated_chloride_current)
p.g_ClCa = 10.0;   % milliS_per_microF (in calcium_activated_chloride_current)
p.g_Cab = 0.0004;   % milliS_per_microF (in calcium_background_current)
p.Bmax = 109.0;   % micromolar (in calcium_concentration)
p.CSQN_tot = 50000.0;   % micromolar (in calcium_concentration)
p.Kd_buffer = 0.6;   % micromolar (in calcium_concentration)
p.Km_CSQN = 630.0;   % micromolar (in calcium_concentration)
p.A_PRyR = -0.2;   % per_millisecond (in calcium_fluxes)
p.B_PRyR = -3.0;   % per_millisecond (in calcium_fluxes)
p.Km_up = 0.4928;   % micromolar (in calcium_fluxes)
p.V_leak = 2.0e-5;   % per_millisecond (in calcium_fluxes)
p.V_rel = 4.5;   % per_millisecond (in calcium_fluxes)
p.Vup = 0.07;   % micromolar_per_millisecond (in calcium_fluxes)
p.nH = 2.0;   % dimensionless (in calcium_fluxes)
p.tau_Ca_jcsl_const = 2.108e-7;   % microlitre_per_millisecond (in calcium_fluxes)
p.tau_ca_slcyt_const = 3.25e-6;   % microlitre_per_millisecond (in calcium_fluxes)
p.tau_tr = 20.0;   % millisecond (in calcium_fluxes)
p.Km_PMCA = 0.3506;   % micromolar (in calcium_pump_current)
p.i_PMCA_max = 1.1231;   % picoA_per_picoF (in calcium_pump_current)
p.Atot = 0.000148;   % cm2 (in cell)
p.Cao = 1000.0;   % micromolar (in cell)
p.Clo = 148400.0;   % micromolar (in cell)
p.Cm = 1.0;   % microF_per_cm2 (in cell)
p.F = 96.5;   % coulomb_per_millimole (in cell)
p.Ko = 5400.0;   % micromolar (in cell)
p.Nao = 140000.0;   % micromolar (in cell)
p.R = 8.314;   % joule_per_mole_kelvin (in cell)
p.T = 310.0;   % kelvin (in cell)
p.TT_ratio = 0.19;   % dimensionless (in cell)
p.VJSR = 7.7e-8;   % microlitre (in cell)
p.VNSR = 2.31e-7;   % microlitre (in cell)
p.Vjc = 2.2e-8;   % microlitre (in cell)
p.Vmyo = 2.2e-5;   % microlitre (in cell)
p.Vsl = 4.4e-7;   % microlitre (in cell)
p.prepulses_number = 1000000.0;   % dimensionless (in cell)
p.stim_amplitude = -50.0;   % picoA_per_picoF (in cell)
p.stim_duration = 3.0;   % millisecond (in cell)
p.stim_offset = 100.0;   % millisecond (in cell)
p.stim_period = 1000.0;   % millisecond (in cell)
p.p_cl = 1.3e-8;   % cm_per_second (in chloride_current_constant_field)
p.g_Na = 16.0;   % milliS_per_microF (in fast_sodium_current)
p.alpha_a_const = 45.0;   % millivolt (in fast_transient_outward_K_I)
p.alpha_b_const = 45.0;   % millivolt (in fast_transient_outward_K_I)
p.g_Kto_f = 0.535;   % milliS_per_microF (in fast_transient_outward_K_I)
p.g_Kss = 0.06;   % milliS_per_microF (in non_inactivating_steady_state_K_I)
p.B1 = 84200.0;   % micromolar (in pH_regulation)
p.B2 = 29400.0;   % micromolar (in pH_regulation)
p.J_acid_flux_in = 0.045;   % micromolar_per_millisecond (in pH_regulation)
p.K_Cl_che = 18000000.0;   % micromolar (in pH_regulation)
p.K_H_nhe = 0.0001778;   % micromolar (in pH_regulation)
p.K_Na_nhe = 21490.0;   % micromolar (in pH_regulation)
p.K_i_nhe = 0.4111;   % micromolar (in pH_regulation)
p.alpha_f_1_nhe = 0.00197;   % per_millisecond (in pH_regulation)
p.alpha_r_2_nhe = 0.000818;   % per_millisecond (in pH_regulation)
p.k_f_1_che = 0.00429;   % per_millisecond (in pH_regulation)
p.k_f_2_che = 0.0681;   % per_millisecond (in pH_regulation)
p.k_f_2_nhe = 0.01415;   % per_millisecond (in pH_regulation)
p.k_r_1_che = 0.25;   % per_millisecond (in pH_regulation)
p.k_r_1_nhe = 1.1724;   % per_millisecond (in pH_regulation)
p.nH_nhe = 2.9053;   % dimensionless (in pH_regulation)
p.pHo = 7.4;   % dimensionless (in pH_regulation)
p.pK1 = 6.03;   % dimensionless (in pH_regulation)
p.pK2 = 7.57;   % dimensionless (in pH_regulation)
p.pK_H_che = 7.95;   % dimensionless (in pH_regulation)
p.g_Kr = 0.0165;   % milliS_per_microF (in rapid_delayed_rectifier_K_I)
p.kb = 0.036778;   % per_millisecond (in rapid_delayed_rectifier_K_I)
p.kf = 0.023761;   % per_millisecond (in rapid_delayed_rectifier_K_I)
p.k_minus_a = 0.07125;   % per_millisecond (in ryanodine_receptors)
p.k_minus_b = 0.965;   % per_millisecond (in ryanodine_receptors)
p.k_minus_c = 0.0008;   % per_millisecond (in ryanodine_receptors)
p.k_plus_a = 6.075e-6;   % micromolar4_per_millisecond (in ryanodine_receptors)
p.k_plus_b = 4.05e-6;   % micromolar3_per_millisecond (in ryanodine_receptors)
p.k_plus_c = 0.009;   % per_millisecond (in ryanodine_receptors)
p.m = 3.0;   % dimensionless (in ryanodine_receptors)
p.n = 4.0;   % dimensionless (in ryanodine_receptors)
p.g_Ks = 0.00575;   % milliS_per_microF (in slow_delayed_rectifier_K_I)
p.g_Nab = 0.0026;   % milliS_per_microF (in sodium_background_current)
p.I_NCX_max = 6.4116;   % picoA_per_picoF (in sodium_calcium_exchange_current)
p.K_mAllo = 0.0;   % micromolar (in sodium_calcium_exchange_current)
p.K_mCai = 3.6;   % micromolar (in sodium_calcium_exchange_current)
p.K_mCao = 1400.0;   % micromolar (in sodium_calcium_exchange_current)
p.K_mNai = 12000.0;   % micromolar (in sodium_calcium_exchange_current)
p.K_mNao = 88000.0;   % micromolar (in sodium_calcium_exchange_current)
p.eta = 0.35;   % dimensionless (in sodium_calcium_exchange_current)
p.k_sat = 0.27;   % dimensionless (in sodium_calcium_exchange_current)
p.Bmax_SL = 1650.0;   % micromolar (in sodium_concentration)
p.Bmax_jct = 3700.0;   % micromolar (in sodium_concentration)
p.koff = 0.001;   % per_millisecond (in sodium_concentration)
p.kon = 1.0e-7;   % per_micromolar_per_millisecond (in sodium_concentration)
p.tau_Na_jcsl_const = 1.39e-6;   % microlitre_per_millisecond (in sodium_concentration)
p.tau_Na_slcyt_const = 4.77e-5;   % microlitre_per_millisecond (in sodium_concentration)
p.Km_Ko = 1500.0;   % micromolar (in sodium_potassium_pump_current)
p.Km_Na_alpha1 = 21000.0;   % micromolar (in sodium_potassium_pump_current)
p.Km_Na_alpha2 = 21000.0;   % micromolar (in sodium_potassium_pump_current)
p.i_NKA_max_alpha1 = 4.4;   % picoA_per_picoF (in sodium_potassium_pump_current)
p.i_NKA_max_alpha2 = 0.9;   % picoA_per_picoF (in sodium_potassium_pump_current)
p.nH_NKAalpha1 = 3.0;   % dimensionless (in sodium_potassium_pump_current)
p.nH_NKAalpha2 = 3.0;   % dimensionless (in sodium_potassium_pump_current)
p.g_K1 = 0.35;   % milliS_per_microF (in time_independent_K_I)
% old i_Kur
% p.g_Kur = 0.45;   % milliS_per_microF (in ultra_rapidly_activating_delayed_rectifier_K_I)
% p.tau_a_const = 2.058;   % millisecond (in ultra_rapidly_activating_delayed_rectifier_K_I)
% p.tau_i_const = 643.0;   % millisecond (in ultra_rapidly_activating_delayed_rectifier_K_I)

% new i_Kur
p.g_Kur1 = 0.25;   % milliS_per_microF (in ultra_rapidly_activating_delayed_rectifier_K_I)
p.g_Kur2 = 0.20;   % milliS_per_microF (in ultra_rapidly_activating_delayed_rectifier_K_I)

%--------- Pace current parameters
p.pacestart = 50.0;
p.pacedur = 0.5;
p.paceamp = -50.0;

%---------Voltage Clamp current

%p.VClamp_n = 1;
p.VClampAmp = [];
p.VClampTimes = [];
p.VClampR = 50; % In Mohm

% Y = [20524.4225224089, 0.39648348828885, 0.1870667990772, 2.51368630484038e-5, 3.39087179979051e-5, 0.975309150371633, 1.12948431110931e-34, 284.399331749045, 320.282688446674, 0.0948356862379185, 0.122423236690487, 0.0975652924644621, 1.91667709962937e-6, -78.3418282996005, 0.000538721428676182, 0.0250882080700854, 1.69319147238874e-5, 2.40166758551435e-5, 0.023628100999744, 0.460857600573785, 0.000507376796911785, 2.3601674345502e-6, 0.0170221330728053, 0.995453949127464, 0.145968004744334, 1.0, 7.11872164143083, 107013.338083988, 0.00133123745592853, 0.000993143374968524, 0.00125492962155179, 0.0103988907946101, 0.144694335439345, 0.0014553162561183, 1.12097329264715e-11, 0.00158518129158872, 9955.89761461143, 9955.72777296044, 1846.43649618642, 9955.88479546724, 823.416316861962, 0.0404038645845194, 0.973465329670617];
% YNames = {'Cli', 'Ijc', 'Isl', 'Ojc', 'Osl', 'y_gate', 'anion_i', 'CaJSR', 'CaNSR', 'Cai', 'Cajc', 'Casl', 'P_RyR', 'V', 'C_Na1', 'C_Na2', 'I1_Na', 'I2_Na', 'IC_Na2', 'IC_Na3', 'IF_Na', 'O_Na', 'ato_f', 'ito_f', 'aKss', 'iKss', 'pH_i', 'Ki', 'C_K1', 'C_K2', 'I_K', 'O_K', 'P_C2', 'P_O1', 'P_O2', 'nKs', 'Nai', 'Najc', 'Najc_buf', 'Nasl', 'Nasl_buf', 'aur', 'iur'};
%--------- end parameters ------------------------------

%---------- Initial conditions----------------------------------------

x0.Cli =15569.7457;
x0.Ijc =0.36275;
x0.Isl =0.16548;
x0.Ojc =3.2893e-05;
x0.Osl =4.3146e-05;
x0.y_gate =0.97475;
x0.anion_i =22565.5576;
x0.CaJSR =51.0285;
x0.CaNSR =67.1904;
x0.Cai =0.12351;
x0.Cajc =0.10913;
x0.Casl =0.10782;
x0.P_RyR =1.4044e-06;
x0.V =-76.6253;
x0.C_Na1 =0.00069288;
x0.C_Na2 =0.026511;
x0.I1_Na =0.00011712;
x0.I2_Na =0.00010156;
x0.IC_Na2 =0.031338;
x0.IC_Na3 =0.50944;
x0.IF_Na =0.0008191;
x0.O_Na =3.8075e-06;
x0.ato_f =0.020083;
x0.ito_f =0.99361;
x0.aKss =0.16532;
x0.iKss =1;
x0.pH_i =6.8533;
x0.Ki =109634.3309;
x0.C_K1 =0.0015002;
x0.C_K2 =0.0013072;
x0.I_K =0.0024688;
x0.O_K =0.018359;
x0.P_C2 =0.19354;
x0.P_O1 =0.0019466;
x0.P_O2 =1.0622e-11;
x0.nKs =0.0028004;
x0.Nai =14422.4689;
x0.Najc =14422.0543;
x0.Najc_buf =2185.4578;
x0.Nasl =14422.4306;
x0.Nasl_buf =974.6053;
% old i_Kur
% x0.aur =0.040404;
% x0.iur =0.97347;

% new i_Kur
x0.aur1 =0.040404;
x0.aur2 =0.5; 
x0.iur1 =0.97347;
x0.iur2 =1.1;

if exist('inputstruct','var')
    if length(inputstruct)>0
        for field=inputfields
          field = field{1};
          if isfield(x0,field)
            x0.(field) = inputstruct.(field);
          end
        end
    end
end

%---------- End Initial conditions----------------------------------------

%-------------- Logging parameters--------------------------------
names.modelname = 'Li 2012 KO model';

names.states = fieldnames(x0)';
% Place the ionexchangers last, ie 'INaCa','INaK'. (For logging purposes)
names.currents = {'i_ClCa','i_CaL_jc','i_CaL_sl','i_CaL','i_anion', ...
                 'i_Cab_jc','i_Cab_sl','i_Cab','i_PMCA_sl','i_NCX_sl','i_PMCA_jc',...
                 'i_NCX_jc','i_PMCA','i_PMCA_proton','i_Stim','i_NCX','i_Na_jc','i_Na_sl',...
                 'i_Na','i_Nab_jc','i_Nab_sl','i_Nab','i_NKA_alpha1_jc','i_NKA_alpha2_jc',...
                 'i_NKA_jc','i_NKA_alpha1_sl','i_NKA_alpha2_sl','i_NKA_sl','i_NKA','i_Kto_f',...
                 'i_K1','i_Ks','i_Kur','i_Kur1','i_Kur2','i_Kss','i_Kr','i_NKA_alpha1','i_NKA_alpha2','IVClamp'};

names.ionfluxes = {'Jche','Jae','J_Cab','JCa_sl','JCa_jc','JCa_slcyt',...
                 'JCa_jcsl','J_rel','J_leak','J_serca','J_tr','J_netup',...
                 'J_PMCA_sl','J_PMCA_jc','J_PMCA_total','J_PMCA_proton','JCO2',...
                 'Jnbc','Jnhe','Jhyd','J_NCX_sl','J_NCX_jc','J_NCX_total',...
                 'JNa_slcyt','JNa_jcsl','J_Na_NKA'};
names.buffers = {'Bsl','Bjc','Bi','BJSR'};
names.rates = {};
names.hiddenstates = {};

% The possible logging variables. Use same order as in Bond.m file.
possible_logvar = {'state_der','currents','ionfluxes','buffers','rates'};%,'hiddenstates'};


loginfo = [];
if nargin > 0 
  argin = cell2cell(varargin);
  % If ioninfo to be logged then add both buffers and currents as logging
  % parameters
  if ismemb('ioninfo',argin)
    if nxargout < 3
      error('If logging ioninfo, then 5 or more output variables are needed.')
    end
    if ~ismemb('buffers',argin)
      argin = [argin {'buffers'}];
    end
    if ~ismemb('currents',argin)
      argin = [argin {'currents'}];
    end
    if ~ismemb('ionfluxes',argin)
      argin = [argin {'ionfluxes'}];
    end
  end
  [tmp,ind] = ismemb(argin,possible_logvar);
  logvar = possible_logvar(sort(ind(ind>0)));
  logging = length(logvar) > 0;
  for i=1:length(logvar)
    if strcmp(logvar{i},'state_der')
      field = 'states';
    else
      field = logvar{i};
    end
    loginfo.(logvar{i}) = length(names.(field));
  end
else
  logging = 0;
  loginfo = [];
end

if ~logging
end

% Explicit set the different ion compartments
ions.types = {'Ca','Na','K'};
ionvalence = {2 1 1};

statenames.Ca = {'Cai','Cajc','Casl','CaJSR','CaNSR'};
statenames.Na = {'Nai','Najc','Nasl',};
statenames.K = {'Ki'};

buffers.Ca = {'Bsl','Bjc','Bi','BJSR'};
buffers.Na = {};
buffers.K = {};

currents.Ca = {'i_CaL','i_Cab','i_PMCA'};
currents.Na = {'i_Na','i_Nab'};
currents.K =  {'i_Kto_f','i_K1','i_Ks','i_Kur','i_Kss','i_Kr'};

% Ion contributions from special currents (i.e. exchangers and pumps)
currents.special.Ca.i_NaCa = -2;
currents.special.Na.i_NaCa = 3;
currents.special.Na.i_NKA = 3;
currents.special.K.i_NKA = -2;

% Iterate through the ions and set the differents parameters
for i = 1:length(ions.types)
  ion = ions.types{i};
  % Setting the indexes to ion compartments
%  ions.ind.(ion).states = indfind(names.states,statenames.(ion));
  % Setting volumes to the corresponding ion compartments
  for n=1:length(statenames.(ion))
    statename = statenames.(ion){n};
    ions.conc.states.(ion).(statename) = p.(volumestr(statename));
  end
  for n=1:length(buffers.(ion))
    buffername = buffers.(ion){n};
    ions.conc.buffers.(ion).(buffername) = p.(volumestr(buffername));
  end
  % Setting ioncurrents together with their corresponding conversion
  % factor
  for n=1:length(currents.(ion))
    curname = currents.(ion){n};
    ions.currents.(ion).(curname) = p.Atot/(ionvalence{i}*p.F);
  end
  special_names = fieldnames(currents.special.(ion));
  for n=1:length(special_names)
    curname = special_names{n};
    ions.currents.(ion).(curname) =currents.special.(ion).(curname)*...
        p.Atot/(ionvalence{i}*p.F);
  end
end

% Index of currents that are super threshold and subthreshold
suprth = [1 5 6 9 10 14];
currents.suprth = suprth;
currents.subth = setdiff([1:1:15],suprth );



% Assigning out variables
if nxargout >= 1
  varargout{1} = loginfo;
end
if nxargout >=2
  varargout{2} = names;
end
if nxargout >= 3 
  varargout{3} = ions;
end
if nxargout >=4
  varargout{4} = currents;
end

%******************************************************
function out = volumestr(instr)
% VOLUMESTR Return a volume corresponding to the compartment name instr.
  substrings = {'jc','JSR','sl','myo'};

  for i=1:length(substrings)
    ind = findstr(substrings{i},instr);
    found = length(ind)>0;
    if found
      break;
    end
  end
  
  if found
    out = ['V' substrings{i}];
  else
    out = 'Vmyo';
  end

