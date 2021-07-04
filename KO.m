%===============================================================================
% CellML file:   C:\Users\Andy\Documents\020911_FF_NaKsubspace_SS_b_for_thesis.cellml
% CellML model:  MyModel_2006
% Date and time: 4/5/2013 at 5:57:21 at AM
%-------------------------------------------------------------------------------
% Conversion from CellML 1.0 to MATLAB (init) was done using COR (0.9.31.1409)
%    Copyright 2002-2013 Dr Alan Garny
%    http://cor.physiol.ox.ac.uk/ - cor@physiol.ox.ac.uk
%-------------------------------------------------------------------------------
% http://www.cellml.org/
%===============================================================================
% Remaining hacks were implemented by Andy Edwards (April 10th 2013) to permit 
% voltage clamp, logging of all calculated variables, and variable pacing protocols
% by the pacemodel and paceone scripts. The initialization file for this 
% model is KOInit.

function [dY,varargout] = KO(time, Y,p, loginfo)

%-------------------------------------------------------------------------------
% Initial conditions
%-------------------------------------------------------------------------------

% Y = [20524.4225224089, 0.39648348828885, 0.1870667990772, 2.51368630484038e-5, 3.39087179979051e-5, 0.975309150371633, 1.12948431110931e-34, 284.399331749045, 320.282688446674, 0.0948356862379185, 0.122423236690487, 0.0975652924644621, 1.91667709962937e-6, -78.3418282996005, 0.000538721428676182, 0.0250882080700854, 1.69319147238874e-5, 2.40166758551435e-5, 0.023628100999744, 0.460857600573785, 0.000507376796911785, 2.3601674345502e-6, 0.0170221330728053, 0.995453949127464, 0.145968004744334, 1.0, 7.11872164143083, 107013.338083988, 0.00133123745592853, 0.000993143374968524, 0.00125492962155179, 0.0103988907946101, 0.144694335439345, 0.0014553162561183, 1.12097329264715e-11, 0.00158518129158872, 9955.89761461143, 9955.72777296044, 1846.43649618642, 9955.88479546724, 823.416316861962, 0.0404038645845194, 0.973465329670617];

% YNames = {'Cli', 'Ijc', 'Isl', 'Ojc', 'Osl', 'y_gate', 'anion_i', 'CaJSR', 'CaNSR', 'Cai', 'Cajc', 'Casl', 'P_RyR', 'V', 'C_Na1', 'C_Na2', 'I1_Na', 'I2_Na', 'IC_Na2', 'IC_Na3', 'IF_Na', 'O_Na', 'ato_f', 'ito_f', 'aKss', 'iKss', 'pH_i', 'Ki', 'C_K1', 'C_K2', 'I_K', 'O_K', 'P_C2', 'P_O1', 'P_O2', 'nKs', 'Nai', 'Najc', 'Najc_buf', 'Nasl', 'Nasl_buf', 'aur', 'iur'};
% YUnits = {'micromolar', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'micromolar', 'micromolar', 'micromolar', 'micromolar', 'micromolar', 'micromolar', 'dimensionless', 'millivolt', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'micromolar', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'dimensionless', 'micromolar', 'micromolar', 'micromolar', 'micromolar', 'micromolar', 'dimensionless', 'dimensionless'};
% YComponents = {'Cl_concentration', 'L_type_calcium_current', 'L_type_calcium_current', 'L_type_calcium_current', 'L_type_calcium_current', 'L_type_calcium_current', 'anion_flux', 'calcium_concentration', 'calcium_concentration', 'calcium_concentration', 'calcium_concentration', 'calcium_concentration', 'calcium_fluxes', 'cell', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_sodium_current', 'fast_transient_outward_K_I', 'fast_transient_outward_K_I', 'non_inactivating_steady_state_K_I', 'non_inactivating_steady_state_K_I', 'pH_regulation', 'potassium_concentration', 'rapid_delayed_rectifier_K_I', 'rapid_delayed_rectifier_K_I', 'rapid_delayed_rectifier_K_I', 'rapid_delayed_rectifier_K_I', 'ryanodine_receptors', 'ryanodine_receptors', 'ryanodine_receptors', 'slow_delayed_rectifier_K_I', 'sodium_concentration', 'sodium_concentration', 'sodium_concentration', 'sodium_concentration', 'sodium_concentration', 'ultra_rapidly_activating_delayed_rectifier_K_I', 'ultra_rapidly_activating_delayed_rectifier_K_I'};

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: Cli (micromolar) (in Cl_concentration)
% 2: Ijc (dimensionless) (in L_type_calcium_current)
% 3: Isl (dimensionless) (in L_type_calcium_current)
% 4: Ojc (dimensionless) (in L_type_calcium_current)
% 5: Osl (dimensionless) (in L_type_calcium_current)
% 6: y_gate (dimensionless) (in L_type_calcium_current)
% 7: anion_i (micromolar) (in anion_flux)
% 8: CaJSR (micromolar) (in calcium_concentration)
% 9: CaNSR (micromolar) (in calcium_concentration)
% 10: Cai (micromolar) (in calcium_concentration)
% 11: Cajc (micromolar) (in calcium_concentration)
% 12: Casl (micromolar) (in calcium_concentration)
% 13: P_RyR (dimensionless) (in calcium_fluxes)
% 14: V (millivolt) (in cell)
% 15: C_Na1 (dimensionless) (in fast_sodium_current)
% 16: C_Na2 (dimensionless) (in fast_sodium_current)
% 17: I1_Na (dimensionless) (in fast_sodium_current)
% 18: I2_Na (dimensionless) (in fast_sodium_current)
% 19: IC_Na2 (dimensionless) (in fast_sodium_current)
% 20: IC_Na3 (dimensionless) (in fast_sodium_current)
% 21: IF_Na (dimensionless) (in fast_sodium_current)
% 22: O_Na (dimensionless) (in fast_sodium_current)
% 23: ato_f (dimensionless) (in fast_transient_outward_K_I)
% 24: ito_f (dimensionless) (in fast_transient_outward_K_I)
% 25: aKss (dimensionless) (in non_inactivating_steady_state_K_I)
% 26: iKss (dimensionless) (in non_inactivating_steady_state_K_I)
% 27: pH_i (dimensionless) (in pH_regulation)
% 28: Ki (micromolar) (in potassium_concentration)
% 29: C_K1 (dimensionless) (in rapid_delayed_rectifier_K_I)
% 30: C_K2 (dimensionless) (in rapid_delayed_rectifier_K_I)
% 31: I_K (dimensionless) (in rapid_delayed_rectifier_K_I)
% 32: O_K (dimensionless) (in rapid_delayed_rectifier_K_I)
% 33: P_C2 (dimensionless) (in ryanodine_receptors)
% 34: P_O1 (dimensionless) (in ryanodine_receptors)
% 35: P_O2 (dimensionless) (in ryanodine_receptors)
% 36: nKs (dimensionless) (in slow_delayed_rectifier_K_I)
% 37: Nai (micromolar) (in sodium_concentration)
% 38: Najc (micromolar) (in sodium_concentration)
% 39: Najc_buf (micromolar) (in sodium_concentration)
% 40: Nasl (micromolar) (in sodium_concentration)
% 41: Nasl_buf (micromolar) (in sodium_concentration)
% 42: aur (dimensionless) (in ultra_rapidly_activating_delayed_rectifier_K_I)
% 43: iur (dimensionless) (in ultra_rapidly_activating_delayed_rectifier_K_I)

%-------------------------------------------------------------------------------
% Computed variables
%-------------------------------------------------------------------------------

% Cjc (dimensionless) (in L_type_calcium_current)
% Csl (dimensionless) (in L_type_calcium_current)
% FVRT (dimensionless) (in L_type_calcium_current)
% FVRT_Ca (dimensionless) (in L_type_calcium_current)
% alpha_m (per_millisecond) (in L_type_calcium_current)
% alpha_p (per_millisecond) (in L_type_calcium_current)
% epsilon_m (per_millisecond) (in L_type_calcium_current)
% epsilon_p (per_micromolar_millisecond) (in L_type_calcium_current)
% expVL (dimensionless) (in L_type_calcium_current)
% i_CaL (picoA_per_picoF) (in L_type_calcium_current)
% i_CaL_jc (picoA_per_picoF) (in L_type_calcium_current)
% i_CaL_sl (picoA_per_picoF) (in L_type_calcium_current)
% y_gate_inf (dimensionless) (in L_type_calcium_current)
% y_gate_tau (millisecond) (in L_type_calcium_current)
% i_anion (picoA_per_picoF) (in anion_flux)
% O_ClCa (dimensionless) (in calcium_activated_chloride_current)
% i_ClCa (picoA_per_picoF) (in calcium_activated_chloride_current)
% E_Ca_jc (millivolt) (in calcium_background_current)
% E_Ca_sl (millivolt) (in calcium_background_current)
% J_Cab (micromolar_per_millisecond) (in calcium_background_current)
% i_Cab (picoA_per_picoF) (in calcium_background_current)
% i_Cab_jc (picoA_per_picoF) (in calcium_background_current)
% i_Cab_sl (picoA_per_picoF) (in calcium_background_current)
% BJSR (dimensionless) (in calcium_concentration)
% Bi (dimensionless) (in calcium_concentration)
% Bjc (dimensionless) (in calcium_concentration)
% Bsl (dimensionless) (in calcium_concentration)
% CaSR (micromolar) (in calcium_concentration)
% JCa_jc (micromolar_per_millisecond) (in calcium_concentration)
% JCa_sl (micromolar_per_millisecond) (in calcium_concentration)
% JCa_jcsl (micromolar_per_millisecond) (in calcium_fluxes)
% JCa_slcyt (micromolar_per_millisecond) (in calcium_fluxes)
% J_leak (micromolar_per_millisecond) (in calcium_fluxes)
% J_netup (micromolar_per_millisecond) (in calcium_fluxes)
% J_rel (micromolar_per_millisecond) (in calcium_fluxes)
% J_serca (micromolar_per_millisecond) (in calcium_fluxes)
% J_tr (micromolar_per_millisecond) (in calcium_fluxes)
% P_open (dimensionless) (in calcium_fluxes)
% tau_Ca_jcsl (millisecond) (in calcium_fluxes)
% tau_Ca_slcyt (millisecond) (in calcium_fluxes)
% J_PMCA_jc (micromolar_per_millisecond) (in calcium_pump_current)
% J_PMCA_proton (micromolar_per_millisecond) (in calcium_pump_current)
% J_PMCA_sl (micromolar_per_millisecond) (in calcium_pump_current)
% J_PMCA_total (micromolar_per_millisecond) (in calcium_pump_current)
% i_PMCA (picoA_per_picoF) (in calcium_pump_current)
% i_PMCA_jc (picoA_per_picoF) (in calcium_pump_current)
% i_PMCA_proton (picoA_per_picoF) (in calcium_pump_current)
% i_PMCA_sl (picoA_per_picoF) (in calcium_pump_current)
% Ajc (cm2) (in cell)
% Asl (cm2) (in cell)
% i_Stim (picoA_per_picoF) (in cell)
% past (millisecond) (in cell)
% i_cl (picoA_per_picoF) (in chloride_current_constant_field)
% C_Na3 (dimensionless) (in fast_sodium_current)
% E_Na_jc (millivolt) (in fast_sodium_current)
% E_Na_sl (millivolt) (in fast_sodium_current)
% alpha_Na11 (per_millisecond) (in fast_sodium_current)
% alpha_Na12 (per_millisecond) (in fast_sodium_current)
% alpha_Na13 (per_millisecond) (in fast_sodium_current)
% alpha_Na2 (per_millisecond) (in fast_sodium_current)
% alpha_Na3 (per_millisecond) (in fast_sodium_current)
% alpha_Na4 (per_millisecond) (in fast_sodium_current)
% alpha_Na5 (per_millisecond) (in fast_sodium_current)
% beta_Na11 (per_millisecond) (in fast_sodium_current)
% beta_Na12 (per_millisecond) (in fast_sodium_current)
% beta_Na13 (per_millisecond) (in fast_sodium_current)
% beta_Na2 (per_millisecond) (in fast_sodium_current)
% beta_Na3 (per_millisecond) (in fast_sodium_current)
% beta_Na4 (per_millisecond) (in fast_sodium_current)
% beta_Na5 (per_millisecond) (in fast_sodium_current)
% i_Na (picoA_per_picoF) (in fast_sodium_current)
% i_Na_jc (picoA_per_picoF) (in fast_sodium_current)
% i_Na_sl (picoA_per_picoF) (in fast_sodium_current)
% E_K (millivolt) (in fast_transient_outward_K_I)
% alpha_a (per_millisecond) (in fast_transient_outward_K_I)
% beta_a (per_millisecond) (in fast_transient_outward_K_I)
% i_Kto_f (picoA_per_picoF) (in fast_transient_outward_K_I)
% itof_iss (dimensionless) (in fast_transient_outward_K_I)
% tau_ito_f (millisecond) (in fast_transient_outward_K_I)
% i_Kss (picoA_per_picoF) (in non_inactivating_steady_state_K_I)
% tau_Kss (millisecond) (in non_inactivating_steady_state_K_I)
% Hi (micromolar) (in pH_regulation)
% Ho (micromolar) (in pH_regulation)
% JCO2 (micromolar_per_millisecond) (in pH_regulation)
% Jae (micromolar_per_millisecond) (in pH_regulation)
% Jche (micromolar_per_millisecond) (in pH_regulation)
% Jexch_che (per_millisecond) (in pH_regulation)
% Jhyd (micromolar_per_millisecond) (in pH_regulation)
% Jnbc (micromolar_per_millisecond) (in pH_regulation)
% Jnhe (micromolar_per_millisecond) (in pH_regulation)
% K_OH_che (micromolar) (in pH_regulation)
% OHi (micromolar) (in pH_regulation)
% OHo (micromolar) (in pH_regulation)
% alpha_f_2_nhe (per_millisecond) (in pH_regulation)
% alpha_r_1_nhe (per_millisecond) (in pH_regulation)
% beta_f_1_che (per_millisecond) (in pH_regulation)
% beta_f_2_che (per_millisecond) (in pH_regulation)
% beta_i_1 (micromolar) (beta_i in pH_regulation)
% beta_r_1_che (per_millisecond) (in pH_regulation)
% beta_r_2_che (per_millisecond) (in pH_regulation)
% k_r_2_che (per_millisecond) (in pH_regulation)
% C_K0 (dimensionless) (in rapid_delayed_rectifier_K_I)
% alpha_a0 (per_millisecond) (in rapid_delayed_rectifier_K_I)
% alpha_a1 (per_millisecond) (in rapid_delayed_rectifier_K_I)
% alpha_i (per_millisecond) (in rapid_delayed_rectifier_K_I)
% beta_a0 (per_millisecond) (in rapid_delayed_rectifier_K_I)
% beta_a1 (per_millisecond) (in rapid_delayed_rectifier_K_I)
% beta_i_2 (per_millisecond) (beta_i in rapid_delayed_rectifier_K_I)
% i_Kr (picoA_per_picoF) (in rapid_delayed_rectifier_K_I)
% P_C1 (dimensionless) (in ryanodine_receptors)
% alpha_n (per_millisecond) (in slow_delayed_rectifier_K_I)
% beta_n (per_millisecond) (in slow_delayed_rectifier_K_I)
% i_Ks (picoA_per_picoF) (in slow_delayed_rectifier_K_I)
% i_Nab (picoA_per_picoF) (in sodium_background_current)
% i_Nab_jc (picoA_per_picoF) (in sodium_background_current)
% i_Nab_sl (picoA_per_picoF) (in sodium_background_current)
% J_NCX_jc (micromolar_per_millisecond) (in sodium_calcium_exchange_current)
% J_NCX_sl (micromolar_per_millisecond) (in sodium_calcium_exchange_current)
% J_NCX_total (micromolar_per_millisecond) (in sodium_calcium_exchange_current)
% i_NCX (picoA_per_picoF) (in sodium_calcium_exchange_current)
% i_NCX_jc (picoA_per_picoF) (in sodium_calcium_exchange_current)
% i_NCX_sl (picoA_per_picoF) (in sodium_calcium_exchange_current)
% JNa_jcsl (micromolar_per_millisecond) (in sodium_concentration)
% JNa_slcyt (micromolar_per_millisecond) (in sodium_concentration)
% dNa_SL_buf (micromolar_per_millisecond) (in sodium_concentration)
% dNa_jct_buf (micromolar_per_millisecond) (in sodium_concentration)
% tau_Na_jcsl (millisecond) (in sodium_concentration)
% tau_Na_slcyt (millisecond) (in sodium_concentration)
% J_Na_NKA (micromolar_per_millisecond) (in sodium_potassium_pump_current)
% f_NKA_alpha1 (dimensionless) (in sodium_potassium_pump_current)
% f_NKA_alpha2 (dimensionless) (in sodium_potassium_pump_current)
% i_NKA (picoA_per_picoF) (in sodium_potassium_pump_current)
% i_NKA_alpha1 (picoA_per_picoF) (in sodium_potassium_pump_current)
% i_NKA_alpha1_jc (picoA_per_picoF) (in sodium_potassium_pump_current)
% i_NKA_alpha1_sl (picoA_per_picoF) (in sodium_potassium_pump_current)
% i_NKA_alpha2 (picoA_per_picoF) (in sodium_potassium_pump_current)
% i_NKA_alpha2_jc (picoA_per_picoF) (in sodium_potassium_pump_current)
% i_NKA_alpha2_sl (picoA_per_picoF) (in sodium_potassium_pump_current)
% i_NKA_jc (picoA_per_picoF) (in sodium_potassium_pump_current)
% i_NKA_sl (picoA_per_picoF) (in sodium_potassium_pump_current)
% sigma (dimensionless) (in sodium_potassium_pump_current)
% i_K1 (picoA_per_picoF) (in time_independent_K_I)
% ass (dimensionless) (in ultra_rapidly_activating_delayed_rectifier_K_I)
% i_Kur (picoA_per_picoF) (in ultra_rapidly_activating_delayed_rectifier_K_I)
% iss (dimensionless) (in ultra_rapidly_activating_delayed_rectifier_K_I)
% tau_aur (millisecond) (in ultra_rapidly_activating_delayed_rectifier_K_I)
% tau_iur (millisecond) (in ultra_rapidly_activating_delayed_rectifier_K_I)

%-------------------------------------------------------------------------------
% Computation
%-------------------------------------------------------------------------------

% time (millisecond)

O_ClCa = 0.2/(1.0+exp(-(Y(14)-46.7)/7.8));
i_ClCa = p.g_ClCa*O_ClCa*Y(10)/(Y(10)+p.Km_Cl)*(Y(14)-p.E_Cl);
FVRT = p.F*Y(14)/(p.R*p.T);

if (abs(FVRT) > 0.00001)
   i_cl = p.p_cl*p.F^2.0*Y(14)/(p.R*p.T*p.Cm)*(Y(1)-p.Clo*exp(p.F*Y(14)/(p.R*p.T)))/(1.0-exp(p.F*Y(14)/(p.R*p.T)));
else
   i_cl = p.p_cl*p.F*0.00001/p.Cm*(Y(1)-p.Clo*exp(0.00001))/(1.0-exp(0.00001));
end;

OHi = 10.0^(-14.0+Y(27))*1000000.0;
K_OH_che = 10.0^-(14.0-p.pK_H_che)*1000000.0;
beta_f_1_che = p.k_f_1_che*p.K_Cl_che*OHi/(K_OH_che*p.K_Cl_che+p.K_Cl_che*OHi+K_OH_che*Y(1));
OHo = 10.0^(-14.0+p.pHo)*1000000.0;
beta_f_2_che = p.k_f_2_che*K_OH_che*p.Clo/(K_OH_che*p.K_Cl_che+p.K_Cl_che*OHo+K_OH_che*p.Clo);
beta_r_1_che = p.k_r_1_che*p.K_Cl_che*OHo/(K_OH_che*p.K_Cl_che+p.K_Cl_che*OHo+K_OH_che*p.Clo);
k_r_2_che = p.k_f_2_che*p.k_f_1_che/p.k_r_1_che;
beta_r_2_che = k_r_2_che*K_OH_che*Y(1)/(K_OH_che*p.K_Cl_che+p.K_Cl_che*OHi+K_OH_che*Y(1));
Jexch_che = (beta_f_1_che*beta_f_2_che-beta_r_1_che*beta_r_2_che)/(beta_f_1_che+beta_r_1_che+beta_f_2_che+beta_r_2_che);
Jche = 1000.0*Jexch_che;
Jae = 0.0;
dY(1, 1) = (i_ClCa+i_cl)*p.Atot*p.Cm/(p.Vmyo*p.F)+Jche+Jae;
FVRT_Ca = 2.0*FVRT;
expVL = exp((Y(14)-p.V_L)/p.delta_V_L);
alpha_p = expVL/(p.t_L*(expVL+1.0));
alpha_m = p.phi_L/p.t_L;
epsilon_p = (expVL+p.a)/(p.tau_L*p.K_L*(expVL+1.0));
epsilon_m = p.b*(expVL+p.a)/(p.tau_L*(p.b*expVL+p.a));
y_gate_inf = 1.0/(1.0+exp((Y(14)+p.LCC_C1)/p.LCC_C2))+p.LCC_C3/(1.0+exp((-Y(14)+p.LCC_C4)/p.LCC_C5));
y_gate_tau = p.LCC_C6+p.LCC_C7/(1.0+exp((Y(14)+p.LCC_C8)/p.LCC_C9));
Cjc = 1.0-Y(4)-Y(2);
Csl = 1.0-Y(5)-Y(3);
dY(4, 1) = alpha_p*Cjc-alpha_m*Y(4);
dY(5, 1) = alpha_p*Csl-alpha_m*Y(5);
dY(2, 1) = epsilon_p*Cjc*Y(11)-epsilon_m*Y(2);
dY(3, 1) = epsilon_p*Csl*Y(12)-epsilon_m*Y(3);
dY(6, 1) = (y_gate_inf-Y(6))/y_gate_tau;
Ajc = p.Atot*p.TT_ratio;

if (abs(FVRT_Ca) > 0.00001)
   i_CaL_jc = -p.P_CaL*2.0*p.F/(Ajc*p.Cm)*Y(4)*Y(6)*FVRT_Ca/(1.0-exp(-FVRT_Ca))*(p.Cao*exp(-FVRT_Ca)-Y(11));
else
   i_CaL_jc = -p.P_CaL*2.0*p.F/(Ajc*p.Cm)*Y(4)*Y(6)*0.00001/(1.0-exp(-0.00001))*(p.Cao*exp(-0.00001)-Y(11));
end

Asl = p.Atot*(1.0-p.TT_ratio);

if (abs(FVRT_Ca) > 0.00001)
   i_CaL_sl = -0.1*p.P_CaL*2.0*p.F/(Asl*p.Cm)*Y(5)*Y(6)*FVRT_Ca/(1.0-exp(-FVRT_Ca))*(p.Cao*exp(-FVRT_Ca)-Y(12));
else
   i_CaL_sl = -0.1*p.P_CaL*2.0*p.F/(Asl*p.Cm)*Y(5)*Y(6)*0.00001/(1.0-exp(-0.00001))*(p.Cao*exp(-0.00001)-Y(12));
end

i_CaL = i_CaL_jc*p.TT_ratio+i_CaL_sl*(1.0-p.TT_ratio);
i_anion = p.p_anion*p.F^2.0*Y(14)/(p.R*p.T*p.Cm)*(Y(7)-p.anion_o*exp(p.F*Y(14)/(p.R*p.T)))/(1.0-exp(p.F*Y(14)/(p.R*p.T)));
dY(7, 1) = i_anion*p.Atot*p.Cm/(p.Vmyo*p.F)+p.J_acid_flux_in;
E_Ca_jc = p.R*p.T/(2.0*p.F)*log(p.Cao/Y(11));
E_Ca_sl = p.R*p.T/(2.0*p.F)*log(p.Cao/Y(12));
i_Cab_jc = p.g_Cab*(Y(14)-E_Ca_jc);
i_Cab_sl = p.g_Cab*(Y(14)-E_Ca_sl);
i_Cab = i_Cab_jc*p.TT_ratio+i_Cab_sl*(1.0-p.TT_ratio);
J_Cab = -i_Cab*p.Atot*p.Cm/(2.0*p.Vmyo*p.F);
i_PMCA_sl = p.i_PMCA_max*Y(12)^2.0/(p.Km_PMCA^2.0+Y(12)^2.0);
i_NCX_sl = 0.725*p.I_NCX_max/(1.0+(p.K_mAllo/Y(12))^2.0)*(Y(40)^3.0*p.Cao*exp(p.eta*Y(14)*p.F/(p.R*p.T))-p.Nao^3.0*Y(12)*exp((p.eta-1.0)*Y(14)*p.F/(p.R*p.T)))/((p.K_mCao*Y(40)^3.0+p.K_mNao^3.0*Y(12)+p.K_mNai^3.0*p.Cao*(1.0+Y(12)/p.K_mCai)+p.K_mCai*p.Nao^3.0*(1.0+(Y(40)/p.K_mNai)^3.0)+Y(40)^3.0*p.Cao+p.Nao^3.0*Y(12))*(1.0+p.k_sat*exp((p.eta-1.0)*Y(14)*p.F/(p.R*p.T))));
JCa_sl = -(i_Cab_sl+i_PMCA_sl-2.0*i_NCX_sl+i_CaL_sl)*Asl*p.Cm/(2.0*p.Vsl*p.F);
i_PMCA_jc = p.i_PMCA_max*Y(11)^2.0/(p.Km_PMCA^2.0+Y(11)^2.0);
i_NCX_jc = 2.175*p.I_NCX_max/(1.0+(p.K_mAllo/Y(11))^2.0)*(Y(38)^3.0*p.Cao*exp(p.eta*Y(14)*p.F/(p.R*p.T))-p.Nao^3.0*Y(11)*exp((p.eta-1.0)*Y(14)*p.F/(p.R*p.T)))/((p.K_mCao*Y(38)^3.0+p.K_mNao^3.0*Y(11)+p.K_mNai^3.0*p.Cao*(1.0+Y(11)/p.K_mCai)+p.K_mCai*p.Nao^3.0*(1.0+(Y(38)/p.K_mNai)^3.0)+Y(38)^3.0*p.Cao+p.Nao^3.0*Y(11))*(1.0+p.k_sat*exp((p.eta-1.0)*Y(14)*p.F/(p.R*p.T))));
JCa_jc = -(i_Cab_jc+i_PMCA_jc-2.0*i_NCX_jc+i_CaL_jc)*Ajc*p.Cm/(2.0*p.Vjc*p.F);
Bsl = (1.0+p.Bmax*p.Kd_buffer/(p.Kd_buffer+Y(12))^2.0)^-1.0;
tau_Ca_slcyt = (p.tau_ca_slcyt_const/p.Vmyo)^-1.0;
JCa_slcyt = (Y(12)-Y(10))/tau_Ca_slcyt;
tau_Ca_jcsl = (p.tau_Ca_jcsl_const/p.Vsl)^-1.0;
JCa_jcsl = (Y(11)-Y(12))/tau_Ca_jcsl;
dY(12, 1) = Bsl*(JCa_sl-JCa_slcyt*p.Vmyo/p.Vsl+JCa_jcsl);
Bjc = (1.0+p.Bmax*p.Kd_buffer/(p.Kd_buffer+Y(11))^2.0)^-1.0;
J_rel = p.V_rel*(Y(34)+Y(35))*(Y(8)-Y(11))*Y(13);
dY(11, 1) = Bjc*(JCa_jc-JCa_jcsl*p.Vsl/p.Vjc+J_rel*p.VJSR/p.Vjc);
Bi = (1.0+p.Bmax*p.Kd_buffer/(p.Kd_buffer+Y(10))^2.0)^-1.0;
J_leak = p.V_leak*(Y(9)-Y(10));
J_serca = p.Vup*Y(10)^p.nH/(p.Km_up^p.nH+Y(10)^p.nH);
dY(10, 1) = Bi*(J_leak+JCa_slcyt-J_serca);
BJSR = (1.0+p.CSQN_tot*p.Km_CSQN/(p.Km_CSQN+Y(8))^2.0)^-1.0;
J_tr = (Y(9)-Y(8))/p.tau_tr;
dY(8, 1) = BJSR*(J_tr-J_rel);
dY(9, 1) = (J_serca-J_leak)*p.Vmyo/p.VNSR-J_tr*p.VJSR/p.VNSR;
CaSR = Y(8)*p.CSQN_tot/(Y(8)+p.Km_CSQN)*p.VJSR/p.Vmyo+Y(9)*p.VNSR/p.Vmyo;
P_open = (Y(34)+Y(35))*Y(13);
J_netup = J_serca-J_leak;
dY(13, 1) = p.A_PRyR*Y(13)+p.B_PRyR*i_CaL/p.i_CaL_max*exp(-(Y(14)-5.0)^2.0/648.0);
i_PMCA = i_PMCA_jc*p.TT_ratio+i_PMCA_sl*(1.0-p.TT_ratio);
J_PMCA_sl = -i_PMCA*Asl*p.Cm/(2.0*p.Vsl*p.F);
J_PMCA_jc = -i_PMCA*Ajc*p.Cm/(2.0*p.Vjc*p.F);
J_PMCA_total = J_PMCA_sl*p.Vsl/p.Vmyo+J_PMCA_jc*p.Vjc/p.Vmyo;
i_PMCA_proton = -i_PMCA;
J_PMCA_proton = -J_PMCA_total*2.0;
% past = floor(time/p.stim_period)*p.stim_period;
% IVClamp = 0;
% i_Stim = 0;
i_NCX = i_NCX_jc*p.TT_ratio+i_NCX_sl*(1.0-p.TT_ratio);
E_Na_jc = p.R*p.T/p.F*log((0.9*p.Nao+0.1*p.Ko)/(0.9*Y(38)+0.1*Y(28)));
i_Na_jc = p.g_Na*Y(22)*(Y(14)-E_Na_jc);
E_Na_sl = p.R*p.T/p.F*log((0.9*p.Nao+0.1*p.Ko)/(0.9*Y(40)+0.1*Y(28)));
i_Na_sl = p.g_Na*Y(22)*(Y(14)-E_Na_sl);
i_Na = i_Na_jc*p.TT_ratio+i_Na_sl*(1.0-p.TT_ratio);
i_Nab_jc = p.g_Nab*(Y(14)-E_Na_jc);
i_Nab_sl = p.g_Nab*(Y(14)-E_Na_sl);
i_Nab = i_Nab_jc*p.TT_ratio+i_Nab_sl*(1.0-p.TT_ratio);
sigma = 1.0/7.0*(exp(p.Nao/67300.0)-1.0);
f_NKA_alpha1 = 1.0/(1.0+0.2946*exp(-0.1*Y(14)*p.F/(p.R*p.T))+0.0164*sigma*exp(-Y(14)*p.F/(p.R*p.T)));
i_NKA_alpha1_jc = 1.276*p.i_NKA_max_alpha1*f_NKA_alpha1*1.0/(1.0+(p.Km_Na_alpha1/Y(38))^p.nH_NKAalpha1)*p.Ko/(p.Ko+p.Km_Ko);
f_NKA_alpha2 = 1.0/(1.0+0.1245*exp(-0.1*Y(14)*p.F/(p.R*p.T))+0.089*sigma*exp(-Y(14)*p.F/(p.R*p.T)));
i_NKA_alpha2_jc = 2.98*p.i_NKA_max_alpha2*f_NKA_alpha2*1.0/(1.0+(p.Km_Na_alpha2/Y(38))^p.nH_NKAalpha2)*p.Ko/(p.Ko+p.Km_Ko);
i_NKA_jc = i_NKA_alpha1_jc+i_NKA_alpha2_jc;
i_NKA_alpha1_sl = 0.931*p.i_NKA_max_alpha1*f_NKA_alpha1*1.0/(1.0+(p.Km_Na_alpha1/Y(40))^p.nH_NKAalpha1)*p.Ko/(p.Ko+p.Km_Ko);
i_NKA_alpha2_sl = 0.535*p.i_NKA_max_alpha2*f_NKA_alpha2*1.0/(1.0+(p.Km_Na_alpha2/Y(40))^p.nH_NKAalpha2)*p.Ko/(p.Ko+p.Km_Ko);
i_NKA_sl = i_NKA_alpha1_sl+i_NKA_alpha2_sl;
i_NKA = i_NKA_jc*p.TT_ratio+i_NKA_sl*(1.0-p.TT_ratio);
E_K = p.R*p.T/p.F*log(p.Ko/Y(28));
i_Kto_f = p.g_Kto_f*Y(23)^3.0*Y(24)*(Y(14)-E_K);
i_K1 = p.g_K1*p.Ko/(p.Ko+210.0)*(Y(14)-E_K)/(1.0+exp(0.0896*(Y(14)-E_K)));
i_Ks = p.g_Ks*Y(36)^2.0*(Y(14)-E_K);

% Old i_Kur
% i_Kur = p.g_Kur*Y(42)*Y(43)*(Y(14)-E_K);

% New i_Kur
i_Kur1 = p.g_Kur1*Y(42)*Y(44)*(Y(14)-E_K);
i_Kur2 = p.g_Kur2*Y(43)*Y(45)*(Y(14)-E_K);
i_Kur = i_Kur1 + i_Kur2;

i_Kss = p.g_Kss*Y(25)*Y(26)*(Y(14)-E_K);
i_Kr = p.g_Kr*Y(32)*(Y(14)-p.R*p.T/p.F*log((0.98*p.Ko+0.02*p.Nao)/(0.98*Y(28)+0.02*Y(37))));

%-------------------- Stimulating current -----------------------
%-------------------- Caffeine ----------------------------------
if isfield(p,'caffeine') && p.caffeine == 1
    dY(13,1) = 0.95*(1-Y(13));
else
    dY(13, 1) = p.A_PRyR*Y(13)+p.B_PRyR*i_CaL/p.i_CaL_max*exp(-(Y(14)-5.0)^2.0/648.0);
    if (time >= p.pacestart && time < p.pacestart+p.pacedur)
    i_Stim  = p.paceamp;
    else
    i_Stim  = 0;
    end
end
%end----------------- Caffeine ----------------------------------
%-------------------- Voltage Clamp current ---------------------
IVClamp = 0;
if ~isempty(p.VClampAmp)
    i_Stim = 0;
    ind = find(p.VClampTimes>time);
    if ~isempty(ind)
    IVClamp = 1e3*(Y(14)-p.VClampAmp(ind(1)))/(p.VClampR*p.Atot);
    end
else
    IVClamp = 0;
    if (time >= p.pacestart && time < p.pacestart+p.pacedur)
    i_Stim  = p.paceamp;
    else
    i_Stim  = 0;
    end
end
%end----------------- Voltage Clamp current ---------------------
%end----------------- Stimulating current -----------------------

if isfield(p,'VClampAmp')
   dY(14, 1) = -(i_cl+i_CaL+i_NCX+i_Cab+i_Na+i_Nab+i_NKA+i_Kto_f+i_K1+i_Ks+i_Kur+i_Kss+i_Kr+i_ClCa+i_Stim+i_anion+IVClamp);
elseif (time < p.stim_period*p.prepulses_number) 
   dY(14, 1) = -(i_cl+i_CaL+i_NCX+i_Cab+i_Na+i_Nab+i_NKA+i_Kto_f+i_K1+i_Ks+i_Kur+i_Kss+i_Kr+i_ClCa+i_Stim+i_anion+IVClamp);
else
   dY(14, 1) = 0.0;
end;

C_Na3 = 1.0-(Y(22)+Y(15)+Y(16)+Y(21)+Y(17)+Y(18)+Y(19)+Y(20));
alpha_Na11 = 3.802/(0.1027*exp(-(Y(14)+2.5)/17.0)+0.2*exp(-(Y(14)+2.5)/150.0));
beta_Na12 = 0.2*exp(-(Y(14)-2.5)/20.3);
alpha_Na3 = 7.0e-7*exp(-(Y(14)+7.0)/7.7);
beta_Na11 = 0.1917*exp(-(Y(14)+2.5)/20.3);
alpha_Na12 = 3.802/(0.1027*exp(-(Y(14)+2.5)/15.0)+0.23*exp(-(Y(14)+2.5)/150.0));
beta_Na3 = 0.0084+0.00002*(Y(14)+7.0);
dY(16, 1) = alpha_Na11*C_Na3+beta_Na12*Y(15)+alpha_Na3*Y(19)-(beta_Na11*Y(16)+alpha_Na12*Y(16)+beta_Na3*Y(16));
beta_Na13 = 0.22*exp(-(Y(14)-7.5)/20.3);
alpha_Na13 = 3.802/(0.1027*exp(-(Y(14)+2.5)/12.0)+0.25*exp(-(Y(14)+2.5)/150.0));
dY(15, 1) = alpha_Na12*Y(16)+beta_Na13*Y(22)+alpha_Na3*Y(21)-(beta_Na12*Y(15)+alpha_Na13*Y(15)+beta_Na3*Y(15));
alpha_Na2 = 1.0/(0.188495*exp(-(Y(14)+7.0)/16.6)+0.393956);
beta_Na2 = alpha_Na13*alpha_Na2*alpha_Na3/(beta_Na13*beta_Na3);
dY(22, 1) = alpha_Na13*Y(15)+beta_Na2*Y(21)-(beta_Na13*Y(22)+alpha_Na2*Y(22));
beta_Na4 = alpha_Na3;
alpha_Na4 = alpha_Na2/1000.0;
dY(21, 1) = alpha_Na2*Y(22)+beta_Na3*Y(15)+beta_Na4*Y(17)+alpha_Na12*Y(19)-(beta_Na2*Y(21)+alpha_Na3*Y(21)+alpha_Na4*Y(21)+beta_Na12*Y(21));
beta_Na5 = alpha_Na3/50.0;
alpha_Na5 = alpha_Na2/95000.0;
dY(17, 1) = alpha_Na4*Y(21)+beta_Na5*Y(18)-(beta_Na4*Y(17)+alpha_Na5*Y(17));
dY(18, 1) = alpha_Na5*Y(17)-beta_Na5*Y(18);
dY(19, 1) = alpha_Na11*Y(20)+beta_Na12*Y(21)+beta_Na3*Y(16)-(beta_Na11*Y(19)+alpha_Na12*Y(19)+alpha_Na3*Y(19));
dY(20, 1) = beta_Na11*Y(19)+beta_Na3*C_Na3-(alpha_Na11*Y(20)+alpha_Na3*Y(20));
itof_iss = 1.0/(1.0+exp((Y(14)+51.4)/5.0));
tau_ito_f = 9.6645+10.9362/(1.0+exp((Y(14)+51.4)/5.0));
alpha_a = 0.18064*exp(0.03577*(Y(14)+p.alpha_a_const));
beta_a = 0.3956*exp(-0.06237*(Y(14)+p.alpha_b_const));
dY(23, 1) = alpha_a*(1.0-Y(23))-beta_a*Y(23);
dY(24, 1) = (itof_iss-Y(24))/tau_ito_f;

%original
%ass = 1.0/(1.0+exp(-(Y(14)+15.0)/20.0));

%new
ass = 1/(1+exp(-(Y(14)+15)/14));

tau_Kss = 39.3*exp(-0.05*Y(14))+13.17;
dY(25, 1) = (ass-Y(25))/tau_Kss;
dY(26, 1) = 0.0;
Ho = 10.0^-p.pHo*1000000.0;
JCO2 = 0.0;
Jnbc = 0.0;
Hi = 10.0^-Y(27)*1000000.0;
alpha_r_1_nhe = Y(37)/p.K_Na_nhe*p.k_r_1_nhe/((1.0+Y(37)/p.K_Na_nhe)*(1.0+Hi/p.K_H_nhe));
alpha_f_2_nhe = Hi/p.K_H_nhe*p.k_f_2_nhe/((1.0+Y(37)/p.K_Na_nhe)*(1.0+Hi/p.K_H_nhe));
Jnhe = 1000.0*Hi^p.nH_nhe/(Hi^p.nH_nhe+p.K_i_nhe^p.nH_nhe)*(p.alpha_f_1_nhe*alpha_f_2_nhe-alpha_r_1_nhe*p.alpha_r_2_nhe)/(p.alpha_f_1_nhe+alpha_f_2_nhe+alpha_r_1_nhe+p.alpha_r_2_nhe);
beta_i_1 = log(10.0)*(10.0^-Y(27)*1000.0+10.0^(p.pK1-Y(27))*p.B1/(1.0+10.0^(p.pK1-Y(27)))^2.0+10.0^(p.pK2-Y(27))*p.B2/(1.0+10.0^(p.pK2-Y(27)))^2.0);
Jhyd = 0.0;
dY(27, 1) = -1.0/beta_i_1*(Jche-Jnhe+p.J_acid_flux_in+Jhyd+J_PMCA_proton);
dY(28, 1) = -(i_Stim+i_Kto_f+i_K1+i_Ks+i_Kss+i_Kur+i_Kr-2.0*i_NKA)*p.Atot*p.Cm/(p.Vmyo*p.F);
C_K0 = 1.0-(Y(29)+Y(30)+Y(32)+Y(31));
beta_a1 = 0.0000689*exp(-0.04178*Y(14));
alpha_a1 = 0.0335*exp(0.0109*Y(14));
dY(30, 1) = p.kf*Y(29)+beta_a1*Y(32)-(p.kb*Y(30)+alpha_a1*Y(30));
alpha_a0 = 0.022348*exp(0.01176*Y(14));
beta_a0 = 0.047002*exp(-0.0631*Y(14));
dY(29, 1) = alpha_a0*C_K0+p.kb*Y(30)-(beta_a0*Y(29)+p.kf*Y(29));
beta_i_2 = 0.006497*exp(-0.03268*(Y(14)+5.0));
alpha_i = 0.0703*exp(0.0287*(Y(14)+5.0));
dY(32, 1) = alpha_a1*Y(30)+beta_i_2*Y(31)-(beta_a1*Y(32)+alpha_i*Y(32));
dY(31, 1) = alpha_i*Y(32)-beta_i_2*Y(31);
P_C1 = 1.0-(Y(33)+Y(34)+Y(35));
dY(34, 1) = p.k_plus_a*Y(11)^p.n*P_C1+p.k_minus_b*Y(35)+p.k_minus_c*Y(33)-(p.k_minus_a*Y(34)+p.k_plus_b*Y(11)^p.m*Y(34)+p.k_plus_c*Y(34));
dY(35, 1) = p.k_plus_b*Y(11)^p.m*Y(34)-p.k_minus_b*Y(35);
dY(33, 1) = p.k_plus_c*Y(34)-p.k_minus_c*Y(33);
alpha_n = 0.00000481333*(Y(14)+26.5)/(1.0-exp(-0.128*(Y(14)+26.5)));
beta_n = 0.0000953333*exp(-0.038*(Y(14)+26.5));
dY(36, 1) = alpha_n*(1.0-Y(36))-beta_n*Y(36);
J_NCX_sl = i_NCX_sl*Asl*p.Cm/(p.Vsl*p.F);
J_NCX_jc = i_NCX_jc*Ajc*p.Cm/(p.Vjc*p.F);
J_NCX_total = J_NCX_sl*p.Vsl/p.Vmyo+J_NCX_jc*p.Vjc/p.Vmyo;
dNa_jct_buf = p.kon*Y(38)*(p.Bmax_jct-Y(39))-p.koff*Y(39);
dNa_SL_buf = p.kon*Y(40)*(p.Bmax_SL-Y(41))-p.koff*Y(41);
dY(39, 1) = dNa_jct_buf;
dY(41, 1) = dNa_SL_buf;
tau_Na_jcsl = (p.tau_Na_jcsl_const/p.Vsl)^-1.0;
tau_Na_slcyt = (p.tau_Na_slcyt_const/p.Vmyo)^-1.0;
JNa_slcyt = (Y(40)-Y(37))/tau_Na_slcyt;
JNa_jcsl = (Y(38)-Y(40))/tau_Na_jcsl;
dY(37, 1) = JNa_slcyt+Jnhe+Jnbc;
dY(40, 1) = -(i_Na_sl+i_Nab_sl+3.0*i_NKA_sl+3.0*i_NCX_sl)*Asl*p.Cm/(p.Vsl*p.F)-JNa_slcyt*p.Vmyo/p.Vsl+JNa_jcsl-dNa_SL_buf;
dY(38, 1) = -(i_Na_jc+i_Nab_jc+3.0*i_NKA_jc+3.0*i_NCX_jc)*Ajc*p.Cm/(p.Vjc*p.F)-JNa_jcsl*p.Vsl/p.Vjc-dNa_jct_buf;
i_NKA_alpha1 = i_NKA_alpha1_jc*p.TT_ratio+i_NKA_alpha1_sl*(1.0-p.TT_ratio);
i_NKA_alpha2 = i_NKA_alpha2_jc*p.TT_ratio+i_NKA_alpha2_sl*(1.0-p.TT_ratio);
J_Na_NKA = i_NKA*3.0*p.Atot*p.Cm/(p.F*p.Vmyo);

%I_Kur and I_Kss state variables

%original
% iss = 1.0/(1.0+exp((Y(14)+42.1)/5.4));
% tau_aur = 0.493*exp(-0.0629*Y(14))+p.tau_a_const;
% dY(42, 1) = (ass-Y(42))/tau_aur;
% tau_iur = p.tau_i_const+1000.0/(1.0+exp((Y(14)+42.1)/5.4));
% dY(43, 1) = (iss-Y(43))/tau_iur;

%new
iss = 1/(1+exp((Y(14)+48)/6.2));
tau_aur1 = 0.95+0.05*exp(-0.08*Y(14));
dY(42) = (ass-Y(42))/tau_aur1;
tau_aur2 = 1+25/(1+exp(-(Y(14)+45)/8))+20*exp(-((Y(14)+35)/25)^2);
dY(43) = (ass-Y(43))/tau_aur2;
tau_iur1 = 400+900*exp(-((Y(14)+55)/16)^2)-250/(1+exp(-(Y(14)+60)/8));
dY(44) = (iss-Y(44))/tau_iur1;
tau_iur2 = 400+900*exp(-((Y(14)+55)/16)^2)+550/(1+exp(-(Y(14)+60)/8));
dY(45) = (iss-Y(45))/tau_iur2;

nout = max(nargout,1)-1; % Number of output parameters in excess of dx
if nout > 0
  if nargin == 4 && length(loginfo)> 0
    logvar = fieldnames(loginfo);
    nlog = length(logvar);
    argout = cell(1,nlog);
    for i=1:nlog
      switch logvar{i}
       case 'state_der'
        argout{i} = dx;
       case 'currents'
        argout{i} = [i_ClCa,i_CaL_jc,i_CaL_sl,i_CaL,i_anion, ...
                 i_Cab_jc,i_Cab_sl,i_Cab,i_PMCA_sl,i_NCX_sl,i_PMCA_jc,...
                 i_NCX_jc,i_PMCA,i_PMCA_proton,i_Stim,i_NCX,i_Na_jc,i_Na_sl,...
                 i_Na,i_Nab_jc,i_Nab_sl,i_Nab,i_NKA_alpha1_jc,i_NKA_alpha2_jc,...
                 i_NKA_jc,i_NKA_alpha1_sl,i_NKA_alpha2_sl,i_NKA_sl,i_NKA,i_Kto_f,...
                 i_K1,i_Ks,i_Kur,i_Kur1,i_Kur2,i_Kss,i_Kr,i_NKA_alpha1,i_NKA_alpha2,IVClamp];
       case 'ionfluxes'
        argout{i} = [Jche,Jae,J_Cab,JCa_sl,JCa_jc,JCa_slcyt,...
                 JCa_jcsl,J_rel,J_leak,J_serca,J_tr,J_netup,...
                 J_PMCA_sl,J_PMCA_jc,J_PMCA_total,J_PMCA_proton,JCO2,...
                 Jnbc,Jnhe,Jhyd,J_NCX_sl,J_NCX_jc,J_NCX_total,...
                 JNa_slcyt,JNa_jcsl,J_Na_NKA];
       case 'buffers'
        argout{i} = [Bsl,Bjc,Bi,BJSR];
       case 'rates'
        argout{i} = [];
      case 'hiddenstates'
       argout{i} = [];
      end
    end
    varargout{1} = argout;
  else
    varargout{1} = [];
  end
end


%===============================================================================
% End of file
%========================