function inits = initmaker(xvals)

Y = xvals;
YNames = {'Cli', 'Ijc', 'Isl', 'Ojc', 'Osl', 'y_gate', 'anion_i', 'CaJSR', 'CaNSR', 'Cai', 'Cajc', 'Casl', 'P_RyR', 'V', 'C_Na1', 'C_Na2', 'I1_Na', 'I2_Na', 'IC_Na2', 'IC_Na3', 'IF_Na', 'O_Na', 'ato_f', 'ito_f', 'aKss', 'iKss', 'pH_i', 'Ki', 'C_K1', 'C_K2', 'I_K', 'O_K', 'P_C2', 'P_O1', 'P_O2', 'nKs', 'Nai', 'Najc', 'Najc_buf', 'Nasl', 'Nasl_buf', 'aur1','aur2','iur1','iur2'};


for i = 1:length(Y)
    inits{i} = strcat('x0.',YNames{i},' = ',num2str(Y(i)));
end