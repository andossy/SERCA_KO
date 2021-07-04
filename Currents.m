clear all
close all

%FF
p.g_Kur1 = 0.15;
p.g_Kur2 = 0.14;
p.g_Kss = 0.13;
p.g_Kto_f = 0.35;
p.P_CaL = 1.2e-7;
p.g_K1 = 0.30*3;

%KO
p2.P_CaL = 0.6e-07; % original 2.0e-7
p2.LCC_C1 = 38;
p2.LCC_C3 = 0.4;
p2.phi_L = 1.5;
p2.LCC_C6 = 40;
p2.g_Kur1 = 0.15;
p2.g_Kur2 = p.g_Kur2*0.75;
p2.g_Kss = p.g_Kss*0.75;
p2.g_Kto_f = 0.30;
p2.g_K1 = p.g_K1;

%ICaL
Voutact = [-40,-30,-20,-10,0,10,20,30,40,50]';
for j = 1: length(Voutact)
    p.VClampAmp = [-70,Voutact(j),-70];
    p.VClampTimes = [100,200,1000];
    p2.VClampAmp = [-70,Voutact(j),-70];
    p2.VClampTimes = [100,200,1000];
    [T_CaL_FF,S_CaL_FF,loginfo_CaL_FF,x1_CaL_FF,ioninfo_CaL_FF] = pacemodel(@FF,1000,1,p,'ioninfo');
    [T_CaL_KO,S_CaL_KO,loginfo_CaL_KO,x1_CaL_KO,ioninfo_CaL_KO] = pacemodel(@KO,1000,1,p2,'ioninfo');
    startind_FF = find(T_CaL_FF{1}<100,1,'last');
    startind_KO = find(T_CaL_KO{1}<100,1,'last');
    endind_FF = find(T_CaL_FF{1}<200,1,'last');
    endind_KO = find(T_CaL_KO{1}<200,1,'last');
    [peak_ICaL_FF(j),peak_ICaL_ind_FF(j)] = min(loginfo_CaL_FF.currents{1}(startind_FF:endind_FF,4)); 
    [peak_ICaL_KO(j),peak_ICaL_ind_KO(j)] = min(loginfo_CaL_KO.currents{1}(startind_KO:endind_KO,4)); 
    peak_ICaL_FF(j)=peak_ICaL_FF(j)-loginfo_CaL_FF.currents{1}(endind_FF-1,4);
    peak_ICaL_KO(j)=peak_ICaL_KO(j)-loginfo_CaL_KO.currents{1}(endind_KO-1,4);
    if j == 5
        subplot(2,1,1), plot(T_CaL_FF{1},S_CaL_FF{1}(:,14));
        hold on
        subplot(2,1,1), plot(T_CaL_KO{1},S_CaL_KO{1}(:,14));
        subplot(2,1,2), plot(T_CaL_FF{1},loginfo_CaL_FF.currents{1}(:,4),'b'); 
        hold on
        subplot(2,1,2), plot(T_CaL_KO{1},loginfo_CaL_KO.currents{1}(:,4),'r');
    end
end
figure
plot(Voutact,peak_ICaL_FF,'b')
hold on
plot(Voutact,peak_ICaL_KO,'r')

%IK short
figure;
Voutact = [-60,-50,-40,-30,-20,-10,0,10,20,30,40];
for j = 1: length(Voutact)
    p.VClampAmp = [-70,Voutact(j),-70];
    p.VClampTimes = [100,400,1000];
    p2.VClampAmp = [-70,Voutact(j),-70];
    p2.VClampTimes = [100,400,1000];
    [T_k_FF,S_k_FF,loginfo_k_FF,x1_k_FF,ioninfo_k_FF] = pacemodel(@FF,1000,1,p,'ioninfo');
    [T_k_KO,S_k_KO,loginfo_k_KO,x1_k_KO,ioninfo_k_KO] = pacemodel(@KO,1000,1,p2,'ioninfo');
    bulk_k_FF = sum(loginfo_k_FF.currents{1}(:,30:33),2)+loginfo_k_FF.currents{1}(:,36)+loginfo_k_FF.currents{1}(:,37);
    bulk_k_KO = sum(loginfo_k_KO.currents{1}(:,30:33),2)+loginfo_k_KO.currents{1}(:,36)+loginfo_k_KO.currents{1}(:,37);
    
    [peak_k_FF(j), t_peak_FF(j)] = max(bulk_k_FF);
    [peak_k_KO(j),t_peak_KO(j)] = max(bulk_k_KO);
    t_end_FF(j) = find(T_k_FF{1}>400,1);
    to_end_FF(j) = bulk_k_FF(t_end_FF(j));
    t_end_KO(j) = find(T_k_KO{1}>400,1);
    to_end_KO(j) = bulk_k_KO(t_end_KO(j));
    subplot(2,1,1), plot(T_k_FF{1},S_k_FF{1}(:,14),'b')
    hold on
    plot(T_k_KO{1},S_k_KO{1}(:,14),'r')
    subplot(2,1,2), plot(T_k_FF{1},bulk_k_FF,'b'); 
    hold on
    subplot(2,1,2), plot(T_k_KO{1},bulk_k_KO,'r');
    plot(T_k_FF{1}(t_peak_FF(j)),peak_k_FF(j),'b*');
    plot(T_k_KO{1}(t_peak_KO(j)),peak_k_KO(j),'r*');
    plot(T_k_FF{1}(t_end_FF(j)),to_end_FF(j),'b*');
    plot(T_k_KO{1}(t_end_KO(j)),to_end_KO(j),'r*');
    
    Ito_FF(j) = peak_k_FF(j)-to_end_FF(j);
    Ito_KO(j) = peak_k_KO(j)-to_end_KO(j);
end
    
figure 
plot(Voutact,Ito_FF,'b');
hold on
plot(Voutact,Ito_KO,'r');

%IK
figure;
Voutact = [-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60];
for j = 1: length(Voutact)
    p.VClampAmp = [-70,Voutact(j),-70];
    p.VClampTimes = [100,600,1000];
    p2.VClampAmp = [-70,Voutact(j),-70];
    p2.VClampTimes = [100,600,1000];
    [T_k_FF,S_k_FF,loginfo_k_FF,x1_k_FF,ioninfo_k_FF] = pacemodel(@FF,1000,1,p,'ioninfo');
    [T_k_KO,S_k_KO,loginfo_k_KO,x1_k_KO,ioninfo_k_KO] = pacemodel(@KO,1000,1,p2,'ioninfo');
    bulk_k_FF = sum(loginfo_k_FF.currents{1}(:,30:33),2)+loginfo_k_FF.currents{1}(:,36)+loginfo_k_FF.currents{1}(:,37);
    bulk_k_KO = sum(loginfo_k_KO.currents{1}(:,30:33),2)+loginfo_k_KO.currents{1}(:,36)+loginfo_k_KO.currents{1}(:,37);
    
    t_SS_FF(j) = find(T_k_FF{1}>600,1)-1;
    SS_k_FF(j) = mean(bulk_k_FF(find(T_k_FF{1}>580,1):find(T_k_FF{1}>600,1)-1));
    t_SS_KO(j) = find(T_k_KO{1}>600,1)-1;
    SS_k_KO(j) = mean(bulk_k_KO(find(T_k_KO{1}>580,1):find(T_k_KO{1}>600,1)-1));
    subplot(2,1,1), plot(T_k_FF{1},S_k_FF{1}(:,14),'b')
    hold on
    plot(T_k_KO{1},S_k_KO{1}(:,14),'r')
    subplot(2,1,2), plot(T_k_FF{1},bulk_k_FF,'b'); 
    hold on
    subplot(2,1,2), plot(T_k_KO{1},bulk_k_KO,'r');
    plot(T_k_FF{1}(t_SS_FF(j)),SS_k_FF(j),'b*');
    plot(T_k_KO{1}(t_SS_KO(j)),SS_k_KO(j),'r*');
    
    IKp_FF(j) = SS_k_FF(j);
    IKp_KO(j) = SS_k_KO(j);
end

figure 
plot(Voutact,IKp_FF,'b');
hold on
plot(Voutact,IKp_KO,'r');



    