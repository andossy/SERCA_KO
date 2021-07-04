clear all
close all

%FF
%Control
p.g_Kur2 = 0.14;
p.g_Kss = 0.13;
p.g_Kto_f = 0.35;
p.g_K1 = 0.30*3;
p.stim_amplitude = -150;

%beta-adrenergic regulation
p.P_CaL = 1.6e-07*3;
p.g_Kur1 = 0.15*1.15;  
p.V_L = -5-9;
p.Km_up = 0.4928*0.4;
p.Km_Na_alpha1 = 21000.0*0.72;
p.Km_Na_alpha2 = 21000.0*0.72;

%KO
%Control

p2.LCC_C1 = 38;
p2.LCC_C3 = 0.4;
p2.phi_L = 1.5;
p2.LCC_C6 = 40;
p2.g_Kur2 = p.g_Kur2*0.75;
p2.g_Kss = p.g_Kss*0.75;
p2.g_Kto_f = 0.30;
p2.g_K1 = p.g_K1;
p2.stim_amplitude = p.stim_amplitude;

%beta-adrenergic regulation
p2.P_CaL = 0.6e-07*3;
p2.g_Kur1 = 0.15*1.15;
p2.Km_up = p.Km_up;
p2.V_L = p.V_L;
p2.Km_Na_alpha1 = p.Km_Na_alpha1;
p2.Km_Na_alpha2 = p.Km_Na_alpha2;

[T1,S1,loginfo1,x1,ioninfo1] = pacemodel(@FF,1000,60,p,'ioninfo');

for i = 2:length(T1)
    subplot(3,1,1), plot([T1{i};2000],[S1{i}(:,14);S1{i}(end,14)]);
    hold on
    subplot(3,1,2),plot([T1{i};2000],[loginfo1.currents{i}(:,19);loginfo1.currents{i}(end,19)],'k')
    hold on
    subplot(3,1,2),plot([T1{i};2000],[loginfo1.currents{i}(:,4);loginfo1.currents{i}(end,4)],'y')
    subplot(3,1,2),plot([T1{i};2000],[loginfo1.currents{i}(:,16);loginfo1.currents{i}(end,16)],'g')
%     subplot(3,1,3),plot([T1{i};2000],[loginfo1.currents{i}(:,30);loginfo1.currents{i}(end,30)],'k')
%     hold on
%     subplot(3,1,3),plot([T1{i};2000],[loginfo1.currents{i}(:,33);loginfo1.currents{i}(end,33)],'y')
%     subplot(3,1,3),plot([T1{i};2000],[loginfo1.currents{i}(:,34);loginfo1.currents{i}(end,34)],'g')  
end

for i = 2:length(T1)
    if i > 2
        tim=cat(1,tim,T1{i}+tim(end));
        mV=cat(1,mV,S1{i}(:,14));
        INa=cat(1,INa,loginfo1.currents{i}(:,19));
        ICa=cat(1,ICa,loginfo1.currents{i}(:,4));
    else
        tim=T1{i};
        mV=S1{i}(:,14);
        Ito=loginfo1.currents{i}(:,30);
        Ikur=loginfo1.currents{i}(:,33);
        Ikss=loginfo1.currents{i}(:,34);
        ICa=loginfo1.currents{i}(:,4);
        INa=loginfo1.currents{i}(:,19);
    end
end

figure 

[T3,S3,loginfo3,x3,ioninfo3] = pacemodel(@KO,1000,60,p2,'ioninfo');

for i = 1:length(T3)
    if i > 1
        tim3=cat(1,tim3,T3{i}+tim3(end));
        mV3=cat(1,mV3,S3{i}(:,14));
        Ito3=cat(1,Ito3,S3{i}(:,30));
        Ikur3=cat(1,Ikur3,S3{i}(:,33));
        Ikss3=cat(1,Ikss3,S3{i}(:,34));
        INa3=cat(1,INa3,loginfo3.currents{i}(:,19));
        ICa3=cat(1,ICa3,loginfo3.currents{i}(:,4));
    else
        tim3=T3{i};
        mV3=S3{i}(:,14);
        INa3=loginfo3.currents{i}(:,19);
        Ito3=loginfo3.currents{i}(:,30);
        Ikur3=loginfo3.currents{i}(:,33);
        Ikss3=loginfo3.currents{i}(:,34);
        ICa3=loginfo3.currents{i}(:,4);   
    end
end
T_FF = T1{1};
V_FF = S1{1}(:,14);
T_KO = T3{1};
V_KO = S3{1}(:,14);
for i = 2:length(T3)
    subplot(3,1,1), plot(T3{i},S3{i}(:,14),'r');
    hold on
    subplot(3,1,2),plot(T3{i},loginfo3.currents{i}(:,19),'m')
    hold on
    subplot(3,1,2),plot(T3{i},loginfo3.currents{i}(:,4),'b')
    subplot(3,1,2),plot(T3{i},loginfo3.currents{i}(:,16),'r') 
    subplot(3,1,3),plot(T3{i},loginfo3.currents{i}(:,30),'m')
    hold on
    subplot(3,1,3),plot(T3{i},loginfo3.currents{i}(:,33),'b')
    subplot(3,1,3),plot(T3{i},loginfo3.currents{i}(:,34),'r') 
    subplot(3,1,3),plot(T3{i},loginfo3.currents{i}(:,35),'b*')
    
    T_FF = cat(1,T_FF,T_FF(end)+T1{i});
    V_FF = cat(1,V_FF,S1{i}(:,14));
    T_KO = cat(1,T_KO,T_KO(end)+T3{i});
    V_KO = cat(1,V_KO,S3{i}(:,14));
end

figure
plot(T_FF,V_FF);

figure
plot(T_KO,V_KO);


