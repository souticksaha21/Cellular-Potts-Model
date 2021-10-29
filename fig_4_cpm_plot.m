clear all

dd1 = '2018_04_02_tau_100_sec_alpha_2_lambda_10_r_1_by_tau_100pts_N_5000_z_100_c_codes_cpm6s/dat/';

param = load([dd1 'cpm5g.param.dat']);
g1(1)=0.01
g1(2)=1
g1(3)=5
g1(4)=50


    for j = 1:4
        v = load([dd1 'cpm5g_g' num2str(j-1) '.v.dat']);
        CI = load([dd1 'cpm5g_g' num2str(j-1) '.CI.dat']);
        CR = load([dd1 'cpm5g_g' num2str(j-1) '.CR.dat']);

        v = v*60*60; % um/h

        Z = length(v);
        vbar1(j) = mean(v);
        dv1(j) = std(v)/sqrt(Z);
        verr1(j)=log10(vbar1(j)+dv1(j))-log10(vbar1(j));
        CIbar1(j) = mean(CI);
        dCI1(j) = std(CI)/sqrt(Z);
        CRbar1(j) = mean(CR);
        dCR1(j) = std(CR)/sqrt(Z);
    end


dd2 = '2018_04_02_tau_100_sec_alpha_2_lambda_10_r_1_by_tau_100pts_N_10000_z_100_c_codes_cpm6s/dat/';

param = load([dd2 'cpm5g.param.dat']);
g1(1)=0.01
g1(2)=1
g1(3)=5
g1(4)=50


    for j = 1:4
        v = load([dd2 'cpm5g_g' num2str(j-1) '.v.dat']);
        CI = load([dd2 'cpm5g_g' num2str(j-1) '.CI.dat']);
        CR = load([dd2 'cpm5g_g' num2str(j-1) '.CR.dat']);

        v = v*60*60; % um/h

        Z = length(v);
        vbar2(j) = mean(v);
        dv2(j) = std(v)/sqrt(Z);
        verr2(j)=log10(vbar2(j)+dv2(j))-log10(vbar2(j));
        CIbar2(j) = mean(CI);
        dCI2(j) = std(CI)/sqrt(Z);
        CRbar2(j) = mean(CR);
        dCR2(j) = std(CR)/sqrt(Z);
    end

dd3 = '2018_04_02_tau_100_sec_alpha_2_lambda_10_r_1_by_tau_100pts_N_15000_z_100_c_codes_cpm6s/dat/';

param = load([dd3 'cpm5g.param.dat']);
g1(1)=0.01
g1(2)=1
g1(3)=5
g1(4)=50


    for j = 1:4
        v = load([dd3 'cpm5g_g' num2str(j-1) '.v.dat']);
        CI = load([dd3 'cpm5g_g' num2str(j-1) '.CI.dat']);
        CR = load([dd3 'cpm5g_g' num2str(j-1) '.CR.dat']);

        v = v*60*60; % um/h

        Z = length(v);
        vbar3(j) = mean(v);
        dv3(j) = std(v)/sqrt(Z);
        verr3(j)=log10(vbar3(j)+dv3(j))-log10(vbar3(j));
        CIbar3(j) = mean(CI);
        dCI3(j) = std(CI)/sqrt(Z);
        CRbar3(j) = mean(CR);
        dCR3(j) = std(CR)/sqrt(Z);
    end
    
dd4 = '2018_04_02_tau_100_sec_alpha_2_lambda_10_r_pt33_by_tau_100pts_N_10000_z_100_c_codes_cpm6s/dat/';

param = load([dd4 'cpm5g.param.dat']);
g1(1)=0.01
g1(2)=1
g1(3)=5
g1(4)=50


    for j = 1:4
        v = load([dd4 'cpm5g_g' num2str(j-1) '.v.dat']);
        CI = load([dd4 'cpm5g_g' num2str(j-1) '.CI.dat']);
        CR = load([dd4 'cpm5g_g' num2str(j-1) '.CR.dat']);

        v = v*60*60; % um/h

        Z = length(v);
        vbar4(j) = mean(v);
        dv4(j) = std(v)/sqrt(Z);
        verr4(j)=log10(vbar4(j)+dv4(j))-log10(vbar4(j));
        CIbar4(j) = mean(CI);
        dCI4(j) = std(CI)/sqrt(Z);
        CRbar4(j) = mean(CR);
        dCR4(j) = std(CR)/sqrt(Z);
    end

    
  dd5 = '2018_04_02_tau_100_sec_alpha_2_lambda_10_r_pt1_by_tau_100pts_N_10000_z_100_c_codes_cpm6s/dat/';

param = load([dd5 'cpm5g.param.dat']);
g1(1)=0.01
g1(2)=1
g1(3)=5
g1(4)=50


    for j = 1:4
        v = load([dd5 'cpm5g_g' num2str(j-1) '.v.dat']);
        CI = load([dd5 'cpm5g_g' num2str(j-1) '.CI.dat']);
        CR = load([dd5 'cpm5g_g' num2str(j-1) '.CR.dat']);

        v = v*60*60; % um/h

        Z = length(v);
        vbar5(j) = mean(v);
        dv5(j) = std(v)/sqrt(Z);
        verr5(j)=log10(vbar5(j)+dv5(j))-log10(vbar5(j));
        CIbar5(j) = mean(CI);
        dCI5(j) = std(CI)/sqrt(Z);
        CRbar5(j) = mean(CR);
        dCR5(j) = std(CR)/sqrt(Z);
    end
    
 hold off
figure(1); 
%plot(log10(g1),log10(vbar1),'-o',log10(g1),log10(vbar2),'-x',log10(g1),log10(vbar3),'-+')
h1 = errorbar(g1,vbar2,dv2,'r-o'); %hold on;
set([h1],'LineWidth',1.2)
set(gca,'xscale','log') 
ylim ([0 40])
xlabel('g(nM/mm)')
ylabel('v(\mum/hr)')
title('v(\mum/hr)')

figure(2); 
%plot(log10(g1),CIbar1,'-o',log10(g1),CIbar2,'-x',log10(g1),CIbar3,'-+')

%h2=errorbar(g1,CIbar1,dCI1,'b-^'); hold on;
h3=errorbar(g1,CIbar2,dCI2,'r-o'); hold on;
%h4=errorbar(g1,CIbar3,dCI3,'g-s');
set([h3],'LineWidth',1.2)
%legend('N = 5000','N = 10000','N = 15000','Location','northwest')
set(gca,'xscale','log')
%set(gca,'xticklabel',arrayfun(@(x) num2str(x),log10(get(gca,'xtick')),'un',0))
xlabel('g(nM/mm)')
ylabel('CI')
title('CI')
hold off

figure(3); 
%plot(log10(g1),CRbar1,'-o',log10(g1),CRbar2,'-x',log10(g1),CRbar3,'-+')

h5=errorbar(g1,CRbar2,dCR2,'r-o'); hold on;
h6=errorbar(g1,CRbar4,dCR4,'m-x'); hold on;
h7=errorbar(g1,CRbar5,dCR5,'k-+');
set([h5 h6 h7],'LineWidth',1.2)
%set(gca,'xticklabel',arrayfun(@(x) num2str(x),log10(get(gca,'xtick')),'un',0))
ylim ([0 1])
legend('r = 1/\tau','r = 1/3\tau','r = 1/10\tau','Location','northwest')
set(gca,'xscale','log')
xlabel('g(nM/mm)')
ylabel('CR')
title('CR')
hold off
