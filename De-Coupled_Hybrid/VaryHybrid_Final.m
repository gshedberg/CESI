n = 100;%# of iterations
m = 0;
if m ==1  %decision on whether to run constant vs varying recovery
    recovery = zeros(n);
    Efficiency = zeros(n);
    Eff_FC = zeros(n);
    Eff_GT = zeros(n);
    W_net = zeros(n);
    W_fc = zeros(n);
    W_gt = zeros(n);
    T_out = zeros(n);       %Initialization of all output matrices
    N_out = zeros(n);
    V_fc = zeros(n);
    R_actual = zeros(n);
    FC_util = zeros(n);
    E0 = zeros(n);
    Qimbalance = zeros(n);
    n_h2 = zeros(n);
    for j= 1:n
        P_ITMperm = linspace(50,50)'; %ITM back pressure in kPa
        V_loss = linspace(.3,.12)'; %Fuel cell voltage
        recovery = linspace(.99*(j/n),.99*(j/n))'; %Fixed value of recovery
        TIT = linspace(1200,1200)';
        Fuel = 1; % 0 for no suplemental fuel into combustor, 1 for fixed % recovery
        Pr = linspace(15,15)'; % Compressor pressure ratio
        vectorLength = max([length(Pr), length(P_ITMperm),length(V_loss)]); %set length of vectors to correspond to given inputs

        TXNin = zeros(vectorLength,9);
        TXNin(:,1) = 300;
        Xin = [0 0 0 0 0 .79 .21];      %initial conditions of temp, x, and flow
        for i =1:1:7
            TXNin(:,i+1) = Xin(i);
        end
        TXNin(:,9) = .5;

        if Fuel==1
            [Efficiency(:,j),Eff_FC(:,j),Eff_GT(:,j),W_net(:,j),W_fc(:,j),W_gt(:,j),T_out(:,j),X8,N_out(:,j),V_fc,R_actual(:,j),recovery(:,j),FC_util(:,j),Qimbalance(:,j),E0(:,j),n_h2(:,j),Wc2(:,j),Rt(:,j)] = hybrid_final(TXNin,Pr, P_ITMperm, V_loss,TIT,recovery);
        else
            [Efficiency(:,j),Eff_FC(:,j),Eff_GT(:,j),W_net(:,j),W_fc(:,j),W_gt(:,j),T_out(:,j),X8,N_out(:,j),V_fc(:,j),R_actual(:,j),recovery(:,j),FC_util(:,j),Qimbalance(:,j),E0(:,j),n_h2(:,j),Wc2(:,j),Rt(:,j)] = hybrid_final(TXNin, Pr, P_ITMperm, V_loss,TIT);
        end
    end
    x = V_fc;
    y = linspace(.0099,.99);
    z = Efficiency;
    ax = gca;
    set(gca,'Fontsize',18)
    contour(x,y,z','showtext','on','linewidth',2)
    
else %original model that outputs all vectors [100,1]
    P_ITMperm = linspace(50,50)'; %ITM back pressure in kPa
    V_loss = linspace(.21,.21)'; %Fuel cell voltage
    recovery = linspace(.51,.51)'; %Fixed value of recovery
    TIT = linspace(1200,1200)';
    Fuel = 0; % 0 for no suplemental fuel into combustor, 1 for fixed % recovery
    Pr = linspace(15,15)'; % Compressor pressure ratio
    vectorLength = max([length(Pr), length(P_ITMperm),length(V_loss)]); %set length of vectors to correspond to given inputs

    TXNin = zeros(vectorLength,9);
    TXNin(:,1) = 300;
    Xin = [0 0 0 0 0 .79 .21];      %initial conditions of temp, x, and flow
    for i =1:1:7
        TXNin(:,i+1) = Xin(i);
    end
    TXNin(:,9) = .5;

    if Fuel==1
        [Efficiency,Eff_FC,Eff_GT,W_net,W_fc,W_gt,T_out,X8,N_out,V_fc,R_actual,recovery,FC_util,Qimbalance,E0,n_h2,Wc2,Rt] = hybrid_final(TXNin,Pr, P_ITMperm, V_loss,TIT,recovery);
    else
        [Efficiency,Eff_FC,Eff_GT,W_net,W_fc,W_gt,T_out,X8,N_out,V_fc,R_actual,recovery,FC_util,Qimbalance,E0,n_h2,Wc2,Rt] = hybrid_final(TXNin, Pr, P_ITMperm, V_loss,TIT);
    end
    figure(1)
    A = find(recovery>.99);
    if ~isempty(A)
       recovery(A) = 1;
    end
    ax = gca;
    set(gca,'FontSize',15)
% %     set(gca,'Xtick',.4:.2:3);
    x = V_fc;
    y = recovery;
    line(x,y,'linewidth',3,'color','r');
%     plot(x,y,'linewidth','3');
% 
% %     hold on
%     x2 = P_ITMperm;
%     y2 = R_actual; 
%   
%     [hAx,hLine1,hLine2] = plotyy(x,y,x2,y2);
% %     ax.XTickLabel = {};
% %     set(gca,'FontSize',18)
% %     set(hLine1,'linewidth',3)
% %     set(hLine2,'linewidth',3)
%     xlabel('Operating Voltage (V)')
%     ylabel(hAx(1),'Theoretical Recovery') % left y-axis
%     ylabel(hAx(2),'Actual Recovery') % left y-axis
%     title('Power Ratio (FC/GT) vs Fuel Cell Voltage')
%     line(x,y,'linewidth',3)
    
    ylabel('Percent Oxygen Recovered (%)')
    %xlabel('Operating Voltage')
    %ylabel('Percent Oxygen Recovered')
    ax.XTickLabel = {'.760','.780','.800','.820','.840','.860','.880','.900','.920'}; %Voltage
%     ax.XTickLabel ={'.625','.651','.678','.707','.739','.774','.812','.853','.899'}; %Utilization
    %ax.XTickLabel ={'.300','.278','.258','.238','.218','.198','.178','.158','.138'};%Current Density

end


% ax1 = gca;
% ax1_pos = ax1.Position; % position of first axes
% hold on
% FUtilLabel = interp1(V_fc,FC_util,(.84:1.04));
% num2str(FUtilLabel);
% ax.xtick = FUtilLabel;
% c = Eff_GT+Eff_GT.*y;
% d = 1+(Eff_GT.*(y./((x./E0).*FC_util)));
% z = (c*d');
% contour(x,y,z,'ShowText','on')
% xlabel('Operating Voltage and Fuel Utilization')
% ylabel('Power Ratio (FC/GT)')
% % constant = ones(100,1)*V_fc(81);
% % plot(constant,W_fc./W_gt)
% 
% % ax2 = axes('Position',ax1_pos,...
% %     'XAxisLocation','top',...
% %     'YAxisLocation','right',...
% %     'Color','none');
% % x1 = FC_util;
% % y1 = y;
% % hline1 = line(x1,y1,'Parent',ax2);
% % ylabel('FC Utility')
% title('Operating Voltage versus Power Ratio and FC Utilization')

% plot(x,y)
% [hAx,hLine1,hLine2] = plotyy(P_ITMperm,W_net,P_ITMperm,W_fc);
% xlabel('Back Pressure of OTM (kpa)')
% ylabel(hAx(1),'Net Power of Hybrid') % left y-axis
% ylabel(hAx(2),'Net Power of Fuel Cell') % left y-axis
% hLine1.LineStyle = '--';
% hLine2.LineStyle = ':';
% 
% figure(2)
% [hAx,hLine1,hLine2] = plotyy(P_ITMperm,W_net,P_ITMperm,W_gt);
% xlabel('Back Pressure of OTM (kpa)')
% ylabel(hAx(1),'Net Power of Hybrid')
% ylabel(hAx(2),'Net Power of Gas Turbine') 
% hLine1.LineStyle = '--';
% hLine2.LineStyle = ':';
% 
% figure(3)
% [hAx,hLine1,hLine2] = plotyy(P_ITMperm,Efficiency,P_ITMperm,Eff_FC);
% xlabel('Back Pressure of OTM (kpa)')
% ylabel(hAx(1),'Efficiency of Hybrid')
% ylabel(hAx(2),'Efficiency of Fuel Cell')
% hLine1.LineStyle = '--';
% hLine2.LineStyle = ':'; 
% 
% figure(4)
% [hAx,hLine1,hLine2] = plotyy(P_ITMperm,Efficiency,P_ITMperm,Eff_GT);
% xlabel('Back Pressure of OTM (kpa)')
% ylabel(hAx(1),'Efficiency of Hybrid')
% ylabel(hAx(2),'Efficiency of Gas Turbine')
% hLine1.LineStyle = '--';
% hLine2.LineStyle = ':';

% figure(5)
% [hAx,hLine1,hLine2] = plotyy(V_fc,Efficiency,V_fc,W_fc);
% xlabel('Fuel Cell Output Voltage')
% ylabel(hAx(1),'Efficiency')
% ylabel(hAx(2),'W_fc')
% hLine1.LineStyle = '--';
% hLine2.LineStyle = ':';
% 
