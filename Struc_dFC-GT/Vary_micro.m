n = 100;%# of iterations
m = 1;
%% Reiterate for large plots
if m ==0  %decision on whether to run constant vs varying recovery
    recovery_m = zeros(n);
    Efficiency_m = zeros(n);
    Eff_FC_m = zeros(n);
    Eff_GT_m = zeros(n);
    W_net_m = zeros(n);
    Wfc_vec_m = zeros(n);
    W_gt_m = zeros(n);
    Wc2_m = zeros(n);
    T_out_m = zeros(n);       %Initialization of all output matrices
%     TFlowOut = [];
%     ReactMix = [];
    V_vec_m = zeros(n);
    R_actual_m = zeros(n);
    Rt_m = zeros(n);
    Qextra_m = zeros(n);
    Utilization = zeros(n);
%     i_array = [];
    recirc_vec_m = zeros(n);
    const_recovery = linspace(.05,.99);
    for k= 1:n
        %ITM back pressure in kPa
        P_ITMperm = linspace(50,50)'; 
        %Fixed value of recovery
        recovery = ones(100,1).*const_recovery(k);
        TIT = linspace(1200,1200)';
        %Mass Flow of GT
        Mflow = linspace(5,5)'; 
        Fuel = 1; % 0 for no suplemental fuel into combustor, 1 for fixed % recovery
        Pr = linspace(3,5)'; % Compressor pressure ratio
        %Average Current Density for FC
        iDen = linspace(.6,.6)';
        vectorLength = max([length(Pr), length(P_ITMperm),length(recovery)]); %set length of vectors to correspond to given inputs
        Tin = zeros(vectorLength,1)+ 300;
        if Fuel==1
            [Eff,Eff_FC,Eff_GT,W_net,Wfc_vec,...
                W_gt,Wc2,T_out,~,~,...
                V_vec,Util,R_actual,Rt,recovery,Qextra,~,recirc_vec]...
                = microturb(Tin,Pr,P_ITMperm,TIT,Mflow,iDen,recovery);
            Efficiency_m(:,k) = Eff;
            Eff_FC_m(:,k) = Eff_FC;
            Eff_GT_m(:,k) = Eff_GT;
            W_net_m(:,k) = W_net;
            Wfc_vec_m(:,k) = Wfc_vec;
            W_gt_m(:,k) = W_gt;
            Wc2_m(:,k) = Wc2;
            T_out_m(:,k) = T_out;  
            V_vec_m(:,k) = V_vec;
            R_actual_m(:,k) = R_actual;
            Rt_m(:,k) = Rt;
            Qextra_m(:,k) = Qextra;
            Utilization(:,k) = Util;
            recirc_vec_m(:,k) = recirc_vec;
%             [Efficiency(:,k),Eff_FC(:,k),Eff_GT(:,k),W_net(:,k),Wfc_vec(:,k),...
%                 W_gt(:,k),Wc2(:,k),T_out(:,k),~,~,...
%                 V_vec(:,k),R_actual(:,k),Rt(:,k),recovery(:,k),Qextra(:,k),~,recirc_vec(:,k)]...
%                 = hybrid_final_struc(Tin,Pr,P_ITMperm,TIT,Mflow,iDen,recovery);
        else
            [Eff(:,k),Eff_FC(:,k),Eff_GT(:,k),W_net(:,k),Wfc_vec(:,k),...
                W_gt(:,k),Wc2(:,k),T_out(:,k),~,~,...
                V_vec(:,k),R_actual(:,k),Rt(:,k),recovery(:,k),Qextra(:,k),~,recirc_vec(:,k)]...
                = microturb(Tin,Pr,P_ITMperm,TIT,Mflow,iDen);            
        end
    end
%     contour(V_vec,const_recovery
%% Original model that outputs all vectors [100,1]
else 
    %ITM back pressure in kPa
    P_ITMperm = linspace(50,50)'; 
    %Fixed value of recovery
    recovery = linspace(.51,.51)'; 
    TIT = linspace(1200,1200)';
    %Mass Flow of GT
    Mflow = linspace(6,6)'; 
    % GT pressure ratio
    Pr = linspace(5,5)';
    %Average Current Density to Determine # of Cells in FC
    iDen = linspace(.5,.5)';
    %Initialize
    Fuel = 0; % 0 for no suplemental fuel into combustor, 1 for fixed % recovery
    vectorLength = max([length(Pr), length(P_ITMperm)]); %set length of vectors to correspond to given inputs
    Tin = zeros(vectorLength,1)+ 300;
    if Fuel==1
       [Eff,Eff_FC,Eff_GT,W_net,Wfc_vec,W_gt,Wc2,T_out,TFlowOut,...
           ReactMix,V_vec,Utilization,R_actual,Rt,recovery,Qextra,i_array,recirc_vec]...
           = microturb(Tin,Pr,P_ITMperm,TIT,Mflow,iDen,recovery);
    else
       [Eff,Eff_FC,Eff_GT,W_net,Wfc_vec,W_gt,Wc2,T_out,TFlowOut,...
           ReactMix,V_vec,Utilization,R_actual,Rt,recovery,Qextra,i_array,recirc_vec]...
           = microturb(Tin,Pr,P_ITMperm,TIT,Mflow,iDen);
    end
    figure(1)
    A = find(recovery>.99);
    if ~isempty(A)
       recovery(A) = 1;
    end
    ax = gca;
    set(gca,'FontSize',15)
% %     set(gca,'Xtick',.4:.2:3);
    x = V_vec;
    y = recovery;
    line(x,y,'linewidth',3,'color','r')
    ylabel('Percent Oxygen Recovered (%)')
    xlabel('Operating Voltage')
%     ax.XTickLabel = {'.760','.780','.800','.820','.840','.860','.880','.900','.920'}; %Voltage
%   ax.XTickLabel ={'.625','.651','.678','.707','.739','.774','.812','.853','.899'}; %Utilization
%   ax.XTickLabel ={'.300','.278','.258','.238','.218','.198','.178','.158','.138'};%Current Density
end

