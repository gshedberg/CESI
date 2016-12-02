n = 100;%# of iterations
m = 0;
%% Reiterate for large plots
if m ==1  %decision on whether to run constant vs varying recovery
    recovery = zeros(n);
    Efficiency = zeros(n);
    Eff_FC = zeros(n);
    Eff_GT = zeros(n);
    W_net = zeros(n);
    Wfc_vec = zeros(n);
    W_gt = zeros(n);
    Wc2 = zeros(n);
    T_out = zeros(n);       %Initialization of all output matrices
    TFlowOut = zeros(n);
    ReactMix = zeros(n);
    V_vec = zeros(n);
    R_actual = zeros(n);
    Rt = zeros(n);
    Qextra = zeros(n);
    i_array = zeros(n);
    recirc_vec = zeros(n);
    for j= 1:n
        %ITM back pressure in kPa
        P_ITMperm = linspace(50,50)'; 
        %Fixed value of recovery
        recovery = linspace(.99*(j/n),.99*(j/n))'; 
        TIT = linspace(1200,1200)';
        %Mass Flow of GT
        Mflow = linspace(15,15); 
        Fuel = 1; % 0 for no suplemental fuel into combustor, 1 for fixed % recovery
        Pr = linspace(15,15)'; % Compressor pressure ratio
        vectorLength = max([length(Pr), length(P_ITMperm),length(recovery)]); %set length of vectors to correspond to given inputs
        Tin = zeros(vectorLength,9)+ 300;
        if Fuel==1
            [Efficiency(:,j),Eff_FC(:,j),Eff_GT(:,j),W_net(:,j),Wfc_vec(:,j),...
                W_gt(:,j),Wc2(:,j),T_out(:,j),TFlowOut(:,j),ReactMix(:,j),...
                V_vec(:,j),R_actual(:,j),Rt(:,j),recovery(:,j),Qextra(:,j),i_array(:,j),recirc_vec(:,j)]...
                = hybrid_final_struc(Tin,Pr,P_ITMperm,TIT,Mflow,recovery);
        else
            [Efficiency(:,j),Eff_FC(:,j),Eff_GT(:,j),W_net(:,j),Wfc_vec(:,j),...
                W_gt(:,j),Wc2(:,j),T_out(:,j),TFlowOut(:,j),ReactMix(:,j),...
                V_vec(:,j),R_actual(:,j),Rt(:,j),recovery(:,j),Qextra(:,j),i_array(:,j),recirc_vec(:,j)]...
                = hybrid_final_struc(Tin,Pr,P_ITMperm,TIT,Mflow);            
        end
    end
%% Original model that outputs all vectors [100,1]
else 
    %ITM back pressure in kPa
    P_ITMperm = linspace(25,25)'; 
    %Fixed value of recovery
    recovery = linspace(.51,.51)'; 
    TIT = linspace(1200,1200)';
    %Mass Flow of GT
    Mflow = linspace(20,20)'; 
    % GT pressure ratio
    Pr = linspace(15,15)';
    iDen = linspace(.25,1)';
    %Initialize
    Fuel = 1; % 0 for no suplemental fuel into combustor, 1 for fixed % recovery
    vectorLength = max([length(Pr), length(P_ITMperm)]); %set length of vectors to correspond to given inputs
    Tin = zeros(vectorLength,1)+ 300;
    if Fuel==1
       [Efficiency,Eff_FC,Eff_GT,W_net,Wfc_vec,W_gt,Wc2,T_out,TFlowOut,...
           ReactMix,V_vec,Utilization,R_actual,Rt,recovery,Qextra,i_array,recirc_vec] = hybrid_final_struc(Tin,Pr,P_ITMperm,TIT,Mflow,recovery,iDen);
    else
       [Efficiency,Eff_FC,Eff_GT,W_net,Wfc_vec,W_gt,Wc2,T_out,TFlowOut,...
           ReactMix,V_vec,Utilization,R_actual,Rt,recovery,Qextra,i_array,recirc_vec] = hybrid_final_struc(Tin,Pr,P_ITMperm,TIT,Mflow);
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
    line(x,y,'linewidth',3,'color','r');
    ylabel('Percent Oxygen Recovered (%)')
    ax.XTickLabel = {'.760','.780','.800','.820','.840','.860','.880','.900','.920'}; %Voltage
%   ax.XTickLabel ={'.625','.651','.678','.707','.739','.774','.812','.853','.899'}; %Utilization
%   ax.XTickLabel ={'.300','.278','.258','.238','.218','.198','.178','.158','.138'};%Current Density
end

