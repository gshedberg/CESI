n = 100;%# of iterations
m = 0;
if m ==1  %decision on whether to run constant vs varying recovery
    recovery = zeros(n);
    Psi_comp1_dest = zeros(n);
    Psi_OTM_dest = zeros(n);
    Psi_comp2_dest = zeros(n);
    Psi_FC_dest = zeros(n);
    Psi_combust_dest = zeros(n);
    Psi_turbine_dest = zeros(n);
    V_fc = zeros(n);
    Tstoich = zeros(n);
    Wfc = zeros(n);
    HinFC = zeros(n);
    Efficiency = zeros(n);
    for j= 1:n
        P_ITMperm = linspace(50,50)'; %ITM back pressure in kPa
        V_loss = linspace(.3,.1)'; %Fuel cell voltage
        recovery = linspace(.99*(j/n),.99*(j/n))'; %Fixed value of recovery
        TIT = linspace(1200,1200)';
        Fuel = 1; % 0 for no suplemental fuel into combustor, 1 for fixed % recovery
        Pr = linspace(4,4)'; % Compressor pressure ratio
        vectorLength = max([length(Pr), length(P_ITMperm),length(V_loss)]); %set length of vectors to correspond to given inputs

        TXNin = zeros(vectorLength,9);
        TXNin(:,1) = 300;
        Xin = [0 0 0 0 0 .79 .21];      %initial conditions of temp, x, and flow
        for i =1:1:7
            TXNin(:,i+1) = Xin(i);
        end
        TXNin(:,9) = .5;

        if Fuel==1
            [Psi_comp1_dest(:,j),Psi_OTM_dest(:,j),Psi_comp2_dest(:,j),Psi_FC_dest(:,j),Psi_combust_dest(:,j),Psi_turbine_dest(:,j),V_fc(:,j),Tstoich(:,j),Wfc(:,j),HinFC(:,j),Efficiency(:,j)] = hybrid_exergy(TXNin,Pr, P_ITMperm, V_loss,TIT,recovery);
        else
            [Psi_comp1_dest(:,j),Psi_OTM_dest(:,j),Psi_comp2_dest(:,j),Psi_FC_dest(:,j),Psi_combust_dest(:,j),Psi_turbine_dest(:,j),V_fc(:,j),Tstoich(:,j),Wfc(:,j),HinFC(:,j),Efficiency(:,j)] = hybrid_exergy(TXNin,Pr, P_ITMperm, V_loss,TIT,recovery);
        end
    end
    x = V_fc;
    y = W_fc;
    z = Efficiency;
    contour(x,y,z','showtext','on','linewidth',2)
    
else %original model that outputs all vectors [100,1]
    P_ITMperm = linspace(50,50)'; %ITM back pressure in kPa
    V_loss = linspace(.3,.1)'; %Fuel cell voltage
    recovery = linspace(.3,.3)'; %Fixed value of recovery
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
        [Psi_comp1_dest,Psi_OTM_dest,Psi_comp2_dest,Psi_FC_dest,Psi_combust_dest,Psi_turbine_dest,V_fc,Tstoich,Wfc,Efficiency] = hybrid_exergy(TXNin,Pr, P_ITMperm, V_loss,TIT,recovery);
    else
        [Psi_comp1_dest,Psi_OTM_dest,Psi_comp2_dest,Psi_FC_dest,Psi_combust_dest,Psi_turbine_dest,V_fc,Tstoich,Wfc,Efficiency] = hybrid_exergy(TXNin,Pr, P_ITMperm, V_loss,TIT,recovery);    
    end
    A = find(recovery>.99);
    if ~isempty(A)
       recovery(A) = 1;
    end
    figure(1)
    x = V_fc;
    y = Psi_FC_dest;
    plot(x,y)
    xlabel('Voltage of Fuel Cell')
    ylabel('Exergy Destroyed')
    figure(2)
    x = V_fc;
    y = Efficiency;
    plot(x,y)
    xlabel('Voltage of Fuel Cel')
    ylabel('First Law Efficiency of the System')
    
end