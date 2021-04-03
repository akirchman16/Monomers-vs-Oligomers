% clearvars -except DNA;
clearvars;
% close all;

% This code takes solution of RAD51 in equilibrium between monomers and
% dimers and performs a Gillespie lattice model with the corresponding
% solution. The concentrations of RAD51 in solution is assumed to remain
% constant. The First-Reaction method is used for the Gillespie Algorithm.

N = 8660;
n = 3;
w = 1;
L_Total = 2;    %total concentration of RAD51
Percent_Monomer = [0,1];   %Percentage of solution which is monomers
k_on = 1;   %kinetic rate constants
k_off = 1;

minIterations = 1000;

%Memory Allocation
EventFractions = zeros(numel(Percent_Monomer),7);
FracCover = zeros(numel(Percent_Monomer),minIterations);
t = zeros(numel(Percent_Monomer),minIterations);
Max_Time = zeros(1,numel(Percent_Monomer));
Equilibrium_Coverage = zeros(1,numel(Percent_Monomer));
L_Monomer = zeros(1,numel(Percent_Monomer));
L_Dimer = zeros(1,numel(Percent_Monomer));

Loops = 0;
for Ratio = Percent_Monomer
    Loops = Loops+1;
    
    L_Monomer(Loops) = Ratio*L_Total;        %Concentration of monomer RAD51
    L_Dimer(Loops) = (1-Ratio)*L_Total;    %Concentration of dimer RAD51

%     DNA(1) = [];    %clearing ends of DNA to compare counting method to the old model
%     DNA(N+1) = [];
    
    DNA = zeros(1,N);
    BoundAtSpot = zeros(1,N);   %records where monomers are bound on lattice

    %Memroy Allocations
    Populations = zeros(minIterations,7);
    a = zeros(minIterations,7);
    Probabilities = zeros(minIterations,7);
    FiringAmounts = zeros(minIterations,7);
    dt = zeros(1,minIterations);
    j = zeros(1,minIterations);
    Location_History = zeros(7,minIterations);

    Equilibrium = 0;
    Events = 0;
    while ~Equilibrium
        Events = Events+1;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Vectors to store locations of each site
        Isolated_M = 0;
        Singly_Contiguous_M = 0;
        Doubly_Contiguous_M = 0;
        Isolated_D = 0;
        Singly_Contiguous_D = 0;
        Doubly_Contiguous_D = 0;
        
        if (DNA(1:1+(n-1)) == 0) & (DNA(1+n) == 0) %monomer isolated at first location
            Isolated_M = [Isolated_M,1];
        elseif (DNA(1:1+(n-1)) == 0) & (DNA(1+n) == 1) %monomer singly contiguous at first location
            Singly_Contiguous_M = [Singly_Contiguous_M,1];
        end
        if (DNA(1:1+(2*n-1)) == 0) & (DNA(1+2*n) == 0) %dimer isolated at first location
            Isolated_D = [Isolated_D,1];
        elseif (DNA(1:1+(2*n-1)) == 0) & (DNA(1+2*n) == 1) %dimer singly contiguous at first location
   
            Singly_Contiguous_D = [Singly_Contiguous_D,1];
        end
        if (DNA(N-(n-1):N) == 0) & (DNA(N-n) == 0) %monomer isolated at last location
            Isolated_M = [Isolated_M,N-(n-1)];
        elseif (DNA(N-(n-1):N) == 0) & (DNA(N-n) == 1) %monomer singly contiguous at last location
            Singly_Contiguous_M = [Singly_Contiguous_M,N-(n-1)];
        end
        if (DNA(N-(2*n-1):N) == 0) & (DNA(N-2*n) == 0) %dimer isolated at last location
            Isolated_D = [Isolated_D,2*n-1];
        elseif (DNA(N-(2*n-1):N) == 0) & (DNA(N-2*n) == 1) %dimer singly contiguous at last location
            Singly_Contiguous_D = [Singly_Contiguous_D,2*n-1];
        end
        for y = 2:N-(n-1)-1    %check for monomer binding sites
            if DNA(y:y+(n-1)) == 0
                if DNA(y-1) == 0 & DNA(y+n) == 0  %isolated monomer site
                    Isolated_M = [Isolated_M,y];
                elseif (DNA(y-1) == 0 & DNA(y+n) == 1) | (DNA(y-1) == 1 & DNA(y+n) == 0)    %singly contiguous monomer site
                    Singly_Contiguous_M = [Singly_Contiguous_M,y];
                elseif DNA(y-1) == 1 & DNA(y+n) == 1   %doubly contiguous monomer site
                    Doubly_Contiguous_M = [Doubly_Contiguous_M,y];
                end
            end
        end
        for z = 2:N-(2*n-1)-1   %check for dimer binding sites
            if DNA(z:z+(2*n-1)) == 0
                if DNA(z-1) == 0 & DNA(z+2*n) == 0 %isolated dimer site
                    Isolated_D = [Isolated_D,z];
                elseif (DNA(z-1) == 0 & DNA(z+2*n) == 1) | (DNA(z-1) == 1 & DNA(z+2*n) == 0) %singly contiguous dimer site
                    Singly_Contiguous_D = [Singly_Contiguous_D,z];
                elseif DNA(z-1) == 1 & DNA(z+2*n) == 1 %doubly contiguous dimer site
                    Doubly_Contiguous_D = [Doubly_Contiguous_D,z];
                end
            end
        end
%         clearing all zeros stored at the beginning of the search process
        Isolated_M(Isolated_M == 0) = [];
        Singly_Contiguous_M(Singly_Contiguous_M == 0) = [];
        Doubly_Contiguous_M(Doubly_Contiguous_M == 0) = [];
        Isolated_D(Isolated_D == 0) = [];
        Singly_Contiguous_D(Singly_Contiguous_D == 0) = [];
        Doubly_Contiguous_D(Doubly_Contiguous_D == 0) = [];
%         organzing list of locations
        Isolated_M = sort(Isolated_M);
        Singly_Contiguous_M = sort(Singly_Contiguous_M);
        Doubly_Contiguous_M = sort(Doubly_Contiguous_M);
        Isolated_D = sort(Isolated_D);
        Singly_Contiguous_D = sort(Singly_Contiguous_D);
        Doubly_Contiguous_D = sort(Doubly_Contiguous_D);        
%         population numbers from search process
        xB_IM = length(Isolated_M);
        xB_SCM = length(Singly_Contiguous_M);
        xB_DCM = length(Doubly_Contiguous_M);
        xB_ID = length(Isolated_D);
        xB_SCD = length(Singly_Contiguous_D);
        xB_DCD = length(Doubly_Contiguous_D);
        xAB = sum(DNA)/n;
        
        Populations(Events,:) = [xB_IM,xB_SCM,xB_DCM,xB_ID,xB_SCD,xB_DCD,xAB];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        a(Events,:) = Populations(Events,:).*[k_on*L_Monomer(Loops),k_on*L_Monomer(Loops)*w,k_on*L_Monomer(Loops)*(w^2),k_on*L_Dimer(Loops),k_on*L_Dimer(Loops)*w,k_on*L_Dimer(Loops)*(w^2),k_off];
        Probabilities(Events,:) = a(Events,:)./sum(a(Events,:));

        r = [rand,rand,rand,rand,rand,rand,rand];
        tau = (1./a(Events,:)).*log(1./r);
        FiringAmounts(Events,:) = a(Events,:).*tau;
        dt(Events) = tau(tau == min(tau));
        j(Events) = find(tau == min(tau));

        if j(Events) == 1        %Isolated monomer binding
            Bind_Spot_I_M = Isolated_M(randi(xB_IM));
            DNA(Bind_Spot_I_M:Bind_Spot_I_M+(n-1)) = 1;
            Location_History(1,Events) = Bind_Spot_I_M;
            BoundAtSpot(Bind_Spot_I_M) = 1;
        elseif j(Events) == 2    %Singly-Contiguous monomer binding
            Bind_Spot_SC_M = Singly_Contiguous_M(randi(xB_SCM));
            DNA(Bind_Spot_SC_M:Bind_Spot_SC_M+(n-1)) = 1;
            Location_History(2,Events) = Bind_Spot_SC_M;
            BoundAtSpot(Bind_Spot_SC_M) = 1;
        elseif j(Events) == 3    %Doubly-Contiguous monomer binding
            Bind_Spot_DC_M = Doubly_Contiguous_M(randi(xB_DCM));
            DNA(Bind_Spot_DC_M:Bind_Spot_DC_M+(n-1)) = 1;
            Location_History(3,Events) = Bind_Spot_DC_M;
            BoundAtSpot(Bind_Spot_DC_M) = 1;
        elseif j(Events) == 4    %Isolated dimer binding
            Bind_Spot_I_D = Isolated_D(randi(xB_ID));
            DNA(Bind_Spot_I_D:Bind_Spot_I_D+(2*n-1)) = 1;
            Location_History(4,Events) = Bind_Spot_I_D;
            BoundAtSpot(Bind_Spot_I_D) = 1;
            BoundAtSpot(Bind_Spot_I_D+n) = 1;
        elseif j(Events) == 5    %Singly-Contiguous dimer binding
            Bind_Spot_SC_D = Singly_Contiguous_D(randi(xB_SCD));
            DNA(Bind_Spot_SC_D:Bind_Spot_SC_D+(2*n-1)) = 1;
            Location_History(5,Events) = Bind_Spot_SC_D;
            BoundAtSpot(Bind_Spot_SC_D) = 1;
            BoundAtSpot(Bind_Spot_SC_D+n) = 1;
        elseif j(Events) == 6    %Doubly-Contiguous dimer binding
            Bind_Spot_DC_D = Doubly_Contiguous_D(randi(xB_DCD));
            DNA(Bind_Spot_DC_D:Bind_Spot_DC_D+(2*n-1)) = 1;
            Location_History(6,Events) = Bind_Spot_DC_D;
            BoundAtSpot(Bind_Spot_DC_D) = 1;
            BoundAtSpot(Bind_Spot_DC_D+1) = 1;
        elseif j(Events) == 7    %Monomer unbinding
            Bound_Locations = find(BoundAtSpot == 1);
            Unbind_Spot = Bound_Locations(length(Bound_Locations));
            DNA(Unbind_Spot:Unbind_Spot+(n-1)) = 0;
            Location_History(7,Events) = Unbind_Spot;
            BoundAtSpot(Unbind_Spot) = 0;
        end

        FracCover(Loops,Events+1) = sum(DNA)/N;
        t(Loops,Events+1) = t(Loops,Events)+dt(Events);

    %     Testing for Equilibrium
        if Events > minIterations && Equilibrium == 0
            FracCoverStates = (FracCover(Loops,(Events-floor(0.25*Events)):1:Events));
            FracCoverChange = abs(diff(FracCoverStates));   %Difference between each state for the last 1/4 of the simulation
            if (mean(FracCoverChange) <= 2*n/N) && (abs(FracCover(Loops,Events-floor(0.25*Events))-FracCover(Loops,Events)) <= 5*n/N)
                Equilibrium = 1;
%                 Irrelevant_t = find(t(Loops,2:end) == 0);
%                 t(Loops,Irrelevant_t) = NaN;
%                 Irrelevant_FracCover = find(FracCover(Loops,2:end) == 0);
%                 FracCover(Loops,Irrelevant_FracCover) = NaN;
            end
        end
    end

    EventFractions(Loops,:) = [numel(find(j==1)),numel(find(j==2)),numel(find(j==3)),numel(find(j==4)),numel(find(j==5)),numel(find(j==6)),numel(find(j==7))]./Events;
    Max_Time(Loops) = max(t(Loops,:));
    
    Equilibrium_Coverage(Loops) = mean(FracCoverStates);
    disp(['(Ratio = ', num2str(Ratio), ') - Equilibrium Saturation = ', num2str(Equilibrium_Coverage(Loops))]);
    
    figure(1);
    hold on;
    if Ratio == 1
        scatter(t(Loops,:),FracCover(Loops,:),3,'cyan','filled');
    elseif Ratio == 0
        scatter(t(Loops,:),FracCover(Loops,:),3,'magenta','filled');
    else
        scatter(t(Loops,:),FracCover(Loops,:),2,'filled');
    end
    xlabel('Time, t');
    ylabel('Fractional Coverage');
    xlim([0 max(t(Loops,:))]);
    ylim([0 1]);
%     yline(Equilibrium_Coverage(Loops),'k',['\rho = ', num2str(Percent_Monomer(Loops))]);
    title('Saturation of DNA Lattice');
end


Total_Events = zeros(1,length(Percent_Monomer));
Legend = cell(length(Percent_Monomer),1);

for c = 1:length(Percent_Monomer)
    Legend{c} = ['\rho = ', num2str(Percent_Monomer(c))];
end
figure(1);
legend(Legend,'location','southeast');