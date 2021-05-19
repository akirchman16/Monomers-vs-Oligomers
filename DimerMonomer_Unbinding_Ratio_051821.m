% clearvars -except DNA;
clearvars;
close all;

% This code takes solution of RAD51 in equilibrium between monomers and
% dimers and performs a Gillespie lattice model with the corresponding
% solution. The concentrations of RAD51 in solution is assumed to remain
% constant. This is the same as the previous ratio model except that both
% monomers and dimers are allowed to unbind. In the previous model only
% monomers were allowed to unbind from the ssDNA lattice. The
% First-Reaction method is used for the Gillespie Algorithm.

N = 8660;
n = 3;  %length of a monomer
w = 1;
L_Total = 2;    %total concentration of RAD51
Percent_Monomer = [0,0.5,1];   %Percentage of solution which is monomers
k_on = 1;   %kinetic rate constants
k_off = 1;

minIterations = 1000;

%Memory Allocation
EventFractions = zeros(numel(Percent_Monomer),8);
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
    
    DNA = zeros(1,N);
    BoundAtSpot = zeros(1,N);   %records where monomers are bound on lattice

    %Memroy Allocations
    Populations = zeros(minIterations,8);
    a = zeros(minIterations,8);
    Probabilities = zeros(minIterations,8);
    FiringAmounts = zeros(minIterations,8);
    dt = zeros(1,minIterations);
    j = zeros(1,minIterations);
    Location_History = zeros(8,minIterations);

    Equilibrium = 0;
    Events = 0;
    while max(t(Loops,:)) < 1.5
%         if max(t(Loops,:)) >= 1.25
%             L_Monomer(Loops) = 0;
%             L_Dimer(Loops) = 0;
%         end
        Events = Events+1;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        searcH_vars = {'Left','Right','Gap_Size','Left_Available_M','Right_Available_M','Left_Available_D','Right_Available_D','Right_BindingSite_M','Right_BindingSite_D','Gaps_L2_M','Gaps_R2_M','Gaps_L2_D','Gaps_R2_D','Gap_Size2_M','Gap_Size2_D','Gap_SizeI_M','Gap_SizeI_D','Left_I_M','Left_I_D','Isolated_M','SinglyContiguous_M','Doubly_Contiguous_M','Isolated_D','SinglyContiguous_D','DoublyContiguous_D','Left_Filaments','Right_Filaments','Filament_Lengths','Left_Dimer_Filaments','Dimer_Filament_Lengths','Monomers_per_Dimer_Filaments','Left_Dimer_Locations'};
        clear searcH_vars;

        Left = find(diff([1 DNA 1]) == -1);
        Right = find(diff([1 DNA 1]) == 1)-1;
        Gap_Size = Right-Left+1;

        Left_Available_M = Left(Gap_Size >= n);
        Right_Available_M = Right(Gap_Size >= n);
        Left_Available_D = Left(Gap_Size >= 2*n);
        Right_Available_D = Right(Gap_Size >= 2*n);
        Right_BindingSite_M = Right_Available_M-(n-1);
        Right_BindingSite_D = Right_Available_D-(2*n-1);

        %Doubly Contiguous Searches
        Doubly_Contiguous_M = Left_Available_M(Left_Available_M == Right_BindingSite_M & 1 < Left_Available_M & Left_Available_M < N-(n-1));
        Doubly_Contiguous_D = Left_Available_D(Left_Available_D == Right_BindingSite_D & 1 < Left_Available_D & Left_Available_D < N-(2*n-1));

        %Singly Contiguous Searches
        Singly_Contiguous_M = unique([Left_Available_M(~ismember(Left_Available_M,Doubly_Contiguous_M)),Right_Available_M(~ismember(Right_Available_M-(n-1),Doubly_Contiguous_M))-(n-1)]);
        Singly_Contiguous_M(Singly_Contiguous_M == N-(n-1) | Singly_Contiguous_M == 1) = [];
        Singly_Contiguous_D = unique([Left_Available_D(~ismember(Left_Available_D,Doubly_Contiguous_D)),Right_Available_D(~ismember(Right_Available_D-(2*n-1),Doubly_Contiguous_D))-(2*n-1)]);
        Singly_Contiguous_D(Singly_Contiguous_D == N-(2*n-1) | Singly_Contiguous_D == 1) = [];

        %Isolated Searches
        Gaps_L2_M = Left(~ismember(Left,Doubly_Contiguous_M));
        Gaps_L2_D = Left(~ismember(Left,Doubly_Contiguous_M) & ~ismember(Left,Doubly_Contiguous_D));
        Gaps_R2_M = Right(~ismember(Right,Doubly_Contiguous_M+(n-1)));
        Gaps_R2_D = Right(~ismember(Right,Doubly_Contiguous_M+(n-1)) & ~ismember(Right,Doubly_Contiguous_D+(2*n-1)));
        Gap_Size2_M = Gaps_R2_M-Gaps_L2_M+1;
        Gap_Size2_D = Gaps_R2_D-Gaps_L2_D+1;

        Gap_SizeI_M = Gap_Size2_M(Gap_Size2_M > n); %size of gaps where isolated binding can happen
        Gap_SizeI_D = Gap_Size2_D(Gap_Size2_D > 2*n);
        Left_I_M = Gaps_L2_M(ismember(Gap_Size2_M,Gap_SizeI_M));    %left position of gaps where isolated binding can happen
        Left_I_D = Gaps_L2_D(ismember(Gap_Size2_D,Gap_SizeI_D));

        Isolated_M = zeros(1,sum(Gap_SizeI_M)-((n+1)*numel(Gap_SizeI_M))+logical(ismember(N,Left_I_M+Gap_SizeI_M-1))+logical(sum(DNA) == 0));   %memory allocation
        Isolated_D = zeros(1,sum(Gap_SizeI_D)-((2*n+1)*numel(Gap_SizeI_D))+logical(ismember(N,Left_I_D+Gap_SizeI_D-1))+logical(sum(DNA) == 0));
        for i = 1:numel(Left_I_M)
            Isolated_M(find(Isolated_M == 0,1):find(Isolated_M == 0,1)+Gap_SizeI_M(i)-(n+1)-1+logical(N == Left_I_M(i)+Gap_SizeI_M(i)-1)+logical(Left_I_M(i) == 1)) = Left_I_M(i)+(1-logical(Left_I_M(i) == 1):Gap_SizeI_M(i)-(n+1)+logical(N == Left_I_M(i)+Gap_SizeI_M(i)-1));
        end
        for k = 1:numel(Left_I_D)
            Isolated_D(find(Isolated_D == 0,1):find(Isolated_D == 0,1)+Gap_SizeI_D(k)-(2*n+1)-1+logical(N == Left_I_D(k)+Gap_SizeI_D(k)-1)+logical(Left_I_D(k) == 1)) = Left_I_D(k)+(1-logical(Left_I_D(k) == 1):Gap_SizeI_D(k)-(2*n+1)+logical(N == Left_I_D(k)+Gap_SizeI_D(k)-1));
        end

        %Search for dimers
        Left_Filaments = find(diff([0 DNA 0]) == 1);    %left ends of all filaments
        Right_Filaments = find(diff([0 DNA 0]) == -1)-1;   %right ends of all filaments
        Filament_Lengths = Right_Filaments-Left_Filaments+1;    %lengths of all filaments
        Left_Dimer_Filaments = Left_Filaments(Filament_Lengths > n);   %left end of filament that contains at least a single dimer
        Dimer_Filament_Lengths = Filament_Lengths(Filament_Lengths > n);    %lengths of alll filaments that contain a dimer
        Monomers_per_Dimer_Filaments = Dimer_Filament_Lengths./n;    %number of monomers in each dimer filament
        Left_Dimer_Locations = [];
        if numel(Dimer_Filament_Lengths) > 0
            for f = 1:numel(Dimer_Filament_Lengths)
                if Dimer_Filament_Lengths(f) == 2*n    %if there's only one dimer, add only the end of the filament
                    Left_Dimer_Locations = [Left_Dimer_Locations, Left_Dimer_Filaments(f)];
                else   %otherwise add at each monomer location
                    for g = 1:Monomers_per_Dimer_Filaments(f)-1
                        Left_Dimer_Locations = [Left_Dimer_Locations, Left_Dimer_Filaments(f)+((g-1)*n)];   %adds location of each possible dimer in a long filament
                    end
                end
            end
        end
        
        %Population Numbers
        xB_IM = length(Isolated_M);
        xB_SCM = length(Singly_Contiguous_M);
        xB_DCM = length(Doubly_Contiguous_M);
        xB_ID = length(Isolated_D);
        xB_SCD = length(Singly_Contiguous_D);
        xB_DCD = length(Doubly_Contiguous_D);
        xAB_M = sum(DNA)/n;
        xAB_D = numel(Left_Dimer_Locations);

        Populations(Events,:) = [xB_IM,xB_SCM,xB_DCM,xB_ID,xB_SCD,xB_DCD,xAB_M,xAB_D];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        a(Events,:) = Populations(Events,:).*[k_on*L_Monomer(Loops),k_on*L_Monomer(Loops)*w,k_on*L_Monomer(Loops)*(w^2),k_on*L_Dimer(Loops),k_on*L_Dimer(Loops)*w,k_on*L_Dimer(Loops)*(w^2),k_off,k_off];
        Probabilities(Events,:) = a(Events,:)./sum(a(Events,:));

        r = [rand,rand,rand,rand,rand,rand,rand,rand];
        tau = (1./a(Events,:)).*log(1./r);
        FiringAmounts(Events,:) = a(Events,:).*tau;
        dt(Events) = min(tau);
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
            BoundAtSpot(Bind_Spot_DC_D+n) = 1;
        elseif j(Events) == 7    %Monomer unbinding
            Bound_Locations = find(BoundAtSpot == 1);
            Unbind_Spot = Bound_Locations(randi(length(Bound_Locations)));
            DNA(Unbind_Spot:Unbind_Spot+(n-1)) = 0;
            Location_History(7,Events) = Unbind_Spot;
            BoundAtSpot(Unbind_Spot) = 0;
        elseif j(Events) == 8   %Dimer unbinding
            Dimer_Unbind_Spot = Left_Dimer_Locations(randi(length(Left_Dimer_Locations)));
            DNA(Dimer_Unbind_Spot:Dimer_Unbind_Spot+((2*n)-1)) = 0;
            Location_History(8,Events) = Dimer_Unbind_Spot;
            BoundAtSpot(Dimer_Unbind_Spot) = 0;
            BoundAtSpot(Dimer_Unbind_Spot+n) = 0;
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
        
        if ~isempty(find(BoundAtSpot == 1 & DNA ~= 1, 1))
            disp('STOP');
            flag = 1;
            break
        end
    end

    EventFractions(Loops,:) = [numel(find(j==1)),numel(find(j==2)),numel(find(j==3)),numel(find(j==4)),numel(find(j==5)),numel(find(j==6)),numel(find(j==7)),numel(find(j==8))]./Events;
    Max_Time(Loops) = max(t(Loops,:));
      
    figure(1);
    hold on;
    if Ratio == 1
        scatter(t(Loops,:),FracCover(Loops,:),5,'cyan','filled');
    elseif Ratio == 0
        scatter(t(Loops,:),FracCover(Loops,:),5,'magenta','filled');
    else
        scatter(t(Loops,:),FracCover(Loops,:),4,'filled');
    end
    xlabel('Time, t');
    ylabel('Fractional Coverage');
    xlim([0 max(t(Loops,:))]);
    ylim([0 1]);
    title('Saturation of DNA Lattice');
    
    if flag == 1
        break
    end
end

Total_Events = zeros(1,length(Percent_Monomer));
Legend = cell(length(Percent_Monomer),1);
for c = 1:length(Percent_Monomer)
    Legend{c} = ['\rho = ', num2str(Percent_Monomer(c))];
end
figure(1);
legend(Legend,'location','southeast');