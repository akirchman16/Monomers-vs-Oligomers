clearvars;
% close all;

% This code is a copy of the "Bidirectional Gillespie Model" that is in the
% "Gillespie Lattic Model" folder. The only difference is that no matter
% what size n is, unbinding will occur in n=3 unbinding units. This will
% make it possible to compare values of the old model between pure monomers
% or pure dimers. It will show whether the difference is causes by an error
% or if the unbinding unit is causing the issue.

N = 8660;    %length of DNA lattice
BindingUnit = [3,3,3,6,6,6];  %length of each protein (3 = monomer | 6 = dimer)
UnbindingUnit = 3;  %unbinding unit size (1 monomer)
k_on = 1;   %kinetic rate constant for binding (units: 1/(nt*s)??? )
k_off = 1;  %kinetic rate constant for unbinding (units: 1/s??? )
L = 2;    %concentration of free proteins (units: Molarity???)
w = 1;    %cooperativity constant

minIterations = 1000;  %minimum number of evnts that will occur
% K = k_on/k_off; %equilibrium constant (units: )

Loops = 0;

for n = BindingUnit
    Loops = Loops+1;

    DNA = zeros(1,N+2);   %empty DNA lattice with dummy zeros
    BoundAtSpot = zeros(1,N+2);  %empty CurrentBound array with dummy zeros

    xAB = zeros(1,minIterations+1); %memory allocation
    t = zeros(1,minIterations+1);
    a = zeros(minIterations+1,4);
    dt = zeros(1,minIterations);
    FracCover = zeros(1,minIterations+1);
    BindHist = zeros(4,minIterations); %1: isolated binding, 2: SC binding, 3: DC binding, 4: unbinding
    xB_I = zeros(1,minIterations);
    xB_SC = zeros(1,minIterations);
    xB_DC = zeros(1,minIterations);
    Length_nm = zeros(1,minIterations);

    t(1) = 0;   %all initial values for variables
    xB_I(1) = N-(n-1);    %initial values for a free lattice
    xB_SC(1) = 0;
    xB_DC(1) = 0;
    xAB(1) = 0;
    j = zeros(1,minIterations);
    FracCover(1) = sum(DNA)/N;
    x_Occupied(1) = 0;
    x_Unoccupied(1) = N;
    Length_nm(1) = 0.34*N;
    FracCoverStates = 0;
    FracCoverChange = 0;
    FilamentEnds_L = 0;
    FilamentEnds_R = 0;
    FilamentCount = 0;
    FilamentLengths = 0;
    Avg_FilamentLength = 0;
    BindCounter = 0;
    BindCounter_I = 0;
    BindCounter_SC = 0;
    BindCounter_DC = 0;
    UnbindCounter = 0;
    Equilibrium = 0;
    EqualTimes = 0;
    Events = 0;

    while (~Equilibrium || ~EqualTimes) && (t(end) < 1.5)
        Isolated = 0;   %preparing arrays for types of available locations
        SinglyContiguous = 0;
        DoublyContiguous = 0;
        Counter_I = 0;  %resets counters to count types of locations
        Counter_SC = 0;
        Counter_DC = 0;
        for x = 2:N-(n-1)+1
            if DNA(x:x+(n-1)) == 0
                if DNA(x-1) == 0 && DNA(x+n) == 0   %records all isolated locations
                    Isolated(Counter_I+1) = x;
                    Counter_I = Counter_I+1;
                elseif (DNA(x-1) == 0 && DNA(x+n) == 1) || (DNA(x-1) == 1 && DNA(x+n) == 0) %records all singly contiguous locations
                    SinglyContiguous(Counter_SC+1) = x;
                    Counter_SC = Counter_SC+1;
                elseif DNA(x-1) == 1 && DNA(x+n) == 1   %records all doubly contiguous locations
                    DoublyContiguous(Counter_DC+1) = x;
                    Counter_DC = Counter_DC+1;
                end
            end
        end
        xB_I(Events+1) = Counter_I;    %amounts of each location, used in propensity functions
        xB_SC(Events+1) = Counter_SC;
        xB_DC(Events+1) = Counter_DC;
        xAB(Events+1) = sum(DNA)/n;

        a(Events+1,:) = [k_on*L*xB_I(Events+1),k_on*L*xB_SC(Events+1)*w,k_on*L*xB_DC(Events+1)*(w^2),k_off*xAB(Events+1)];
        Rand = [rand,rand,rand,rand];
        tau = (1./a(Events+1,:)).*log(1./Rand);
        dt(Events+1) = min(tau);
        j(Events+1) = find(tau == min(tau));

        if j(Events+1) == 1 %tests for isolated binding event
            pos_I = randi(length(Isolated));
            SpotB_I = Isolated(pos_I); %chooses a random location for binding to occur
            DNA(SpotB_I:SpotB_I+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_I = BindCounter_I+1;
            BindHist(1,Events+1) = SpotB_I;
            BoundAtSpot(SpotB_I) = 1; %shows there's currently a protein bound at said location
            if n == 6
                BoundAtSpot(SpotB_I+(n/2)) = 1;   %additional location where monomer is bound from dimer binding
            end
        elseif j(Events+1) == 2    %tests for singly contiguous binding event
            pos_SC = randi(length(SinglyContiguous));
            SpotB_SC = SinglyContiguous(pos_SC);  %chooses a random location for binding
            DNA(SpotB_SC:SpotB_SC+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_SC = BindCounter_SC+1;
            BindHist(2,Events+1) = SpotB_SC;
            BoundAtSpot(SpotB_SC) = 1;    %shows there's currently a protein bound at location
            if n == 6
                BoundAtSpot(SpotB_SC+(n/2)) = 1;   %additional location where monomer is bound from dimer binding
            end
        elseif j(Events+1) == 3    %tests for doubly contiguous binding event
            pos_DC = randi(length(DoublyContiguous));
            SpotB_DC = DoublyContiguous(pos_DC);  %chooses a random location for binding
            DNA(SpotB_DC:SpotB_DC+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_DC = BindCounter_DC+1;
            BindHist(3,Events+1) = SpotB_DC;
            BoundAtSpot(SpotB_DC) = 1;    %shows there's currently a protein bound at location
            if n == 6
                BoundAtSpot(SpotB_DC+(n/2)) = 1;   %additional location where monomer is bound from dimer binding
            end
        elseif j(Events+1) == 4    %otherwise an unbinding event occurs
                CurrentBound = find(BoundAtSpot == 1);
                pos = randi(length(CurrentBound));
                SpotU = CurrentBound(pos);  %selects an already bound protein to unbind

                DNA(SpotU:SpotU+(UnbindingUnit-1)) = 0; %unbinds protein
                BoundAtSpot(SpotU) = 0;
                UnbindCounter = UnbindCounter+1;
                BindHist(4,Events+1) = SpotU;
        end

        FracCover((Events+1)+1) = sum(DNA)/N;  %fractional coverage

        if Events > minIterations && Equilibrium == 0
            FracCoverStates = (FracCover((Events-floor(0.25*Events)):1:Events));   %every 10 states for the last 20% of the model
            FracCoverChange = abs(diff(FracCoverStates));   %difference among each of those states
            if ((sum(FracCoverChange)/numel(FracCoverChange)) <= 2*n/N) && (abs(FracCover(Events-floor(0.25*Events))-FracCover(Events)) <= 3*n/N)
                Equilibrium = 1;
            else
                Equilibrium = 0;
            end
        end

        t((Events+1)+1) = t(Events+1)+dt(Events+1);
        Events = Events+1;

        if Loops > 1       %makes sure each cooperativity value runs apporximately the same amount of time
            if max(t) >= max(TotalTime)
                EqualTimes = 1;
            else
                EqualTimes = 0;
            end
        elseif Loops == 1
            EqualTimes = 1;
        end

        x_Occupied(Events) = length(find(DNA ~= 0));
        x_Unoccupied(Events) = length(find(DNA(2:N+1) == 0));
        Length_nm(Events) = (0.34*x_Unoccupied(Events))+(0.51*x_Occupied(Events));   %length of molecule in nm

        Probabilities(Events,:) = a(Events,:)./sum(a(Events));

        FilamentEnds_L = 1+find(diff(DNA)==1);  %locations of the left end of each filament
        FilamentEnds_R = 1+find(diff(DNA)==-1); %locations of the right end of each filament
        FilamentCount = length(FilamentEnds_R); %number of filaments on the lattice
        FilamentLengths(Events,1:FilamentCount) = FilamentEnds_R-FilamentEnds_L;    %length of filaments
        Avg_FilamentLength(Events) = sum(FilamentLengths(Events,1:FilamentCount))/length(FilamentLengths(Events,1:FilamentCount)); %average filament length for each event
    end

    Length_um = Length_nm/1e+3; %length of molecule in um
    Equilibrium_Coverage(Loops) = sum(FracCoverStates)/length(FracCoverStates);
%     disp(['(Binding Unit = ', num2str(n),') - Equilibrium Saturation = ', num2str(Equilibrium_Coverage(Loops))]);

    TotalTime(Loops) = max(t);
    LoopLengths = sort(TotalTime,'descend');

    EventFractions(Loops,:) = [length(find(j==1)),length(find(j==2)),length(find(j==3)),length(find(j==4))]./Events;

    figure(1);

%     Plot of fractional coverage of the DNA

%     subplot(2,1,1);
    if n == 3
        scatter(t,FracCover,2,'b','filled');
    elseif n == 6
        scatter(t,FracCover,2,'r','filled');
    end
    hold on;
    % plot(X,Y_exp_fit,'b');  %bounded exponential fit
    % plot(X,Y_sig_fit,'k');  %sigmoidal curve fit
    xlabel('Time, t');
    xlim([0 1.25]);
    ylabel('Saturation Level');
    ylim([0 1]);
    title('RAD51 Saturation of DNA');
    box on;

%     plot of Length of Molecule over time (using values of 0.34nm and
%     0.51nm from van der Heijden paper)

% %     subplot(2,1,2);
% %     scatter(t,Length_um,2,'filled');
% %     hold on;
% %     xlim([0 max(TotalTime)]);
% %     xlabel('Time, t (s)');
% %     ylabel('Length (\mum)');
% %     title('Length of DNA Molecule');

%     plot of average filament length over time

%     subplot(2,1,2);
%     scatter(t(2:length(t)),Avg_FilamentLength,2,'filled');
%     hold on;
%     xlabel('Time, t');
%     xlim([0 max(TotalTime)]);
%     ylabel('Average Filament Length (Sites)');
%     ylim([0 inf]);
%     title('Average Filament Length');

%     creates histogram of filament lengths in final state of model

%     figure(2);
%     subplot(1,length(Cooperativity),Loops);
%     histogram(FilamentLengths(Events,1:FilamentCount),15);
%     xlabel(['Mean = ', num2str(round(mean(FilamentLengths(Events,1:FilamentCount)),0))]);
%     ylabel('Occurence');
%     title('Filament Lengths of Final State');
end

% SortedEquilibrium = zeros(2,length(BindingUnit));
% if ~isempty(find(BindingUnit == 3, 1))
%     SortedEquilibrium(1,1:length(find(BindingUnit == 3))) = Equilibrium_Coverage(BindingUnit == 3);    %Equilibrium values for monomer only
%     Mean_Equilibrium_M = sum(SortedEquilibrium(1,:))/numel(find(BindingUnit == 3));  %Avg Equilibrium values for monomer only
%     yline(Mean_Equilibrium_M,'--k', ['\rho = 1 (', num2str(round(Mean_Equilibrium_M,3)), ')'],'LineWidth',1);
% end
% if ~isempty(find(BindingUnit == 6,1))
%     SortedEquilibrium(2,1:length(find(BindingUnit == 3))) = Equilibrium_Coverage(BindingUnit == 6);    %Equilibrium values for dimer only
%     Mean_Equilibrium_D = sum(SortedEquilibrium(2,:))/numel(find(BindingUnit == 6));  %Avg Equilibrium values for dimer only
%     yline(Mean_Equilibrium_D,'--k', ['\rho = 0 (', num2str(round(Mean_Equilibrium_D,3)), ')'],'LineWidth',1);
% end

figure(1);
% subplot(2,1,1);
Legend = cell(length(BindingUnit),1);
for c = 1:length(BindingUnit)
    Legend{c} = ['Binding Unit = ', num2str(BindingUnit(c))];
end
legend(Legend,'location','southeast');