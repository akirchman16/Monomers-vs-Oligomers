clearvars;
% close Figure 1;

% This code will model reactions between free proteins and
% a one-dimensional DNA lattice. Dynamical Monte Carlo method (Gillespie
% Algorithm - First Reaction Method) will be used to model the individual reactions between the
% proteins and DNA. The fractional coverage of the DNA lattice will be
% plotted over time until equilibirum is reached.

N = 8660;    %length of DNA lattice
n = 3;  %length of each protein
k_on = 1;   %kinetic rate constant for binding (units: 1/(nt*s)??? )
k_off = 1;  %kinetic rate constant for unbinding (units: 1/s??? )
L = 2;    %concentration of free proteins (units: Molarity???)

minIterations = 1000;  %minimum number of evnts that will occur
K = k_on/k_off; %equilibrium constant (units: )

Cooperativity = [1,1,1]; %cooperativity values for the simulation
Loops = 0;

for w = Cooperativity
    Loops = Loops+1;

    DNA = zeros(1,N+2);   %empty DNA lattice with dummy zeros
    BoundAtSpot = zeros(1,N+2);  %empty CurrentBound array with dummy zeros

    t = zeros(1,minIterations+1);
    a = zeros(minIterations+1,4);
    dt = zeros(1,minIterations);
    FracCover = zeros(1,minIterations+1);
    BindHist = zeros(4,minIterations); %1: isolated binding, 2: SC binding, 3: DC binding, 4: unbinding
    xB_I = zeros(1,minIterations);
    xB_SC = zeros(1,minIterations);
    xB_DC = zeros(1,minIterations);
    xAB = zeros(1,minIterations);
    Length_nm = zeros(1,minIterations);

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

    while ~Equilibrium || ~EqualTimes
%         Isolated = 0;   %preparing arrays for types of available locations
%         SinglyContiguous = 0;
%         DoublyContiguous = 0;
%         Counter_I = 0;  %resets counters to count types of locations
%         Counter_SC = 0;
%         Counter_DC = 0;
%         for x = 2:N-(n-1)+1
%             if DNA(x:x+(n-1)) == 0
%                 if DNA(x-1) == 0 && DNA(x+n) == 0   %records all isolated locations
%                     Isolated(Counter_I+1) = x;
%                     Counter_I = Counter_I+1;
%                 elseif (DNA(x-1) == 0 && DNA(x+n) == 1) || (DNA(x-1) == 1 && DNA(x+n) == 0) %records all singly contiguous locations
%                     SinglyContiguous(Counter_SC+1) = x;
%                     Counter_SC = Counter_SC+1;
%                 elseif DNA(x-1) == 1 && DNA(x+n) == 1   %records all doubly contiguous locations
%                     DoublyContiguous(Counter_DC+1) = x;
%                     Counter_DC = Counter_DC+1;
%                 end
%             end
%         end
%         xB_I(Events+1) = Counter_I;    %amounts of each location, used in propensity functions
%         xB_SC(Events+1) = Counter_SC;
%         xB_DC(Events+1) = Counter_DC;
%         xAB(Events+1) = sum(DNA)/n;

%   NEW SEARCH METHOD TO COMPARE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Left = find(diff([1 DNA(2:N+1) 1]) == -1)+1;
        Right = find(diff([1 DNA(2:N+1) 1]) == 1);

        Gaps_L = Left(Right-Left+1 >= n);
        Gaps_R = Right(Right-Left+1 >= n);
        Gap_Size = Gaps_R-Gaps_L+1;

        % Doubly Contiguous Search
        DoublyContiguous = Left(Right-Left+1 == n & Left>2 & Right<N+1-(n-1));
        % End Doubly Contiguous Search

        Gaps_L2 = Gaps_L(~ismember(Gaps_L,DoublyContiguous));
        Gaps_R2 = Gaps_R(~ismember(Gaps_R,DoublyContiguous+(n-1)));
        Gap_Size2 = Gaps_R2-Gaps_L2+1;

        % Singly Contiguous Search
        SinglyContiguous = unique([Gaps_L2(Gaps_L2>2),Gaps_R2(Gaps_R2<N+1-(n-1))-(n-1)]);
        % End Singly Contiguous Search

        Gap_SizeI = Gap_Size2(Gap_Size2 > n+1);
        Gaps_LI = Gaps_L2(ismember(Gap_Size2,Gap_SizeI));
        Real_Gap_SizeI = [Gap_SizeI(1:end-1),Gap_SizeI(end)+(numel(find(Gaps_R2(ismember(Gap_Size2,Gap_SizeI)) == N+1)))];

        % Isolated Search
        Isolated = zeros(1,sum(Real_Gap_SizeI)-((n+1)*(numel(Gap_SizeI)))+logical(sum(DNA) == 0));
        for i = 1:numel(Gaps_LI)
            Isolated(find(Isolated == 0, 1 ):find(Isolated == 0, 1)+Real_Gap_SizeI(i)-(n+1)-1+logical(sum(DNA) == 0)) = Gaps_LI(i)+(1:Real_Gap_SizeI(i)-(n+1)+logical(sum(DNA) == 0));
        end
        Isolated(Isolated == 0) = [];
        % End Isolated Search
        
        % Population Numbers
        xB_I(Events+1) = length(Isolated);
        xB_SC(Events+1) = length(SinglyContiguous);
        xB_DC(Events+1) = length(DoublyContiguous);
        xAB(Events+1) = sum(DNA)/n;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        elseif j(Events+1) == 2    %tests for singly contiguous binding event
            pos_SC = randi(length(SinglyContiguous));
            SpotB_SC = SinglyContiguous(pos_SC);  %chooses a random location for binding
            DNA(SpotB_SC:SpotB_SC+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_SC = BindCounter_SC+1;
            BindHist(2,Events+1) = SpotB_SC;
            BoundAtSpot(SpotB_SC) = 1;    %shows there's currently a protein bound at location
        elseif j(Events+1) == 3    %tests for doubly contiguous binding event
            pos_DC = randi(length(DoublyContiguous));
            SpotB_DC = DoublyContiguous(pos_DC);  %chooses a random location for binding
            DNA(SpotB_DC:SpotB_DC+(n-1)) = 1; %binds protein
            BindCounter = BindCounter+1;
            BindCounter_DC = BindCounter_DC+1;
            BindHist(3,Events+1) = SpotB_DC;
            BoundAtSpot(SpotB_DC) = 1;    %shows there's currently a protein bound at location
        elseif j(Events+1) == 4    %otherwise an unbinding event occurs
            CurrentBound = find(BoundAtSpot == 1);
            pos = randi(length(CurrentBound));
            SpotU = CurrentBound(pos);  %selects an already bound protein to unbind

            DNA(SpotU:SpotU+(n-1)) = 0; %unbinds protein
            BoundAtSpot(SpotU) = 0;
            UnbindCounter = UnbindCounter+1;
            BindHist(4,Events+1) = SpotU;
        end

        FracCover((Events+1)+1) = sum(DNA)/N;  %fractional coverage

        if Events > minIterations && Equilibrium == 0
            FracCoverStates = (FracCover((Events-floor(0.25*Events)):1:Events));   %every 10 states for the last 20% of the model
            FracCoverChange = abs(diff(FracCoverStates));   %difference among each of those states
            if (mean(FracCoverChange) <= 2*n/N) && (abs(FracCover(Events-floor(0.25*Events))-FracCover(Events)) <= 3*n/N)
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
        Avg_FilamentLength(Events) = mean(FilamentLengths(Events,1:FilamentCount)); %average filament length for each event
    end

    Length_um = Length_nm/1e+3; %length of molecule in um
    Equilibrium_Coverage(Loops) = mean(FracCoverStates);
    disp(['(w = ', num2str(w),') - Equilibrium Saturation = ', num2str(Equilibrium_Coverage(Loops))]);

    TotalTime(Loops) = max(t);
    LoopLengths = sort(TotalTime,'descend');
    
    EventFractions(Loops,:) = [length(find(j==1)),length(find(j==2)),length(find(j==3)),length(find(j==4))]./Events;

    figure(1);
    
%     Plot of fractional coverage of the DNA
    
%     subplot(2,1,1);
    scatter(t,FracCover,2,'b','filled');
    hold on;
    % plot(X,Y_exp_fit,'b');  %bounded exponential fit
    % plot(X,Y_sig_fit,'k');  %sigmoidal curve fit
    xlabel('Time, t');
    xlim([0 max(LoopLengths)]);
    ylabel('Fractional Coverage');
    ylim([0 1]);
%     yline(Equilibrium_Coverage(Loops),'k',['\omega = ', num2str(w)]);
    title('Saturation of DNA Lattice');

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

figure(1);
% subplot(2,1,1);
Legend = cell(length(Cooperativity),1);
for c = 1:length(Cooperativity)
    Legend{c} = ['\omega = ', num2str(Cooperativity(c))];
end
legend(Legend,'location','southeast');