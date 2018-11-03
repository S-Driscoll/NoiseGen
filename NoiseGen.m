function [X,COV,COVM] = NoiseGen(chan,reps,varargin)
%%
%                       NoiseGen   v. 1.1
%
% This function simulates multiple types of instrument noise commonly
% observed in analytical chemistry using a 6-term model of the error covariance
% matrix. (See PAPER REFERENCE). The explicit inputs passed to the function are:
%
%   chan   -    the length of the output signal.
%
%   reps   -    the number of signals to generate.
%
% Additional variable length input arguments are:
%   'BO'  - Random baseline offset noise
%    [X, COV, COVM]=NoiseGen(chan,reps,'BO',SD);
%    SD (1 x 1)
%
%   'MO'  - Multiplicative offset noise
%    [X, COV, COVM]=NoiseGen(chan,reps,'MO',RSD,ReferenceSignal);
%    [X, COV, COVM]=NoiseGen(chan,reps,'MO',RSD,ReferenceSignal,Proportionality);
%    RSD (1 x 1) or (1 x reps)
%    ReferenceSignal (1 x chan) or (reps x chan)
%
%   'IID'  - White noise
%    [X, COV, COVM]=NoiseGen(500,400,'IID',SD);
%    SD (1 x 1)
%
%   'MN'  - Multiplicative noise
%    [X, COV, COVM]=NoiseGen(500,400,'MN',0.1,y,1);
%
%   'PL'  - Power law noise (1/f^alpha noise)
%    [X, COV, COVM]=NoiseGen(500,400,'PL',0.1,y,1);
%
%   'PPL' - Proportional power law noise (1/f^alpha proportional noise)
%    [X, COV, COVM]=NoiseGen(500,400,'PPL',0.1,y,1);
%
%
% The returned variables are:
%
%   X    -    the (reps x chan) matrix of simulated noise signals.
%
%   COV  -    the (chan x chan) covariance matrix for X (X'*X/(reps-1)).
%
%   COVM -    the (chan x chan) theoretical covariance matrix
%             If different ECM for each row, then it has dimension (chan x
%             chan x reps)
%
% Author: Steve Driscoll
% Trace Analysis Research Centre, Department of Chemistry,
% Dalhousie University, Halifax, Nova Scotia, Canada B3H 4J3
% Email: Stephen.Driscoll@dal.ca
% Website: http://groupwentzell.chemistry.dal.ca/
% April 2017; Last revision: Nov 2018
% V1.1 change log: - Added functionality to add multiple noise sources of
% the same type
%%

% Comb varargin to get parameter structure
Inpt=varargin; % Ressign input
% Check for any invalid paramter calls
s2={'MO','BO','IID','MN','PL','PPL'};
for i=1:length(Inpt)
    if iscellstr(Inpt(i))==1
        P_check=strcmpi(Inpt{i},s2);
        if isempty(find(P_check==1,1))==1
            tempT=['Unrecognized input string "', (Inpt{i}), '" (Possibilities are MO, BO, IID, MN, PL, and PPL)'];
            error(tempT)
        end
    end
end
% Build boolean vector of Inpt structure wrt valid parameters
cellfind=@(string)(@(cell_contents)(strcmpi(string,cell_contents))); % Inline function for strcmpi --> boolean
% Define cypher for paramater list [1=MO 2=BO 3=IID 4=MN 5=PN 6=PPN]
BooPa=cellfun(cellfind('MO'),Inpt);
BooPa=BooPa+2*cellfun(cellfind('BO'),Inpt);
BooPa=BooPa+3*cellfun(cellfind('IID'),Inpt);
BooPa=BooPa+4*cellfun(cellfind('MN'),Inpt);
BooPa=BooPa+5*cellfun(cellfind('PL'),Inpt);
BooPa=BooPa+6*cellfun(cellfind('PPL'),Inpt);
% Now search boolean vector for model paramters according to cypher
Pos=find(BooPa); % Find where noise types are in input
PosID=BooPa(Pos); % Find what noise types are in input
k=2; % Counter for # of params associated with each noise type after param string (2 because (1) is the noise type string)
ECM_count=1; % Counter for # of ECMs to be summed
ECM_MO_Row_Flag=0; % Flag for multiplicative offset noise per row
ECM_MN_Row_Flag=0; % Flag for multiplicative noise per row
ECM_PN_Row_Flag=0; % Flag for power law noise per row
Mcount_BO=1; % Counter for multiple BO offset contributions
Mcount_MO=1; % Counter for multiple MO offset contributions
Mcount_PN=1; % Counter for multiple PL offset contributions
Mcount_PPN=1; % Counter for multiple PPL offset contributions
Mcount_IID=1; % Counter for multiple IID offset contributions
Mcount_MN=1; % Counter for multiple MN offset contributions
ECM=zeros(chan,chan,1); % Pre-allocating
for i=1:length(BooPa)
    if BooPa(i)~=0 % a 0 in BooPa means we are not interested in element value
        if BooPa(i)==1 % MO noise begin
            RefSig=zeros(1,chan); % Force default
            c=1; % Counter for hopping along parameters (variable length)
            if length(Pos)==1
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==1 && numel(find(PosID==1))>1 && Mcount_MO<numel(find(PosID==1))
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
                Mcount_MO=Mcount_MO+1;
            elseif PosID(end)==1 && numel(find(PosID==1))>1 && Mcount_MO==numel(find(PosID==1))
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==1
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            else
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            end
            if length(P1{1})==1 && (isrow(P1{2})==1 || iscolumn(P1{2})==1) % One RSD, One ref vector
                if length(P1)==3
                    Contr=P1{1}; % Contribution
                    RefSig=P1{2}; % Reference signal
                    Prop=P1{3}; % Proprtionality
                elseif length(P1)==2
                    Contr=P1{1}; % Contribution
                    RefSig=P1{2}; % Reference signal
                    Prop=1; % Default
                else
                    error('MO noise requires a reference signal (1 x chan)')
                end
                if isrow(RefSig)==1
                    RefSig=RefSig';
                end
                if length(RefSig)~=chan
                    error('Length of the reference signal should be equal to chan (1 x chan)')
                end
                ECM(:,:,ECM_count)=[RefSig.^(2*(Prop))]*[Contr^2]*[RefSig.^(2*(Prop))]'; % Build ECM
            elseif length(P1{1})==1 && isrow(P1{2})==0 && iscolumn(P1{2})==0 % MO prop each row single prop
                ECM_MO_Row_Flag=1;
                if length(P1)==3
                    Contr=P1{1}; % Contribution
                    RefSig=P1{2}; % Reference signal
                    Prop=P1{3}; % Proprtionality
                elseif length(P1)==2
                    Contr=P1{1}; % Contribution
                    RefSig=P1{2}; % Reference signal
                    Prop=1; % Default
                else
                    error('MO noise requires a reference signal (1 x chan)')
                end
                [rowe, cole]=size(RefSig);
                if rowe~=reps && cole~=chan
                    error('Expecting reference signal matrix either 1 x chan or m x chan')
                end
                if length(RefSig)~=chan
                    error('Length of the reference signal should be equal to chan (1 by chan)')
                end
                for RR=1:reps
                    ECM_Row_MO(:,:,RR)=[RefSig(RR,:).^(2*(Prop))]'*[Contr^2]*[RefSig(RR,:).^(2*(Prop))]; % Build ECM
                end
            elseif length(P1{1})>1  % MO prop each row different prop
                ECM_MO_Row_Flag=1;
                if isrow(P1{2})==1 || iscolumn(P1{2})==1
                    error('Expecting reference signal matrix of size reps by chan')
                end
                if length(P1)==3
                    Contr=P1{1}; % Contribution
                    if length(Contr)~=reps
                        error('Length of contributions should be equal to reps')
                    end
                    RefSig=P1{2}; % Reference signal (stacked m by n)
                    if length(RefSig)~=chan
                        error('Length of reference signal should be equal to chan')
                    end
                    Prop=P1{3}; % Proprtionality
                elseif length(P1)==2
                    Contr=P1{1}; % Contribution
                    if length(Contr)~=reps
                        error('Length of contributions should be equal to reps')
                    end
                    RefSig=P1{2}; % Reference signal
                    Prop=1; % Default
                else
                    error('Vector of contributions detected (reps by 1), expecting reference signal matrix reps by chan')
                end
                for RR=1:reps
                    ECM_Row_MO(:,:,RR)=[RefSig(RR,:).^(2*(Prop))]'*[Contr(RR)^2]*[RefSig(RR,:).^(2*(Prop))];% Build ECM
                end
            end
            k=k+1;
            ECM_count=ECM_count+1;
        end
        if BooPa(i)==2 % BO noise start
            c=1;
            if length(Pos)==1
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==2 && numel(find(PosID==2))>1 && Mcount_BO<numel(find(PosID==2))
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
                Mcount_BO=Mcount_BO+1;
            elseif PosID(end)==2 && numel(find(PosID==2))>1 && Mcount_BO==numel(find(PosID==2))
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==2
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            else
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            end
            if length(P1{1})==1
                Contr=P1{1};
                ECM(:,:,ECM_count)=[ones(chan,1)]*[Contr^2]*[ones(chan,1)]';
            else
                error('Expecting scalar for BO contribution')
            end
            k=k+1;
            ECM_count=ECM_count+1;
        end
        if BooPa(i)==3 % IID
            c=1;
            if length(Pos)==1
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==3 && numel(find(PosID==3))>1 && Mcount_IID<numel(find(PosID==3))
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
                Mcount_IID=Mcount_IID+1;
            elseif PosID(end)==3 && numel(find(PosID==3))>1 && Mcount_IID==numel(find(PosID==3))
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==3
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            else
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            end
            if length(P1{1})==1
                Contr=P1{1};
                ECM(:,:,ECM_count)=Contr^2*diag(ones(1,chan));
            else
                error('Expecting scalar for IID contribution')
            end
            k=k+1;
            ECM_count=ECM_count+1;
        end
        if BooPa(i)==4 % MN
            RefSig=zeros(1,chan); % Force default
            c=1; % Counter for hopping along parameters
            if length(Pos)==1
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==4 && numel(find(PosID==4))>1 && Mcount_MN<numel(find(PosID==4))
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
                Mcount_MN=Mcount_MN+1;
            elseif PosID(end)==4 && numel(find(PosID==4))>1 && Mcount_MN==numel(find(PosID==4))
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==4
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            else
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            end
            if length(P1{1})==1 && (isrow(P1{2})==1 || iscolumn(P1{2})==1) % One RSD, One ref vector
                if length(P1)==3
                    Contr=P1{1}; % Contribution
                    RefSig=P1{2}; % Reference signal
                    Prop=P1{3}; % Proprtionality
                elseif length(P1)==2
                    Contr=P1{1}; % Contribution
                    RefSig=P1{2}; % Reference signal
                    Prop=1; % Default
                else
                    error('MN noise requires a reference signal (1 x chan)')
                end
                if isrow(RefSig)==1
                    RefSig=RefSig';
                end
                if length(RefSig)~=chan
                    error('Length of the reference signal should be equal to chan (1 by chan)')
                end
                ECM(:,:,ECM_count)=Contr^2*diag(RefSig.^(2*(Prop))); % Build ECM
            elseif length(P1{1})==1 && isrow(P1{2})==0 && iscolumn(P1{2})==0 % MN prop each row single prop
                ECM_MN_Row_Flag=1;
                if length(P1)==3
                    Contr=P1{1}; % Contribution
                    RefSig=P1{2}; % Reference signal
                    Prop=P1{3}; % Proprtionality
                elseif length(P1)==2
                    Contr=P1{1}; % Contribution
                    RefSig=P1{2}; % Reference signal
                    Prop=1; % Default
                else
                    error('MN noise requires a reference signal (1 x chan)')
                end
                [rowe, cole]=size(RefSig);
                if rowe~=reps && cole~=chan
                    error('Expecting reference signal matrix either 1 by chan or m by chan')
                end
                if length(RefSig)~=chan
                    error('Length of the reference signal should be equal to chan (1 by chan)')
                end
                for RR=1:reps
                    ECM_Row_MN(:,:,RR)=Contr^2*diag(RefSig(RR,:).^(2*(Prop))); % Build ECM
                end
            elseif length(P1{1})>1  % MN prop each row different prop
                ECM_MN_Row_Flag=1;
                
                if isrow(P1{2})==1 || iscolumn(P1{2})==1
                    error('Expecting reference signal matrix of size reps by chan')
                end
                if length(P1)==3
                    Contr=P1{1}; % Contribution
                    if length(Contr)~=reps
                        error('Length of contributions should be equal to reps')
                    end
                    RefSig=P1{2}; % Reference signal (stacked m by n)
                    if length(RefSig)~=chan
                        error('Length of reference signal should be equal to chan')
                    end
                    Prop=P1{3}; % Proprtionality
                elseif length(P1)==2
                    Contr=P1{1}; % Contribution
                    if length(Contr)~=reps
                        error('Length of contributions should be equal to reps')
                    end
                    RefSig=P1{2}; % Reference signal
                    Prop=1; % Default
                else
                    error('Vector of contributions detected (reps by 1), expecting reference signal matrix reps by chan')
                end
                for RR=1:reps
                    ECM_Row_MN(:,:,RR)=Contr(RR)^2*diag(RefSig(RR,:).^(2*(Prop))); % Build ECM
                end
            end
            k=k+1;
            ECM_count=ECM_count+1;
        end
        if BooPa(i)==5 % PN
            c=1; % Counter for hopping along parameters
            if length(Pos)==1
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==5 && numel(find(PosID==5))>1 && Mcount_PN<numel(find(PosID==5))
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
                Mcount_PN=Mcount_PN+1;
            elseif PosID(end)==5 && numel(find(PosID==5))>1 && Mcount_PN==numel(find(PosID==5))
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==5
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            else
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            end
            if length(P1)==3
                Contr=P1{1}; % Contribution
                Alpha=P1{2}; % Alpha
                Rho=P1{3}; % Rho
            elseif length(P1)==2
                Contr=P1{1};
                Alpha=P1{2};
                Rho=1;
            elseif length(P1)==1
                Contr=P1{1};
                Alpha=1;
                Rho=1;
            end
            ECM(:,:,ECM_count)=fgen(chan,Alpha,Contr,Rho);
            k=k+1;
            ECM_count=ECM_count+1;
        end
        if BooPa(i)==6 % PPN
            c=1; % Counter for hopping along parameters
            if length(Pos)==1
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==6 && numel(find(PosID==6))>1 && Mcount_PPN<numel(find(PosID==6))
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
                Mcount_PPN=Mcount_PPN+1;
            elseif PosID(end)==6 && numel(find(PosID==6))>1 && Mcount_PPN==numel(find(PosID==6))
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            elseif PosID(end)==6
                for j=i+1:length(BooPa)
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            else
                for j=i+1:Pos(k)-1
                    P1(c)=Inpt(j);
                    c=c+1;
                end
            end
            if length(P1{1})==1 && (isrow(P1{2})==1 || iscolumn(P1{2})==1) % One RSD, One ref vector
                if length(P1)==4
                    Contr=P1{1};
                    RefSig=P1{2};
                    Alpha=P1{3};
                    Rho=P1{4};
                elseif length(P1)==3
                    Contr=P1{1};
                    RefSig=P1{2};
                    Alpha=P1{3};
                    Rho=1;
                elseif length(P1)==2
                    Contr=P1{1};
                    RefSig=P1{2};
                    Alpha=1;
                    Rho=1;
                end
                if isrow(RefSig)==1
                    RefSig=RefSig';
                end
                if length(RefSig)~=chan
                    error('Length of the reference signal should be equal to chan (1 by chan)')
                end
                ECM=fgen(chan,Alpha,1,Rho);
                fcovp=ECM.*(RefSig*RefSig');
                ECM(:,:,ECM_count)=Contr^2*fcovp; % Build ECM
            elseif length(P1{1})==1 && isrow(P1{2})==0 && iscolumn(P1{2})==0 % PN prop each row single prop
                ECM_PN_Row_Flag=1;
                if length(P1)==4
                    Contr=P1{1};
                    RefSig=P1{2};
                    Alpha=P1{3};
                    Rho=P1{4};
                elseif length(P1)==3
                    Contr=P1{1};
                    RefSig=P1{2};
                    Alpha=P1{3};
                    Rho=1;
                elseif length(P1)==2
                    Contr=P1{1};
                    RefSig=P1{2};
                    Alpha=1;
                    Rho=1;
                end
                [rowe, cole]=size(RefSig);
                if rowe~=reps && cole~=chan
                    error('Expecting reference signal matrix either 1 by chan or m by chan')
                end
                if length(RefSig)~=chan
                    error('Length of the reference signal should be equal to chan (1 by chan)')
                end
                for RR=1:reps
                    ECM=fgen(chan,Alpha,Contr,Rho);
                    fcovp=ECM.*(RefSig(RR,:)'*RefSig(RR,:));
                    ECM_Row_PPN(:,:,RR)=fcovp; % Build ECM
                end
                
            elseif length(P1{1})>1  % PN prop each row different prop
                ECM_PN_Row_Flag=1;
                
                if isrow(P1{2})==1 || iscolumn(P1{2})==1
                    error('Expecting reference signal matrix of size reps by chan')
                end
                if length(P1)==4
                    Contr=P1{1}; % Contribution
                    Alpha=P1{3};
                    Rho=P1{4};
                    if length(Contr)~=reps
                        error('Length of contributions should be equal to reps')
                    end
                    RefSig=P1{2}; % Reference signal (stacked m by n)
                    if length(RefSig)~=chan
                        error('Length of reference signal should be equal to chan')
                    end
                elseif length(P1)==3
                    Contr=P1{1}; % Contribution
                    Alpha=P1{3};
                    Rho=1;
                    if length(Contr)~=reps
                        error('Length of contributions should be equal to reps')
                    end
                    RefSig=P1{2}; % Reference signal
                elseif length(P1)==2
                    Contr=P1{1}; % Contribution
                    Alpha=1;
                    Rho=1;
                    if length(Contr)~=reps
                        error('Length of contributions should be equal to reps')
                    end
                    RefSig=P1{2}; % Reference signal
                else
                    error('Vector of contributions detected (reps by 1), expecting reference signal matrix reps by chan')
                end
                for RR=1:reps
                    ECM=fgen(chan,Alpha,Contr(RR),Rho);
                    fcovp=ECM.*(RefSig(RR,:)'*RefSig(RR,:));
                    ECM_Row_PPN(:,:,RR)=fcovp; % Build ECM
                end
            end
            k=k+1;
            ECM_count=ECM_count+1;
        end
    end
end
% All done building ECMs, now check if we should generate per row
if ECM_MO_Row_Flag>0 || ECM_MN_Row_Flag>0 || ECM_PN_Row_Flag>0
    COVMt=sum(ECM,3);
    for i=1:reps
        if ECM_MO_Row_Flag>0
            COVM_Ri=ECM_Row_MO(:,:,i);
        else
            COVM_Ri=zeros(chan,chan);
        end
        
        if ECM_MN_Row_Flag>0
            COVM_Ri2=ECM_Row_MN(:,:,i);
        else
            COVM_Ri2=zeros(chan,chan);
        end
        
        if ECM_PN_Row_Flag>0
            COVM_Ri3=ECM_Row_PN(:,:,i);
        else
            COVM_Ri3=zeros(chan,chan);
        end
        COVM(:,:,i)=COVM_Ri + COVM_Ri2 + COVM_Ri3 + COVMt;
        [~,S,V]=svd(COVM(:,:,i));
        eiid=randn(1,chan);
        esig=eiid*sqrt(S)*V';
        X(i,:)=esig;
        
    end
    COV=(X'*X)/(reps-1);
else % if not, then sum the ECMs and generate noise
    COVM=sum(ECM,3);
    [~,S,V]=svd(COVM);
    eiid=randn(reps,chan);
    esig=eiid*sqrt(S)*V';
    COV=(esig'*esig)/(reps-1);
    X=esig;
end
end

%% Function to generate 1/f^\alpha ECM
function [ECMout]=fgen(chan,a,p,rho)
w=round(chan/2);
m=rho*2*w+1; % Filter width
f=[1/(3*m) 1/(2*w):1/(2*w):0.5]; % Frequency vector
H=abs(1./(f.^(a/2)));
H=[H fliplr(H(2:end))]; % Mask
Z=fft(H);
c=sqrt(real(Z).^2+imag(Z).^2); % Coefficients
c=[c(w+1:end) c(1:w)];
win=0.54-0.46*cos(2*pi*[1:length(c)]/length(c)); % Hamming window
c=c.*win;
c=c/norm(c); % Normalize

% ECM via filter matrix
npts=chan+2*w;
F=zeros(npts,chan);
for i=1:chan
    indx=i+2*w;
    F(i:indx,i)=c';
end
ECMout=p^2*(F'*F); % Construct ECM with filter

% % ECM via self convolution of c (might be slightly slower, but don't have to store F)
% ECMout=zeros(chan,chan);
% sc2=conv(c,c);
% for i = 1:chan
%     k=1;
%     for j = i:chan
%         ECMout(i,j)=sc2(k+m);
%         k=k+1;
%     end
% end
% ECMout=triu(ECMout)+triu(ECMout,1)';
end