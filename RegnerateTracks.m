
%--------------------------------------------------------------------------
% This file takes in delta X of different states, generates a markov chain, 
% chooses deltaX based on the state sequence, calculate the accumulative sum
% to get the actual x, y position of tracks. 

% Code written by:
%       Ben Guan 
%--------------------------------------------------------------------------

file1=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/725/dX1.txt');
file2=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/725/dX2.txt');
file3=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/725/dX3.txt');
file4=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/725/dX4.txt');
file5=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/725/dX5.txt');
file6=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/725/dX6.txt');
file7=('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/725/dX7.txt');
file = {file1,file2, file3, file4, file5, file6,file7};


%%
formatspec='%f';
AllX=[];
AllDeltaXRAW=[];

for j=1:length(file) %loop through each file, change from 1:14 to 1:10?
    fid=fopen(file{j},'r');
    data=textscan(fid,formatspec); %cell containing each column from the text file
    DeltaX1=data(1);
    DX=DeltaX1{1}(3:3000);
    X=repmat(DX,10,1);
    AllDeltaXRAW=[AllDeltaXRAW, {X}];
end
%%
% split Delta X into unit of 10
for j=1:length(AllDeltaXRAW)
    DeltaX=AllDeltaXRAW(j);
    [X,splitIndex] = SplitTracks(DeltaX,10);
    AllDeltaX{j}=X;
end
% [XSP,splitIndex] = SplitTracks(AllDeltaXRAW,splitLength);
% est_stateSeq=[ones(1000,1);ones(1000,1)*2;ones(1000,1)*3;ones(1000,1)*4;ones(1000,1)*5];
% [MSD,Vacf]=MSD_Vacf(XSP,est_stateSeq,1,length(file));

%%
%set up Parameters
splitLength=10;

%population fraction
Pindex=[0.2 0.2 0.06 0.2 0.174 0.123 0.03];

%transition matrix

% p=0.04;
% transmat=[1-4*p p p p p;
%     p 1-4*p p p p;
%     p p 1-4*p p p;
%     p  p p 1-4*p p;
%     p p p p 1-4*p];

transmat=[0.9948 0.0023 0.00047 0.00057 0.000028 0.0009 0.0012;
    0.0062 0.9746 0.000438 0.007 0.0017 0.0018 0.0084;
    0.0011 0.0035 0.9760 0.0002 0.00006 0.0005 0.0185;
    0.0005 0.009 0.0006 0.98 0.003 0.0023 0.0046;
    0.0021 0.002 0.00013 0.01 0.9777 0.0057 0.0028;
    0.0009 0.0031 0.0026 0.007 0.0077 0.977 0.0013;
    0.0058 0.0576 0.018 0.0488 0.021 0.0018 0.846];

numTracks=200;
N = 100;  % length of particle tracks (length of movie)
%%
%Simulation markov chain
numStates=7;
markovStateSeq = SimulateMarkovState(numTracks,N,Pindex,transmat);
l=length(transmat);
real_stateSeq=cat(1,markovStateSeq{:});
% real_stateSeq=real_stateSeq(1:splitLength:end);

%%

%pick delta X from each states and calculate the culmulatice sum aseemble them into tracks
DXraw=[];
for k=1:length(markovStateSeq)
    stateseq=markovStateSeq{k};
    DXraw=[];
    
    for g=1:length(stateseq)
        state=stateseq(g);
        pool=AllDeltaX(state); %pick a pool of dx based on the state sequence 
        finalDX=pool{1,1}(10*(k-1)+g); %choose dx 
        DXDY=[finalDX{1} finalDX{1}];
        DXraw=[DXraw;DXDY];
        
    end
    Xacc=cumsum(DXraw,1);
    foo{k}=Xacc;
end

%%
[Xsplit,splitIndex] = SplitTracks(foo,splitLength); %feed in x y positions for each track, splitlength,
real_stateSeq=cat(1,markovStateSeq{:});

%%
[MSD,Vacf]=MSD_Vacf(Xsplit,real_stateSeq,1,length(file));

%%
savename='Bessel_Re_Tracks_uniform';
save(strcat('/Users/Jguan/Desktop/Yale/pEMv3/simulations_new/data/718/',savename),'foo','transmat','MSD','Vacf','real_stateSeq');
disp('tracks info saved');

