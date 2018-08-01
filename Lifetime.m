function Lifetime(stateSeq,SPL)
%--------------------------------------------------------------------------
% This function calculates the lifetime of different states 
% note: I implement this to include only two states, feel free to modify
% to include more states 
% Code written by:
%       Ben Guan 
%--------------------------------------------------------------------------
tic
state=[];


for i=1:length(stateSeq)
    states=stateSeq{i,1};
    for m=1:length(states)
        l=length(states);
        foo=diff(states); %find the difference of each state
        goo=[transpose(find(foo)) l]; % find index of transition
        goo1=diff(goo);% find spacing between transition
    end
    for k=1:length(goo1)
        index=goo(k);
        state=states(index);
        if state==1
            state{1}=[state1 goo1(k)];
            
        elseif state==2
            state{2}=[state2 goo1(k)];
        end
        
    end
end
figure;
bin=0:SPL:400;

histogram(state1,bin);
title('Histogram of lifetime for state 1','fontsize',20);

xlabel('Life time','fontsize',20);
ylabel('Number of States','fontsize',20);
figure;
histogram(state2,bin);
title('Histogram of lifetime for state2','fontsize',20);
xlabel('Life time','fontsize',20);
ylabel('Number of States','fontsize',20);

avg1=mean(state{1});
avg2=mean(state{2});
disp(['Avg lifetime for state 1 and state 2 are: ' num2str([avg1 avg2])]);
end