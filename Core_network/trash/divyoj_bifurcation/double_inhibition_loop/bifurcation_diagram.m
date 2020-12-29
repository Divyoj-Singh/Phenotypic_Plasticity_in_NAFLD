clear all
no_of_initial_conditions=1000;
plotting_array=double.empty();
index_no=0;

%% Running loop over different parameter:
for bifurcation_parameter= 0:0.1:5
%% Parameter Set Specification:
index_no=index_no+1;
bifurcation_parameter
plotting_array(index_no,1)=bifurcation_parameter;
%% Time Domain:
% mention the time domain for which u want to run the ODE. 
domain = [0 1000];
s_state=[];
%% Starting the loop for different inital conditions:
parfor j=1:no_of_initial_conditions
% picking random initial condition for the species:    
% here we picked a random number in the ra4nge of 1-20 and converted it to log2 scale.(as done in RACIPE) 
IA = 2^(20*rand(1));
IB  = 2^(20*rand(1));
%% Calling ODE function:
[t, x] = ode45(@(t,x) interactions(t,x,bifurcation_parameter),domain,[IA;IB]);

% now converting the values back to log 
x=log2(x);
end_point=round(x(end,1),1);
s_state=[s_state,end_point];
end
steady_state=[];
[a,b] = histc(s_state,unique(s_state));
y = a(b);
for i=1:length(s_state)
    if y(i)> 10
    steady_state=[steady_state s_state(i)];
    end
end
steady_state=unique(steady_state);

for i=1:length(steady_state)
    plotting_array(index_no,i+1)=steady_state(i);
end

end
%% Plotting the bifurcation Diagram
figure
plotting_array(plotting_array==0)=NaN;
plot(plotting_array(:,1), plotting_array(:,2:end)', '.')
yrtix(gca)