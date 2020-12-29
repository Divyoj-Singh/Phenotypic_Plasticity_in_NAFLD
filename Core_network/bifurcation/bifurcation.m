% This code plots the bifurcation diagram. It needs the MATCONT module.
% init.m from the MATCONT module should be run before running this code.

clear all
% Specify the number of nodes here:
no_of_nodes=4;
%% Calling the function interactions which in turn calls equations.
[x,v,s,h,f] = interactions; 
% order of genes: 1. HNF4A
%		  2. HNF1A
%	          3. PPARG
%                 4. SREBF1		
a = x(5,:); % bifurcation parameter(should be set to no_of_nodes+1)
b = x(3,:); % gene which you want to plot against the bifurcation parameter.
c = a;

%% Based on eigenvalues to judge stable vs. unstable states
ind = zeros(1,4);
snum = size(f);
num = snum(2);
j = 1;
n = 1;

for n = 1:1:(num-1)
    x1 = find(f(:,n) > 0);
    x2 = find(f(:,n+1) > 0);
    if isempty(x1) && ~isempty(x2)
        ind(j) = n + 1;
        j = j + 1;
    elseif ~isempty(x1) && isempty(x2)
        ind(j) = n + 1;
        j = j + 1;
    end
end

%% Plotting:
figure1 = figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);
% plotting the stable vs unstable in different colors red and blue:
plot(c(1:ind(1)),b(1:ind(1)),'b');
hold on
plot(c(ind(1)+1:ind(2)),b(ind(1)+1:ind(2)),'r');
plot(c(ind(2)+1:end),b(ind(2)+1:end),'b');

xlim([0 3]);
xlabel('Input signal to PPARg');
% mention the gene which is being plotted:
ylabel('PPARg levels');
