clear all
no_of_nodes=4;
[x,v,s,h,f] = nafld_wild_type_bifur; 
a = x(5,:); %bifurcation parameter
%b = x(1,:); %x(1):HNF4A; x(2):HNF1A; x(3):PPARG; x(4):SREBF1
c = a;%./1000;

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

%%
figure1 = figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1]);

b = log2(x(2,:));

plot(c(1:ind(1)),b(1:ind(1)),'b');
hold on
plot(c(ind(1)+1:ind(2)),b(ind(1)+1:ind(2)),'r');
% plot(c(ind(2)+1:ind(3)),b(ind(2)+1: ind(3)),'b');
% plot(c(ind(3)+1:ind(4)),b(ind(3)+1:ind(4)),'r');
plot(c(ind(2)+1:end),b(ind(2)+1:end),'b');

xlim([0 3]);
%ylim([0 15]);
xlabel('PPARg degradation rates');
ylabel('PPARg levels');
