function yrtix(axhand)
xt = get(gca, 'XTick');
Yu = unique(fix(xt));
xtix = linspace(min(Yu), max(Yu),  numel(Yu));
set(gca, 'XTick', xtix)
end
