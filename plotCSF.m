function ahandle = plotCSF(spatFreqs,S)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spatFreqs  = frequencies to plot in cycles per visual degree
% S          = log(sensitivities) corresponding to spatFreqs
%
% JPO 2024, Hamburg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fh = figure;
plot(log2(spatFreqs),S,'o-','Color',[0 0 0],'LineWidth',2,'MarkerFaceColor',[1 0 0],'MarkerSize',6);                        % log in base 2 is usualy used in graphs
axis([log2(.1) log2(35) log10(1) log10(250)])       
xticks = [.5 1 2:2:10 14:4:30];
yticks = [1 5 10 50 100 200];
ahandle = gca;
set(ahandle,'XTick',log2(xticks),'XTickLabel',xticks,...
    'YTick',log10(yticks),'YTickLabel',yticks)
xlabel('Spatial frequency cyc/deg')
ylabel('Sensitivity')
title('CSF truncated log-parabola model')