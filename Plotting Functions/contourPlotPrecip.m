% Remember to intput modified X Y when using m_pcolor.
function h=contourPlotPrecip(X,Y,z,axisMin,axisMax,step,tick,labels,cScheme,conType)
    h=figure('units','centimeters','pos',[0 0 18 6]);
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + 1.5*ti(1);
    bottom = outerpos(2) + 1.5*ti(2);
    ax_width = outerpos(3) - 2.5*ti(1) - 2.5*ti(3);
    ax_height = outerpos(4) - 2*ti(2) - 2*ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    steps=.5.^(2:-1:-2);
    
    % Set color palette.
    if nargin<=8 || cScheme==0
        cMap=(cbrewer('seq','Reds',length(axisMin:step:axisMax)-1,'pchip'));
    elseif cScheme==1
        cMap=flipud(cbrewer('div','RdBu',length(axisMin:step:axisMax)-1,'pchip'));
    elseif cScheme==2
        cMap=flipud(cbrewer('qual','Set1',length(axisMin:step:axisMax)-1,'pchip'));
    elseif cScheme==3
        cMap=(cbrewer('seq','Blues',length(axisMin:step:axisMax)-1,'pchip'));
    elseif cScheme==4
        cMap=cbrewer('div','RdBu',length(axisMin:step:axisMax)-1,'pchip');
    elseif cScheme==5
        cMap=cbrewer('seq','Blues',length(steps)+2,'pchip');
        cMap=.95*cMap(3:end,:);
    end
    tMap=flipud(cbrewer('seq','Greys',5,'pchip'));
    
    if nargin<=9
        conType=0;
    end
    
    if conType==1
        X=X-.125;
        Y=Y-.125;
        m_pcolor(X,Y,z'); 
        shading flat;
    elseif conType==2
        [C, h]=m_contour(X,Y,log2(z'),log2(steps),'LineWidth',1.5);
    else
        [~, ch]=m_contourf(X,Y,log2(z'),log2(steps)); 
        set(ch,'edgecolor','none');
        clear ch;
    end
    colormap(cMap);
    c=colorbar;
    caxis([log2(steps(1))-0.5, log2(steps(end))+0.5]);
    set(c,'YTick',log2(steps),'YTickLabel',steps,'FontSize',12,'FontName','Times New Roman')
    ylabel(c,'mm/h','FontSize',12,'FontName','Times New Roman');
    
    if labels==true
        c=colorbar;
        caxis([log2(steps(1))-0.5, log2(steps(end))+0.5]);
        set(c,'YTick',log2(steps),'YTickLabel',steps,'FontSize',12,'FontName','Times New Roman')
        ylabel(c,'m/s','FontSize',12,'FontName','Times New Roman'); 
    else
        %c=colorbar;
    end
    
%     m_gshhs_i('color','black','LineWidth',0.25);
    % Plot the grid
    if labels==true
        m_grid('xtick',0:tick:360,'ytick',-90:tick:90,'FontSize',12,'FontName','Times New Roman','glinewidth',.5);
    else
        m_grid('xtick',0:tick:360,'ytick',-90:tick:90,'xticklabels',[],'yticklabels',[],'glinewidth',.5);
    end
    % Plot topography contour, note MC elevation is 4509m in New Guinea
%     m_tbase('contour',1000:1000:2000,'color',tMap(1,:),'LineWidth',0.25);
%     m_tbase('contour',2000:1000:3000,'color',tMap(2,:),'LineWidth',0.25);
%     m_tbase('contour',3000:1000:4000,'color',tMap(3,:),'LineWidth',0.25);
%     m_tbase('contour',4000:1000:5000,'color',tMap(4,:),'LineWidth',0.25);

end
