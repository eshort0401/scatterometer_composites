% Remember to intput modified X Y when using m_pcolor.
function h=contourPlot(X,Y,z,axisMin,axisMax,step,tick,labels,div,pcolour)
    h=figure('units','centimeters','pos',[0 0 18 6]);
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + 1.5*ti(1);
    bottom = outerpos(2) + 1.5*ti(2);
    ax_width = outerpos(3) - 2.5*ti(1) - 2.5*ti(3);
    ax_height = outerpos(4) - 2*ti(2) - 2*ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    % Set color palette.
    if nargin<=8 || div==0
        cMap=(cbrewer('seq','Reds',length(axisMin:step:axisMax)-1,'pchip'));
    elseif div==1
        cMap=flipud(cbrewer('div','RdBu',length(axisMin:step:axisMax)-1,'pchip'));
    elseif div==2
        cMap=flipud(cbrewer('qual','Set1',length(axisMin:step:axisMax)-1,'pchip'));
    elseif div==3
        cMap=flipud(cbrewer('seq','Blues',length(axisMin:step:axisMax)-1,'pchip'));
    elseif div==4
        cMap=cbrewer('div','RdBu',length(axisMin:step:axisMax)-1,'pchip');
    elseif div==5
        cMap=flipud(cbrewer('seq','Greys',length(axisMin:step:axisMax),'pchip'));
    elseif div==6
        cMap=cbrewer('seq','Greys',length(axisMin:step:axisMax)-1,'pchip');
    end
        
    tMap=flipud(cbrewer('seq','Greys',5,'pchip'));
    
    if nargin<=9
        pcolour=0;
    end
    
    if pcolour
        X=X-.125;
        Y=Y-.125;
        m_pcolor(X,Y,z'); 
        shading flat;
    else
        [cc, ch]=m_contourf(X,Y,z',axisMin:step:axisMax);
%         clabel(cc,ch,'FontSize',12,'FontName','Times New Roman','LabelSpacing',2400);

        set(ch,'edgecolor','none');
        clear ch;
    end
    colormap(cMap);

    caxis([axisMin,axisMax]);
    
    if labels==true
        c=colorbar;
        set(c,'YTick',axisMin:step:axisMax,'FontSize',12,'FontName','Times New Roman');
        ylabel(c,'m/s','FontSize',12,'FontName','Times New Roman'); 
    else
        %c=colorbar;
    end
    
    m_gshhs_i('color','black','LineWidth',0.25);
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
