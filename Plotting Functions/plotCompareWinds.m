function plotCompareWinds(X,Y,u1,v1,u2,v2,axisMin,step,axisMax,windSpeed,quiverSubSample,tick)
    figure('units','centimeters','pos',[0 0 16 8]);
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + 1.5*ti(1);
    bottom = outerpos(2) + 1.5*ti(2);
    ax_width = outerpos(3) - 2.5*ti(1) - 2.5*ti(3);
    ax_height = outerpos(4) - 2*ti(2) - 2*ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    cMap=flipud(cbrewer('div','RdBu',length(axisMin:step:axisMax)-1,'pchip'));
    tMap=flipud(cbrewer('seq','Greys',5,'pchip'));
    
    dLon=X(1,2)-X(1,1);
    arrowScale=(8/8)*quiverSubSample*dLon;
%     arrowScale=1.5*quiverSubSample*dLon;


    xQ=X(1,floor(quiverSubSample/2):quiverSubSample:end)';
    yQ=Y(floor(quiverSubSample/2):quiverSubSample:end,1)';

    u1Q=u1(floor(quiverSubSample/2):quiverSubSample:end,floor(quiverSubSample/2):quiverSubSample:end);
    v1Q=v1(floor(quiverSubSample/2):quiverSubSample:end,floor(quiverSubSample/2):quiverSubSample:end);
    
    u2Q=u2(floor(quiverSubSample/2):quiverSubSample:end,floor(quiverSubSample/2):quiverSubSample:end);
    v2Q=v2(floor(quiverSubSample/2):quiverSubSample:end,floor(quiverSubSample/2):quiverSubSample:end);
    
    Xq=repmat(xQ,1,length(yQ));
    Yq=repmat(yQ,length(xQ),1);
    arrows=~isnan(u1Q) & ~isnan(v1Q) & ~isnan(u2Q) & ~isnan(v2Q); 

%     arrows=~isnan(u1Q) & ~isnan(v1Q) & ~isnan(u2Q) & ~isnan(v2Q) & ...
%         abs(u1Q)<4 & abs(v1Q)<4 & abs(u2Q)<4 & abs(v2Q)<4; 
    
    % Plot wind speeds
    [~,ch]=m_contourf(X,Y,windSpeed',axisMin:step:axisMax);
    set(ch,'edgecolor','none');
    clear ch;

    caxis([axisMin,axisMax]);
    colormap(cMap);
    c=colorbar;
    set(c,'YTick',axisMin:step:axisMax,'FontSize',12,'FontName','Times New Roman')
    ylabel(c,'m/s','FontSize',12,'FontName','Times New Roman');
    hold on;
    
    % Convert to units per deg...
    mag1=max(max(max(sqrt(u1Q(arrows).^2+v1Q(arrows).^2))));
    mag2=max(max(max(sqrt(u2Q(arrows).^2+v2Q(arrows).^2))));
    
    if mag1<mag2
        
        u1Q=arrowScale*u1Q./mag2;
        v1Q=arrowScale*v1Q./mag2;
        
        u2Q=arrowScale*u2Q./mag2;
        v2Q=arrowScale*v2Q./mag2;
        
        m_quiver(Xq(arrows),Yq(arrows),u1Q(arrows),v1Q(arrows),0,'color',[.8 0 0],'LineWidth',0.5);
        m_quiver(Xq(arrows),Yq(arrows),u2Q(arrows),v2Q(arrows),0,'color',[0 0 .8],'LineWidth',0.5);
    else
        
        u1Q=arrowScale*u1Q./mag1;
        v1Q=arrowScale*v1Q./mag1;
        
        u2Q=arrowScale*u2Q./mag1;
        v2Q=arrowScale*v2Q./mag1;
        
        m_quiver(Xq(arrows),Yq(arrows),u1Q(arrows),v1Q(arrows),0,'color',[.8 0 0],'LineWidth',0.5);
        m_quiver(Xq(arrows),Yq(arrows),u2Q(arrows),v2Q(arrows),0,'color',[0 0 .8],'LineWidth',0.5);
    end
        
    m_gshhs_i('linewidth',.25,'color','black');
    
    % Plot topography contour, note MC elevation is 4509m in New Guinea
    m_tbase('contour',1000:1000:2000,'color',tMap(1,:),'LineWidth',0.25);
    m_tbase('contour',2000:1000:3000,'color',tMap(2,:),'LineWidth',0.25);
    m_tbase('contour',3000:1000:4000,'color',tMap(3,:),'LineWidth',0.25);
    m_tbase('contour',4000:1000:5000,'color',tMap(4,:),'LineWidth',0.25);
    
    m_grid('xtick',0:tick:360,'ytick',-90:tick:90,'FontSize',12,'FontName','Times New Roman');
    hold off;

    clear ch Xq Yq u1Q v1Q u2Q v2Q arrows;

end