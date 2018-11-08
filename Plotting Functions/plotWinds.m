function h=plotWinds(X,Y,u,v,axisMin,step,axisMax,windSpeed,plotArgs)
    
    if size(X,1)==1 || size(X,2)==1
        [X, Y]=meshgrid(X,Y);
    end
    
    h=figure('units','centimeters','pos',[0 0 plotArgs.figureX plotArgs.figureY]);
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    if plotArgs.labels
        left = outerpos(1) + 1.5*ti(1);
        bottom = outerpos(2) + 1.5*ti(2);
        ax_width = outerpos(3) - 2.5*ti(1) - 2.5*ti(3);
        ax_height = outerpos(4) - 2*ti(2) - 2*ti(4);
        ax.Position = [left bottom ax_width ax_height];
    else
        left = outerpos(1) + .5*ti(1);
        bottom = outerpos(2) + .5*ti(2);
        ax_width = outerpos(3) - .5*ti(1) - .5*ti(3);
        ax_height = outerpos(4) - 1*ti(2) - 1*ti(4);
        ax.Position = [left bottom ax_width ax_height];
    end
    
    cMap=cbrewer('seq','Reds',length(axisMin:step:axisMax),'pchip');
    cMap=cMap(2:end,:);
    tMap=flipud(cbrewer('seq','Greys',5,'pchip'));
    
    dLon=X(1,2)-X(1,1);
    arrowScale=(7/8)*plotArgs.quiverSubSample*dLon;
    
    xQ=X(1,ceil(plotArgs.quiverSubSample/2):plotArgs.quiverSubSample:end)';
    yQ=Y(ceil(plotArgs.quiverSubSample/2):plotArgs.quiverSubSample:end,1)';

    uQ=u(ceil(plotArgs.quiverSubSample/2):plotArgs.quiverSubSample:end,ceil(plotArgs.quiverSubSample/2):plotArgs.quiverSubSample:end);
    vQ=v(ceil(plotArgs.quiverSubSample/2):plotArgs.quiverSubSample:end,ceil(plotArgs.quiverSubSample/2):plotArgs.quiverSubSample:end);

    Xq=repmat(xQ,1,length(yQ));
    Yq=repmat(yQ,length(xQ),1);
    arrows=~isnan(uQ) & ~isnan(vQ);
    
    mag=max(max(max(sqrt(uQ(arrows).^2+vQ(arrows).^2))));
    uQ=arrowScale*uQ./mag;
    vQ=arrowScale*vQ./mag;

    % Plot wind speeds
    [~,ch]=m_contourf(X,Y,windSpeed',axisMin:step:axisMax);
    set(ch,'edgecolor','none');
    clear ch;

    caxis([axisMin,axisMax]);
    colormap(cMap);
    if plotArgs.labels==true
        c=colorbar;
        set(c,'YTick',axisMin:step:axisMax,'FontSize',12,'FontName','Times New Roman');
        ylabel(c,'m/s','FontSize',12,'FontName','Times New Roman'); 
    else
        %c=colorbar;
    end

    hold on;
    
%     Plot political boundaries
%     M=m_shaperead('STE_2016_AUST'); 
% %     M=m_shaperead('VIC_LOCALITY_POLYGON_shp');
%     for k=1:size(M.ncst,1)
%         m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color',[0 0 0],'LineWidth',0.25); 
%     end 

    % Plot arrows
    m_quiver(Xq(arrows),Yq(arrows),uQ(arrows),vQ(arrows),0,'color','k','LineWidth',0.5);
    m_gshhs_i('LineWidth',0.25,'color','black');
    colormap(cMap);
    % Plot topography contour, note MC elevation is 4509m in New Guinea
    m_tbase('contour',1000:1000:2000,'color',tMap(1,:),'LineWidth',0.25);
    m_tbase('contour',2000:1000:3000,'color',tMap(2,:),'LineWidth',0.25);
    m_tbase('contour',3000:1000:4000,'color',tMap(3,:),'LineWidth',0.25);
    m_tbase('contour',4000:1000:5000,'color',tMap(4,:),'LineWidth',0.25);
    
    % Plot the grid
    if plotArgs.labels==true
        m_grid('xtick',0:plotArgs.tick:360,'ytick',-90:plotArgs.tick:90,'FontSize',12,'FontName','Times New Roman','glinewidth',.5);
    else
        m_grid('xtick',0:plotArgs.tick:360,'ytick',-90:plotArgs.tick:90,'xticklabels',[],'yticklabels',[],'glinewidth',.5);
    end
    
    clear ch Xq Yq uQ vQ arrows;

end