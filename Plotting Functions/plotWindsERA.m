function mag=plotWindsERA(X,Y,u,v,time,timeRange,level,contour3,minC3,stepC3,maxC3,contour1,minC1,stepC1,maxC1,contour2,minC2,stepC2,maxC2,quiverSubSample,tick,labels,units)
        
    close all;

    if nargin<20
        labels=false;
        units='m/s';
        quiverSubSample=4;
        tick=2;
    end
    
    savePlots=true;
    
    % Standard ranges
    % Temp
    % max 310 min 280 step 1
    % cape 
    % max 1400 min 0 step 100
    % tcc
    % 0 1 .1
    
    if length(size(u))>3
        u=u(:,:,level,:);
        v=v(:,:,level,:);
        contour1=squeeze(contour1(:,:,level,:));
        contour2=squeeze(contour2(:,:,level,:));
        contour3=squeeze(contour3(:,:,level,:));
    end
    
    m_proj('equidistant cylindrical','longitudes',[143 149], ...
    'latitudes',[-38 -34]);
    [X, Y]=meshgrid(X,Y);
    
    [~,tInd]=sort(time);
    time=time-time(1);

    windSpeed=sqrt(u.^2+v.^2);
    maxSpeed=ceil(max(max(max(windSpeed))));
    
    if maxSpeed<5
        step=0.5;
    elseif maxSpeed<10
        step=1;
    else
        step=2;
    end

    mag=12;

    for i=1:length(timeRange)
        
        titleString=sprintf('%u hours',time(tInd(timeRange(i))));
        % Since 11:00 AEDT 20/03/2013
%         plotWinds(X,Y,u(:,:,tInd(timeRange(i))),v(:,:,tInd(timeRange(i))),0,step,maxSpeed,windSpeed(:,:,tInd(timeRange(i))),quiverSubSample,tick,labels);
        contourPlot(X,Y,contour2(:,:,tInd(timeRange(i))),minC2,maxC2,stepC2,tick,units);
        hold on;
        
        title(titleString,'FontSize',12,'FontName','Times New Roman','FontWeight','normal');
        
        %Plot wind vectors
        dLon=X(1,2)-X(1,1);
        arrowScale=(4/8)*quiverSubSample*dLon;

        xQ=X(1,ceil(quiverSubSample/2):quiverSubSample:end)';
        yQ=Y(ceil(quiverSubSample/2):quiverSubSample:end,1)';

        uQ=u(ceil(quiverSubSample/2):quiverSubSample:end,ceil(quiverSubSample/2):quiverSubSample:end,tInd(timeRange(i)));
        vQ=v(ceil(quiverSubSample/2):quiverSubSample:end,ceil(quiverSubSample/2):quiverSubSample:end,tInd(timeRange(i)));

        Xq=repmat(xQ,1,length(yQ));
        Yq=repmat(yQ,length(xQ),1);
        arrows=~isnan(uQ) & ~isnan(vQ);

        uQ=arrowScale*uQ./mag;
        vQ=arrowScale*vQ./mag;
    
         % Plot political boundaries
        M=m_shaperead('STE_2016_AUST'); 
    %     M=m_shaperead('VIC_LOCALITY_POLYGON_shp');
        for k=1:size(M.ncst,1)
            m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color',[0 0 0],'LineWidth',0.25); 
        end 

        % Plot arrows
        m_quiver(Xq(arrows),Yq(arrows),uQ(arrows),vQ(arrows),0,'color','k','LineWidth',0.5);
        m_gshhs_i('LineWidth',0.25,'color','black');
        
        % Plot the grid
        if labels==true
            m_grid('xtick',0:tick:360,'ytick',-90:tick:90,'FontSize',12,'FontName','Times New Roman','glinewidth',.5);
        else
            m_grid('xtick',0:tick:360,'ytick',-90:tick:90,'xticklabels',[],'yticklabels',[],'glinewidth',.5);
        end
        
        [c, h]=m_contour(X,Y,contour1(:,:,tInd(timeRange(i)))',minC1:stepC1:maxC1,'color',[0 0 0]);
        clabel(c,h,'FontSize',12,'FontName','Times New Roman','LabelSpacing',2400);
        [c, h]=m_contour(X,Y,contour3(:,:,tInd(timeRange(i)))',minC3:stepC3:maxC3,'--k');
        clabel(c,h,'FontSize',12,'FontName','Times New Roman','LabelSpacing',800);
       
        m_line(146.0011,-35.9819,'marker','o','markersize',5,'color',[0 0 0],'MarkerFaceColor',[0 0 0]);
%         m_text(145.7,-36.2,'Mulwala','FontSize',12,'FontName','Times New Roman')
        
        if savePlots
            print(gcf,'-dpng',[num2str(i) '.png'],'-r800');
        end
    end
    
hold off;

end