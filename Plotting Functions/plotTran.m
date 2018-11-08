function plotTran(varargin)
   
    % Input the coords for the map projection showing projection lines, 
    % then input a variable number of tran structures.
    
    close all;
    
    % Personal Macbook.
    if ismac
        folderLoc='/Volumes/Ewan''s Hard Drive/Figures/';
    end

    % Uni Unix box machines.
    if isunix && not(ismac)
        username=char(java.lang.System.getProperty('user.name'));
        folderLoc=['/media/' username '/Ewan''s Hard Drive/Figures/'];
        clear username;
    end

    startLat=-13;
    endLat=5;
    startLon=94;
    endLon=152;

    plotVarBars=false;
    plotPbar=false;
    plotMap=true;
    plotLegend=0;
    savePlots=1;
    tick=10;
    sigLevel=0.05;
    
    if length(varargin)>8
        cMap=cbrewer('qual','Dark2',length(varargin),'pchip');
    else
        cMap=cbrewer('qual','Dark2',8,'pchip');
    end
    
    if plotMap
        
        figure('units','centimeters','pos',[0 0 14 14]);
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + 1.5*ti(1);
        bottom = outerpos(2) + 1.5*ti(2);
        ax_width = outerpos(3) - 2*ti(1) - 2*ti(3);
        ax_height = outerpos(4) - 2*ti(2) - 2*ti(4);
        ax.Position = [left bottom ax_width ax_height];
        
        m_proj('equidistant cylindrical','longitudes',[94 152], ...
        'latitudes',[-13 6]);

        tMap=flipud(cbrewer('seq','Greys',5,'pchip'));
        
        hold on;

        % Plot lines

        m_gshhs_i('LineWidth',0.25,'color','black');
        
        for i=1:length(varargin)
            
            colourVar='k';
%             colourVar=cMap(i,:);
            
            % First plot a line indicating "coast"
            m_plot([varargin{i}.startLon varargin{i}.endLonCoast], ...
                    [varargin{i}.startLat varargin{i}.endLatCoast],...
                    'Color',colourVar,'LineWidth',.25);
            
            dxC=(varargin{i}.endLonCoast-varargin{i}.startLon)...
                /(varargin{i}.nTrans-1);
            dyC=(varargin{i}.endLatCoast-varargin{i}.startLat)...
                /(varargin{i}.nTrans-1);

            skip=ceil(varargin{i}.nTrans/25);
            
            for j=1:skip:varargin{i}.nTrans
                m_plot([varargin{i}.startLon+(j-1)*dxC varargin{i}.endLonTran+(j-1)*dxC], ...
                    [varargin{i}.startLat+(j-1)*dyC varargin{i}.endLatTran++(j-1)*dyC],...
                    'Color',colourVar,'LineWidth',.25);
            end
            
        end
        
        % Plot topography contour, note MC elevation is 4509m in New Guinea
%         m_tbase('contour',1000:1000:2000,'color',tMap(1,:),'LineWidth',0.25);
%         m_tbase('contour',2000:1000:3000,'color',tMap(2,:),'LineWidth',0.25);
%         m_tbase('contour',3000:1000:4000,'color',tMap(3,:),'LineWidth',0.25);
%         m_tbase('contour',4000:1000:5000,'color',tMap(4,:),'LineWidth',0.25);

        % Plot the grid
        m_grid('xtick',0:tick:360,'ytick',-90:tick:90,'FontSize',12,'FontName','Times New Roman','glinewidth',5);
        if savePlots
            print(gcf,'-dsvg',[folderLoc,'/map ',num2str(varargin{1}.time)]);
        end
        
        hold off;

        clear ch Xq Yq uQ vQ arrows;
        
    end
    
    % Set figure window
    figure('units','centimeters','pos',[0 0 7 5.5]);
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + 2*ti(1);
    bottom = outerpos(2) + 2*ti(2);
    ax_width = outerpos(3) - 2.25*ti(1) - 2.25*ti(3);
    ax_height = outerpos(4) - 2.375*ti(2) - 2.375*ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    hold on 
    
    legendInfo=cell(1,length(varargin));
    
    % Plot dotted line to indicate 0
    plot([0 1000],[0 0],'--','Color',[0 0 0],'LineWidth',1);
    
    for i=1:length(varargin)
        h(i)=plot(varargin{i}.distance,varargin{i}.pertProj,'Color',cMap(i,:),'LineWidth',1);
        
        % Find first statistically insignificant p-Value.
        k=find(varargin{i}.pProj>sigLevel,1,'first');
        % Plot vertical bar to show where data becomes insignificant.
        if ~isempty(k) && plotPbar
            plot([varargin{i}.distance(k) varargin{i}.distance(k)],[-2 2],'--','Color',[.4 .4 .4],'LineWidth',1);
        end
            
        if plotVarBars
            plot(varargin{i}.distance,varargin{i}.pertProj+sqrt(abs(varargin{i}.varProj)),'Color',cMap(i,:),'LineWidth',1,'LineStyle','--');
            plot(varargin{i}.distance,varargin{i}.pertProj-sqrt(abs(varargin{i}.varProj)),'Color',cMap(i,:),'LineWidth',1,'LineStyle','--');
        end
        
        legendInfo{i}=varargin{i}.label;
        axis([0 1000 -2 2]);
        yticks(-10:.5:10);
        xticks(0:200:1000);
    end
    
    set(gca,'FontSize',12,'FontName','Times New Roman');
    xlabel('km','FontSize',12,'FontName','Times New Roman');
    ylabel('m/s','FontSize',12,'FontName','Times New Roman');
    
    if plotLegend
        legend(h(1:end),legendInfo,'FontSize',10,'FontName','Times New Roman');
    end
    
    if savePlots
        print(gcf,'-dsvg',[folderLoc,'/plot ',varargin{1}.label,'.svg']);
    end
        
    hold off;

end