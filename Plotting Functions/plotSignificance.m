function pStar=plotSignificance(binAv,pValue,alpha)
      
    h=figure('units','centimeters','pos',[0 0 20 6]);
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 

    left = outerpos(1) + .5*ti(1);
    bottom = outerpos(2) + .5*ti(2);
    ax_width = outerpos(3) - .5*ti(1) - .5*ti(3);
    ax_height = outerpos(4) - 1*ti(2) - 1*ti(4);
    ax.Position = [left bottom ax_width ax_height];

    m_proj('equidistant cylindrical','longitudes',[binAv.x(1) binAv.x(end)], ...
    'latitudes',[binAv.y(1) binAv.y(end)]);

    tMap=flipud(cbrewer('seq','Greys',5,'pchip'));
    
    [~, pStar]=FDRsignificance(pValue,alpha);
    [iInd, jInd]=find(pValue<=pStar);
    m_plot(binAv.x(iInd),binAv.y(jInd),'.','markers',3,'Color',[.1,.1,0.9]);

    m_gshhs_i('LineWidth',0.25,'color','black');
    % Plot topography contour, note MC elevation is 4509m in New Guinea
    m_tbase('contour',1000:1000:2000,'color',tMap(1,:),'LineWidth',0.25);
    m_tbase('contour',2000:1000:3000,'color',tMap(2,:),'LineWidth',0.25);
    m_tbase('contour',3000:1000:4000,'color',tMap(3,:),'LineWidth',0.25);
    m_tbase('contour',4000:1000:5000,'color',tMap(4,:),'LineWidth',0.25);
    

    m_grid('xtick',0:10:360,'ytick',-90:10:90,'xticklabels',[],'yticklabels',[],'glinewidth',.5);

    
    clear ch Xq Yq uQ vQ arrows;

end