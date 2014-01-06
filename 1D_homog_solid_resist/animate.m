function animate(t,cpcs,clv,plv,ffvec,vvec,disc,fig,output)

close all
tlen = max(size(t));
ss = disc.ss;
steps = disc.steps;
L = disc.L;

if fig=='e'
    % Plot the electrolyte concentration
    figure
    
    for i=1:tlen
        plot(clv,cpcs(i,1:ss+steps),'LineWidth',1)
        axis([clv(1) clv(end) 0 4])
        ylabel('Dimensionless Electrolyte Concentration')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
    end   
    
elseif fig=='s'
    figure
    % Construct solid concentration vector
    cst = cpcs(:,2*(ss+steps)+1:2*ss+3*steps);
    for i=1:tlen
        plot(plv,cst(i,:),'LineWidth',1)
        axis([plv(1) plv(end) 0 1])
        ylabel('Dimensionless Particle Concentration')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
    end  
elseif strcmp(fig,'mp')
    figure
    cst = cpcs(:,2*ss+3*steps+1:end-1);
    for i=1:tlen
        plot(plv,cst(i,:),'LineWidth',1)
        axis([plv(1) plv(end) 0 max(max(cst))])
        ylabel('Dimensionless Metal Potential')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
    end
elseif strcmp(fig,'po')
    figure
    ep = cpcs(:,ss+steps+1:2*(ss+steps));
    mp = cpcs(:,2*ss+3*steps+1:end-1);
    for i=1:tlen
        hold on
        plot(plv,mp(i,:),'LineWidth',1)
        plot(clv,ep(i,:),'LineWidth',1)
        hold off
        axis([clv(1) clv(end) 0 max(max([ep,mp]))])
        ylabel('Dimensionless Metal Potential')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
        cla
    end    
elseif fig=='d'
    scrsz = get(0,'ScreenSize');  %(1 1 width height)
    % Position -> Left Bottom Width Height
    figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2])
    axes('FontName','Arial','FontSize',12)
    % Plot the voltage and filling fraction on the same plot
    cst = cpcs(:,2*(ss+steps)+1:2*ss+3*steps);
    for i=1:tlen      
        subplot(1,2,1)
        plot(ffvec,vvec,'LineWidth',2)
        hold on
        plot(ffvec(i),vvec(i),'ro','MarkerSize',10','MarkerFaceColor','r')
        xlabel('Filling Fraction')
        ylabel('Voltage')
        axis([0 1 min(vvec) max(vvec)])
        hold off
        subplot(1,2,2)
        plot(plv,cst(i,:),'LineWidth',2)
        axis([plv(1) plv(end) 0 1])
        ylabel('Dimensionless Particle Concentration')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
    end
elseif strcmp(fig,'se')
    scrsz = get(0,'ScreenSize');  %(1 1 width height)
    % Position -> Left Bottom Width Height
    figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2])
    axes('FontName','Arial','FontSize',12)
    % Plot the voltage and filling fraction on the same plot
    cst = cpcs(:,2*(ss+steps)+1:2*ss+3*steps);
     for i=1:tlen      
        subplot(1,2,1)
        plot(clv,cpcs(i,1:ss+steps),'LineWidth',2)
        axis([clv(1) clv(end) 0 4])
        ylabel('Dimensionless Electrolyte Concentration')
        xlabel('Dimensionless Electrode Length')
        subplot(1,2,2)
        plot(plv,cst(i,:),'LineWidth',2)
        axis([plv(1) plv(end) 0 1])
        ylabel('Dimensionless Particle Concentration')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
     end
elseif strcmp(fig,'seo')
    figure
    % Construct solid concentration vector
    cst = cpcs(:,2*(ss+steps)+1:2*ss+3*steps);
    for i=1:tlen
        hold on
        plot(plv,cst(i,:),'b-','LineWidth',2)
        axis([clv(1) clv(end) 0 3])
        plot(clv,cpcs(i,1:ss+steps),'r-','LineWidth',2)
        ylabel('ND Concentration')
        xlabel('ND Electrode Length')
        %legend('Solid Conc.','Elect. Conc.','Location','Northeast')
        hold off
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
        cla
    end  
end


if output == 1
    
    %movie(M,1,fps);
    movie2avi(M,'C:\Users\ViN\Desktop\movie.avi','fps',15)
    %movie2avi(M,'D:\movie.avi','fps',15)
end


return;