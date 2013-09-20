function animate(t,cpcs,clv,plva,plvc,ffvec,vvec,disc,fig,output)

dimlen = 1.23;

close all
tlen = max(size(t));
ss = disc.ss;
asteps = disc.asteps;
csteps = disc.csteps;
len = ss+asteps+csteps;
anlen = disc.anlen;
catlen = disc.catlen;
seplen = disc.seplen;

abeta0 = disc.abeta;
cbeta0 = disc.cbeta;
aporos = disc.aporos;
cporos = disc.cporos;

dx = ones(1,len);
dx(1:asteps) = anlen/asteps;
dx(asteps+1:asteps+ss) = seplen/ss;
dx(asteps+ss+1:len) =  catlen/csteps;

if fig=='e'
    % Plot the electrolyte concentration
    figure
    for i=1:tlen
        hold on
        plot(clv,cpcs(i,1:len),'LineWidth',1)
        vline(0,'-r');
        vline(seplen,'r');
        hold off
        axis([clv(1) clv(end) 0 4])
        ylabel('Dimensionless Electrolyte Concentration')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
        % ELECTROLYTE INTEGRAL
        ctot = dx.*cpcs(i,1:len);
        ctot(1:asteps) = ctot(1:asteps)*aporos;
        ctot(asteps+ss+1:end) = ctot(asteps+ss+1:end)*cporos;
%         disp(sum(ctot))
        cla
    end 
elseif strcmp(fig,'dseco')
    scrsz = get(0,'ScreenSize');  %(1 1 width height)
    % Position -> Left Bottom Width Height
    figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2])
    % Construct solid concentration vector
    acst = cpcs(:,2*len+1:2*len+asteps);
    ccst = (cpcs(:,2*len+asteps+1:2*len+asteps+csteps)+cpcs(:,2*len+asteps+csteps+1:2*len+asteps+2*csteps))/2;
    clv = clv*dimlen;           % Scale the dimensionless length by the real length
    % For each time step, we need to find the gold/red/blue areas and fill
    for i=1:tlen/4
        % First the electrolyte concentration
        subplot(1,2,1)
        hold on
        plot(clv,cpcs(i,1:len),'-b','LineWidth',2)
        plot(0,0)
        vline(0,'r');
        vline(seplen*dimlen,'r');
        hold off
        axis([clv(1) clv(end) 0 4])
        ylabel('Dimensionless Electrolyte Concentration','FontSize',14)
        xlabel('Electrode Length (mm)','FontSize',14)
        set(gca,'FontSize',14)
        % Now the solid concentration
        % Need to plot the fill color - requires cycling through concentrations
        subplot(1,2,2)
        len2 = size(ccst);
        ymax = 1;
        rmax = 2;
        for j=1:len2(2)
            if ccst(i,j) > .597
                ymax = j;
            end
            if ccst(i,j) > .08
                rmax = j;
            end            
        end
        hold on
        plot(plvc*dimlen,ccst(i,:),'-k','LineWidth',.5)
        plot(plva*dimlen,acst(i,:),'-b','LineWidth',2)
        area(plvc(1:ymax)*dimlen,ccst(i,1:ymax),'FaceColor',[.668 .625 0],'LineStyle','none')
        area(plvc(ymax:rmax)*dimlen,ccst(i,ymax:rmax),'FaceColor',[.617 .043 .059],'LineStyle','none')
        area(plvc(rmax:end)*dimlen,ccst(i,rmax:end),'FaceColor',[0 .203 .441],'LineStyle','none')
        plot(clv,cpcs(i,1:len),'--k','LineWidth',2)
        vline(0,'r');
        vline(seplen*dimlen,'r');
        hold off
%         axis([clv(1) clv(end) 0 1])
        axis([dimlen 3 0 1])
        ylabel('Dimensionless Particle Concentration','FontSize',14)
        xlabel('Electrode Length (mm)','FontSize',14)
        set(gca,'FontSize',14)
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
        cla
        subplot(1,2,1)
        cla
        
    end  
elseif strcmp(fig,'dse')
    scrsz = get(0,'ScreenSize');  %(1 1 width height)
    % Position -> Left Bottom Width Height
    figure('Position',[1 scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2])
    % Construct solid concentration vector
    acst = cpcs(:,2*len+1:2*len+asteps);
    ccst = (cpcs(:,2*len+asteps+1:2*len+asteps+csteps)+cpcs(:,2*len+asteps+csteps+1:2*len+asteps+2*csteps))/2;
    clv = clv*dimlen;           % Scale the dimensionless length by the real length
    % For each time step, we need to find the gold/red/blue areas and fill
    for i=1:tlen/2
        % First the electrolyte concentration
        subplot(1,2,1)
        hold on
        plot(clv,cpcs(i,1:len),'-b','LineWidth',2)
        plot(0,0)
        vline(0,'r');
        vline(seplen,'r');
        hold off
        axis([clv(1) clv(end) 0 4])
        ylabel('Dimensionless Electrolyte Concentration','FontSize',14)
        xlabel('Electrode Length (mm)','FontSize',14)
        set(gca,'FontSize',14)
        % Now the solid concentration
        % Need to plot the fill color - requires cycling through concentrations
        subplot(1,2,2)
        len2 = size(ccst);
        ymax = 1;
        rmax = 2;
        for j=1:len2(2)
            if ccst(i,j) > .597
                ymax = j;
            end
            if ccst(i,j) > .08
                rmax = j;
            end            
        end
        hold on
        plot(plvc,ccst(i,:),'-k','LineWidth',.5)
        plot(plva,acst(i,:),'-b','LineWidth',2)
        area(plvc(1:ymax),ccst(i,1:ymax),'FaceColor',[.668 .625 0],'LineStyle','none')
        area(plvc(ymax:rmax),ccst(i,ymax:rmax),'FaceColor',[.617 .043 .059],'LineStyle','none')
        area(plvc(rmax:end),ccst(i,rmax:end),'FaceColor',[0 .203 .441],'LineStyle','none')
        vline(0,'r');
        vline(seplen,'r');
        hold off
%         axis([clv(1) clv(end) 0 1])
        axis([1 3 0 1])
        ylabel('Dimensionless Particle Concentration','FontSize',14)
        xlabel('Electrode Length (mm)','FontSize',14)
        set(gca,'FontSize',14)
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
        cla
        subplot(1,2,1)
        cla
        
    end     
elseif fig=='p'
    % Plot the electrolyte concentration
    figure
    
    for i=1:tlen
        hold on
        plot(clv,cpcs(i,len+1:2*len),'LineWidth',1)
        vline(0,'r');
        vline(seplen,'r');
        hold off
        axis([clv(1) clv(end) 0 4])
        ylabel('Dimensionless Electrolyte Concentration')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
        cla
    end   
    
elseif fig=='s'
    figure
    % Construct solid concentration vector
    acst = cpcs(:,2*len+1:2*len+asteps);
    ccst = (cpcs(:,2*len+asteps+1:2*len+asteps+csteps)+cpcs(:,2*len+asteps+csteps+1:2*len+asteps+2*csteps))/2;
    for i=1:tlen
        hold on
        plot(plvc,ccst(i,:),'LineWidth',1)
        plot(plva,acst(i,:),'LineWidth',1)
        vline(0,'r');
        vline(seplen,'r');
        hold off
        axis([clv(1) clv(end) 0 1])
        ylabel('Dimensionless Particle Concentration')
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
    acst = cpcs(:,2*len+1:2*len+asteps);
    ccst = (cpcs(:,2*len+asteps+1:2*len+asteps+csteps)+cpcs(:,2*len+asteps+csteps+1:2*len+asteps+2*csteps))/2;

    for i=1:tlen      
        subplot(1,2,1)
        plot(ffvec,vvec,'LineWidth',2)
        hold on
        plot(ffvec(i),vvec(i),'ro','MarkerSize',10','MarkerFaceColor','r')
        xlabel('Filling Fraction')
        ylabel('Voltage')
        axis([0 1 0 .5])
        hold off
        subplot(1,2,2)
        hold on
        plot(plvc,ccst(i,:),'LineWidth',1)
        plot(plva,acst(i,:),'LineWidth',1)
        vline(0,'r');
        vline(seplen,'r');
        hold off
        axis([clv(1) clv(end) 0 1])
        ylabel('Dimensionless Particle Concentration')
        xlabel('Dimensionless Electrode Length')
        set(gcf,'Renderer','zbuffer')       % Fix for Windows 7
        M(i) = getframe(gcf);
        cla
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
%     movie2avi(M,'C:\Users\trf\Desktop\movie.avi','fps',30)
    movie2avi(M,'D:\movie.avi','fps',20)
end


return;