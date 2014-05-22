function animate(t,cs,ffvec,vvec,fig)

close all

cs = real(cs);
ffvec = real(ffvec);
vvec = real(vvec);

csz = size(cs);
steps = csz(2)-1;
tlen = max(size(t));
xvec = linspace(0,1,steps);

if(strcmp(fig,'s'))
    figure
    for i=1:tlen
        plot(xvec,cs(i,1:end-1))
        axis([0 1 0 1])
        M(i) = getframe(gcf);
    end

elseif(strcmp(fig,'d'))
    scrsz = get(0,'ScreenSize');
    figure('OuterPosition',[3/4 scrsz(4)/2 5*scrsz(3)/8 scrsz(4)/2])
    for i=1:tlen
        subplot(1,2,1)
        plot(xvec,cs(i,1:end-1))
        xlabel('Particle Length')
        ylabel('Concentration')
        axis([0 1 0 1])
        subplot(1,2,2)
        plot(ffvec,vvec)
        axis([0 1 3 3.6])
        hold on
        plot(ffvec(i),vvec(i),'ro')
        xlabel('Filling Fraction')
        ylabel('Voltage')
        M(i) = getframe(gcf);
        cla
    end

end









return;
