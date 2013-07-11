function animate_csmat(t,cpcs,csmat,disc,part,ffvec,vvec,fig,output)

close all
tlen = max(size(t));


figure
dims = size(csmat);

if strcmp(fig,'s')
    for i=1:tlen
        csmattmp = squeeze(csmat(i,:,:));
        surf(csmattmp);    
        axis([0 dims(3) 0 dims(2) 0 1])
        M(i) = getframe(gcf);
    end
elseif strcmp(fig,'c')
    for i=1:tlen
        plot(cpcs(i,1:disc.ss+disc.steps));
        axis([0 disc.ss+disc.steps 0 2])
        M(i) = getframe(gcf);
    end
end

if(output)
    movie2avi(M,'C:\Users\trf\Desktop\movie.avi')    
end


return;
