function animate_dreyer(t,cs,disc,output)

close all
tlen = max(size(t));


figure

for i=1:tlen
    csmat = reshape(cs(i,1:end-1),disc.Ny,disc.Nx);
    surf(csmat);    
    axis([0 disc.Nx 0 disc.Ny 0 1])
    M(i) = getframe(gcf);
end

if(output)
    movie2avi(M,'C:\Users\trf\Desktop\movie.avi')    
end


return;