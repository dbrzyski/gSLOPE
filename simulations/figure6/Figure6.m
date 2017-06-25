
mainpath = 'D:\projects\gSLOPE'; %needs to be updated
cd(mainpath)
cd('.\figures\figure6\plots')

n1 = 200;
n2 = 200;
n3 = 200;

[y,x,z] = ndgrid(linspace(-3,3,n1),linspace(-3,3,n2),linspace(-3,3,n3));

%--------------------------------------------------
V = zeros(n1, n2, n3);
lambda = [2;0];

for ii = 1: n1
    for jj = 1:n2
        for kk =1:n3
            V(ii, jj, kk) = lambda'*sort(sqrt(2)*[norm(  [y(ii, jj, kk); x(ii, jj, kk)]  ); abs(z(ii, jj, kk))], 'descend');
        end
    end
end

cla
p = patch(isosurface(x,y,z,V,1));
p.FaceColor = 'cyan';
p.EdgeColor = 'none';
view(3);
camlight
axis equal
lighting gouraud
box on
%---------------------------------------------------
print('Figure6a.eps','-depsc','-tiff')
%---------------------------------------------------
%---------------------------------------------------

%--------------------------------------------------
V = zeros(n1, n2, n3);
lambda = [2;1];

for ii = 1: n1
    for jj = 1:n2
        for kk =1:n3
            V(ii, jj, kk) = lambda'*sort(sqrt(2)*[norm(  [y(ii, jj, kk); x(ii, jj, kk)]  ); abs(z(ii, jj, kk))], 'descend');
        end
    end
end

cla
p = patch(isosurface(x,y,z,V,1));
p.FaceColor = 'cyan';
p.EdgeColor = 'none';
view(3);
camlight
axis equal
lighting gouraud
box on
%---------------------------------------------------
print('Figure6b.eps','-depsc','-tiff')
%---------------------------------------------------
%---------------------------------------------------
%--------------------------------------------------
V = zeros(n1, n2, n3);
lambda = [2;2];

for ii = 1: n1
    for jj = 1:n2
        for kk =1:n3
            V(ii, jj, kk) = lambda'*sort(sqrt(2)*[norm(  [y(ii, jj, kk); x(ii, jj, kk)]  ); abs(z(ii, jj, kk))], 'descend');
        end
    end
end

cla
p = patch(isosurface(x,y,z,V,1));
p.FaceColor = 'cyan';
p.EdgeColor = 'none';
view(3);
camlight
axis equal
lighting gouraud
box on
%---------------------------------------------------
print('Figure6c.eps','-depsc','-tiff')
%---------------------------------------------------
%---------------------------------------------------




