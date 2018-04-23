mask = ImageBS>120;
AbsImAmpVM2 = AbsImAmpVM1.*mask.*[zeros(20,512);ones(430,512);zeros(512-450,512)];
pixel = 1:512;
interpRe = zeros(512,512);
for i = 1:512
    noneZeroEle = find(AbsImAmpVM2(:,i));
    if length(noneZeroEle)<=2
        continue
    else
    interpRe(:,i) = interp1(noneZeroEle,AbsImAmpVM2(noneZeroEle,i),1:512,'linear',0);
    end
end
figure;imagesc(AbsImAmpVM2,[0 70]);colormap(jet);colorbar
figure;imagesc(interpRe,[0 70]);colormap(jet);colorbar
strain = zeros(510,512);
m=4;
for j=1:512
    for i=1:512-m+1
        strain(i,j)=(6*(2*(i-i)-m+1))*interpRe(i,j)/(m^3-m)+(6*(2*(i+1-i)-m+1))...
            *interpRe(i+1,j)/(m^3-m)+(6*(2*(i+2-i)-m+1))*interpRe(i+2,j)/(m^3-m)...
            +(6*(2*(i+3-i)-m+1))*interpRe(i+3,j)/(m^3-m);%+(6*(2*(i+4-i)-m+1))*AbsImAmpVM1(i+4,j)/(m^3-m)*(i+4);
    end
end
Youngs=1./abs(strain+1e-7);
Youngs = Youngs.*mask(1:510,:).*[zeros(20,512);ones(430,512);zeros(510-450,512)];


Y3=Youngs;



figure;imagesc(medfilt2(Y3.*[zeros(50,512);ones(300,512);zeros(510-350,512)],[5 5]),[0 1500]);colormap(jet);title('after strain stry CompAfter');colorbar

