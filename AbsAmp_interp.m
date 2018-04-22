mask = 20*log10(ImageBS)>50;
AbsImAmpVM2 = AbsImAmpVM1.*mask;
pixel = 1:512;
interpRe = zeros(512,512);
for i = 1:512
    zeroEle = find(AbsImAmpVM2(:,i)==0);
    noneZeroEle = find(AbsImAmpVM2(:,i));
    interpRe(:,i) = interp1(noneZeroEle,AbsImAmpVM2(noneZeroEle,i),1:512,'linear',0);
end
figure;imagesc(AbsImAmpVM2,[0 70]);colormap(jet);colorbar
figure;imagesc(interpRe,[0 70]);colormap(jet);colorbar
strain = zeros(510,512);
for j=1:512
    for i=1:512-m+1
        strain(i,j)=(6*(2*(i-i)-m+1))*interpRe(i,j)/(m^3-m)*i+(6*(2*(i+1-i)-m+1))...
            *interpRe(i+1,j)/(m^3-m)*(i+1)+(6*(2*(i+2-i)-m+1))*interpRe(i+2,j)/(m^3-m)*(i+2);...
            %+(6*(2*(i+3-i)-m+1))*AbsImAmpVM1(i+3,j)/(m^3-m)*(i+3)+(6*(2*(i+4-i)-m+1))*AbsImAmpVM1(i+4,j)/(m^3-m)*(i+4);
    end
end
Youngs=1./abs(strain);

Y3=193*50*medfilt2(Youngs,[2 2]);

% mask = (20*log10(abs(squeeze(ComplexFrames(:,:,100))))>40).*(Snr>50);
% gra = gradient(AbsImAmpVM1,5);
% figure;imagesc(medfilt2(abs(1./gra)*80,[3 3]),[0 300]);colormap(jet);title('gradient')
figure;imagesc(Y3,[0 300]);colormap(jet);title('after strain stry CompAfter');colorbar
% figure;imagesc(medfilt2(AbsImAmpVM1,[1 1]),[0 3000]);colormap(jet);title('Vibration Amp CompAfter')
% Elas=interp1(K(1:510),Y3,KES(1:510),'linear','extrap'); %Dispersion Compensation
% a = mean(abs(ComplexFrames(1:510,:,:)),3);
% b = medfilt2(1e4*Youngs./a,[2 2]);
% figure;imagesc(1e4*b,[0 80]);colormap(jet);title('after comp');
toc