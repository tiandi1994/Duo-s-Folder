% mask = 20*log10(ImageBS)>48;
mask = Snr>50;
AbsImAmpVM2 = AbsImAmpVM1.*mask;
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

% interpRe = medfilt2(interpRe,[6 1]);
aa = zeros(512,512);
for i =1:512
   aa(:,i) =  -cwt((interpRe(:,i)),5,'gaus1');
end
% aa = movmean(aa,10,1);
figure;imagesc(aa,[0 300]);colormap(jet)