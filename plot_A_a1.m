% figure;
% while(1)
% for i=1:512
% % imagesc((A_a1(:,:,i)),[-1 1]);colormap(jet)
% plot(squeeze(Ap1(i,250,:)));
% ylim([0 80])
% title(num2str(i))
% drawnow
% % pause(0.5)
% end
% end

% figure
% for i = 3:512
%     imshow(angle(ComplexFrames(:,:,i))-angle((ComplexFrames(:,:,2))));title(i)
%     drawnow
% end



aa = squeeze(Phframes(156,250,:));
aaa = aa-(aa>1.5)*pi+(aa<-1.5)*pi;
figure;plot(aaa)
figure;plot(abs(fft(aaa)))