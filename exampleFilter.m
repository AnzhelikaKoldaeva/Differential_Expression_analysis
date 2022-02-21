
function exampleFilter(str_matr)
% function exampleFilter
%
% Example of how the Kalman filter performs for real data

%load raw

raw = double(str_matr);

%Apply filter
k=Kalman_Stack_Filter(raw);
k75=Kalman_Stack_Filter(raw,0.75);

figure();

%Set up image window
minMax1=[min(raw(:)),max(raw(:))];
minMax2=[min(k(:)),max(k(:))/2];
minMax3=[min(k75(:)),max(k75(:))/4];
clf, colormap  default  %gray

subplot(1,3,1)
imagesc(raw(:,:,1))
title('original')
set(gca,'clim',minMax1), axis off equal

subplot(1,3,2)
imagesc(k(:,:,1))
title('filtered, gain=0.5')
set(gca,'clim',minMax2), axis off equal

subplot(1,3,3)
imagesc(k75(:,:,1))
title('filtered, gain=0.75')
set(gca,'clim',minMax3), axis off equal


%Loop movie 
disp('crtl-c to stop movie')
while 1
    for i=1:size(k,3)
        
        subplot(1,3,1)
        set(get(gca,'children'),'CData',raw(:,:,i))
        
        subplot(1,3,2)
        set(get(gca,'children'),'CData',k(:,:,i))
        
        subplot(1,3,3)
        set(get(gca,'children'),'CData',k75(:,:,i))
        
        
        pause(0.05)
        drawnow
    end
end
