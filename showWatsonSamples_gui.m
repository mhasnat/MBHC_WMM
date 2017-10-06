function showWatsonSamples_gui(hax, wmmSample, labels, numClass, mu)

%% Visualize data
color_ = {'b'; 'c'; 'm'; 'y'; 'k'; 'r'; 'g';};

colInd = 1;
for h=1:2:numClass*2
    scatter3(hax.axes1, wmmSample(find(labels==h),1), wmmSample(find(labels==h),2), wmmSample(find(labels==h),3), 'fill', color_{colInd}); hold on;
    scatter3(hax.axes1, wmmSample(find(labels==h+1),1), wmmSample(find(labels==h+1),2), wmmSample(find(labels==h+1),3), 'fill', color_{colInd}); hold on;
    
    colInd = colInd + 1;
    
    if(colInd>7)
        colInd = mod(colInd,7)
    end
end

hold on; [x,y,z] = sphere; surface(x,y,z,'FaceColor', 'none'); hold off;
% title(strcat('Samples from Watson mixture model')); xlabel('x'); ylabel('y'); zlabel('z');

