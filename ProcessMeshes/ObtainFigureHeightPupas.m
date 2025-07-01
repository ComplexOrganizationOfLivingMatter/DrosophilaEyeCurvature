function ObtainFigureHeightPupas(pupasPath,outputPath)

% This code will obtain the height of the pupas in pupasPath.

pupasDir = dir(strcat(pupasPath,'*','mat'));

% Initialize arrays to hold file names and iteration numbers
fileNames = {pupasDir.name};
iterationNumbers = zeros(size(fileNames));

% Extract iteration numbers using a regular expression
for i = 1:length(fileNames)
    % Extract the number after 'Pressure_' (or any other unique identifier)
    match = regexp(fileNames{i}, 'Pressure_(\d+)', 'tokens');
    if ~isempty(match)
        iterationNumbers(i) = str2double(match{1}{1});
    else
        iterationNumbers(i) = NaN; % Handle files that don't match the pattern
    end
end

% Sort the files based on iteration numbers:

[~, sortedIndices] = sort(iterationNumbers);
sortedFileList = pupasDir(sortedIndices);

pupaHeightIterations = [];

for i = 1:size(sortedFileList,1)

    pupaName = sortedFileList(i).name;
    load(strcat(pupasPath,pupaName),'coordinates')
    pupaHeight = max(coordinates(:,3))-min(coordinates(:,3));
    pupaHeightIterations = [pupaHeightIterations,pupaHeight];

end

% Finally we save the result and the figure

name2savefig = strsplit(pupasPath,'/');
figure
plot(pupaHeightIterations)
xlabel('Iterations')
ylabel('Height (\mum)')
figTitle = strcat(name2savefig{end-2},' Height along RAW iterations');
title(figTitle)
savefig(strcat(outputPath,name2savefig{end-2},'_Height_Iterations_RAW.fig'));
save(strcat(outputPath,name2savefig{end-2},'_Height_Iterations_RAW.mat'),'pupaHeightIterations','-v7.3')

