function [n,minDifference,minDifferenceIterations,pupaMeanDifferenceIterations] = ObtainMinimumDifference(pupasPath,adultPath,ouputPath)

% This code will compare the Gaussian Distribution info of the different
% iterations of the pupa and will compare them with the adult eyes in
% adultPath. Finally, we will save the coordinates with the minimum
% difference in outputPath.

pupasDir = dir(strcat(pupasPath,'*','mat'));
adultsDir = dir(strcat(adultPath,'*','mat'));

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

minDifferenceIterations = [];
pupaMeanDifferenceIterations = [];

for i = 1:size(sortedFileList,1)

    pupaName = sortedFileList(i).name;
    load(strcat(pupasPath,pupaName),'GaussianCurvaturePercentiles')
    GaussianCurvaturePercentiles1 = GaussianCurvaturePercentiles;
    pupaMeanDifference = [];

    for j = 1:size(adultsDir,1)

        adultName = adultsDir(j).name;
        load(strcat(adultPath,adultName),'GaussianCurvaturePercentiles')
        GaussianCurvaturePercentiles2 = GaussianCurvaturePercentiles;
        pupaMeanDifference = [pupaMeanDifference,1 - abs(corr(GaussianCurvaturePercentiles1(2:11)',GaussianCurvaturePercentiles2(2:11)',"type","Pearson"))];

    end

    pupaMeanDifferenceIterations = [pupaMeanDifferenceIterations;pupaMeanDifference];
    minDifferenceIterations = [minDifferenceIterations,mean(pupaMeanDifference)];

end

% So, the min difference and the iteration is:

minDifference = min(minDifferenceIterations);
n = find(minDifferenceIterations==minDifference);

% Finally, we save the results:

copyfile(strcat(pupasPath,sortedFileList(n).name),strcat(ouputPath,sortedFileList(n).name))