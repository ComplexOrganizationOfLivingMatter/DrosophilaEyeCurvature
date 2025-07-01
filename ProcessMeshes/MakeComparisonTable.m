function MakeComparisonTable(inputPath,outputPath)

% This code will use the .mat in inputPath to obtain the comparison tables
% wich will be saved in outputPath:

filesDir = dir(strcat(inputPath,'*','.mat'));
similarityTable = zeros(size(filesDir,1));

for i = 1:size(filesDir,1)

    fileName1 = filesDir(i).name;
    load(strcat(inputPath,fileName1),'GaussianCurvaturePercentiles');
    GaussianCurvaturePercentiles1 = GaussianCurvaturePercentiles;

    for j = (i+1):size(filesDir,1)

        fileName2 = filesDir(j).name;
        load(strcat(inputPath,fileName2),'GaussianCurvaturePercentiles');
        GaussianCurvaturePercentiles2 = GaussianCurvaturePercentiles;
        similarityTable(j,i) = 1 - abs(corr(GaussianCurvaturePercentiles1(2:end)',GaussianCurvaturePercentiles2(2:end)',"type","Pearson"));

    end

end

% Now that we have the similarity table, it's symmetric so we transpose and
% we sum the matrices:

similarityTableTransposed = transpose(similarityTable);
similarityTableSymmetric = similarityTable + similarityTableTransposed;

% Finally we write the table:

similarityTableSymmetric = array2table(similarityTableSymmetric);
similarityTable = array2table(similarityTable);
tableColumnNames = cell(1,size(filesDir,1));

for i = 1:size(filesDir,1)

    eyeNameForTable = filesDir(i).name;
    eyeNameForTable = erase(eyeNameForTable,'Eye_42_44hr_');
    eyeNameForTable = erase(eyeNameForTable,'Eye_40_42hr_');
    eyeNameForTable = erase(eyeNameForTable,'.mat');
    tableColumnNames{i} = eyeNameForTable;

end

% Now we add the columns and rows name and write the table:

similarityTable.Properties.VariableNames = tableColumnNames;
similarityTable.Properties.RowNames = tableColumnNames;
similarityTableSymmetric.Properties.VariableNames = tableColumnNames;
similarityTableSymmetric.Properties.RowNames = tableColumnNames;
tableNameFull = strcat('Full_Gaussian_Similarity_Table_Pupas.xls');
tableNameSymmetric= strcat('Symmetric_Gaussian_Similarity_Table_Pupas.xls');
writetable(similarityTable,strcat(outputPath,tableNameFull),'WriteRowNames',true);
writetable(similarityTableSymmetric,strcat(outputPath,tableNameSymmetric),'WriteRowNames',true);
