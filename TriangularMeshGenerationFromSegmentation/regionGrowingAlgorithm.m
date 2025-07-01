function [validTriangles, inverseTriangles] = regionGrowingAlgorithm(Tri,initialTriangle)

%% Inputs

% Tri: Connectivity matrix of the total mesh that corresponds to the set of
%      validTriangles and inverseTriangles.

% initialTriangle: Index referred to Tri of the initial triangle.

%% Outputs

% validTriangles: Connectivity matrix that represents the triangles that
%                 correspond to the set of valid triangles.

% inverseTriangles: Connectivity matrix that represents the triangles that
%                   correspond to the set of invalid triangles.

%% Code explanation

% This code will apply a region-growing algorithm to obtain the set of
% validTriangles and inverseTriangles. It will start from the
% initialTriangle to obtain the triangles that satisfy the conditions of
% being validTriangles and inverseTriangles.

%% Algorithm

% The algorith will be the following, identify the bad triangles around
% GoodTriangleIndex (2 shared vertices) and the good ones (1 sahred
% vertices). Inside GoodTriangles there are 3 that is bad (the one that
% share the top vertex) so we need to intersect them to elimite that.

% We first initialize the matrix of GoodTriangles and BadTriangles.

BadTriangles =[];

GoodTriangles = [Tri(initialTriangle,:)];

% The algorith will be the following, identify the bad triangles around
% GoodTriangleIndex (2 shared vertices) and the good ones (1 sahred
% vertices). Inside GoodTriangles there are 3 that is bad (the one that
% share the top vertex) so we need to intersect them to elimite that.

% We first initialize the matrix of GoodTriangles and BadTriangles.

for i = 1:NumRowsTri
    C0 = intersect(Tri(GoodTriangleIndex,:),Tri(i,:));
    if length(C0)==2
        BadTriangles = [BadTriangles;Tri(i,:)];
    end 
end

[NumRowsBad,NumColumnsBad] = size(BadTriangles);

Tri(ismember(Tri,BadTriangles,'rows'),:)=[];
Tri(ismember(Tri,GoodTriangles,'rows'),:)=[];

[NumRowsTri,NumColumnsTri] = size(Tri);

for c = 1:60
    disp(c)
    a = GoodTriangles;
    
    for i = 1:NumRowsTri
        
        intersectGood = arrayfun(@(x) intersect(Tri(i, :),x), GoodTriangles, 'UniformOutput', false);
        intersectBad= arrayfun(@(x) intersect(Tri(i, :),x), BadTriangles, 'UniformOutput', false);
        
        sumIntG= arrayfun(@(x) (~isempty(x{1})), intersectGood, 'UniformOutput', false);
        compareG = sum(cell2mat(sumIntG), 2);

        sumIntBad= arrayfun(@(x) (~isempty(x{1})), intersectBad, 'UniformOutput', false);
        compareBad = sum(cell2mat(sumIntBad), 2);

        if any(compareG == 1) && any(compareBad==2)
            a = [a;Tri(i,:)];
            Tri(i,:)=0;  
        end

    end
    GoodTriangles = a;

    Tri( ~any(Tri,2), : ) = [];
    Tri(ismember(Tri,GoodTriangles,'rows'),:)=[];
    [NumRowsTri,NumColumnsTri] = size(Tri);
    [NumRowsGood,NumColumnsGood] = size(GoodTriangles);    

    for l = 1:NumRowsTri
        intersectGood = arrayfun(@(x) intersect(Tri(l, :),x), GoodTriangles, 'UniformOutput', false);
        sumIntG= arrayfun(@(x) (~isempty(x{1})), intersectGood, 'UniformOutput', false);
        compareG = sum(cell2mat(sumIntG), 2);

        if any(compareG==2)
            BadTriangles = [BadTriangles;Tri(l,:)];
        end
        [NumRowsGood,NumColumsGood] = size(GoodTriangles);
    end

    [NumRowsBad,NumColumnsBad] = size(BadTriangles);
    Tri(ismember(Tri,BadTriangles,'rows'),:)=[];
    [NumRowsTri,NumColumnsTri] = size(Tri);

    if isempty(Tri)
        break
    end
    
end

inverseTriangles = BadTriangles;
validTriangles = GoodTriangles;
