function obj = GetListInputUncertainVariables(obj)

%This function creates a list of the positions in the input structure that contain
%uncertain variables together with the uncertain variables

% InputPosList(nn).UV Gives the value of the UV in the n-th variable
% InputPosList(nn).Pos Gives the position of the n-th UV in the i-th index of the vector
% InputPosList(nn).StructPos Gives the structure position of the n-th UV in the i-th index of the vector, if there is no structure this is empty

UVList = [];
nn = 1;
for ii = 1:length(obj.Input)
    if ~isstruct(obj.Input{ii}{2})
        if strcmp(class(obj.Input{ii}{2}),'UncertainVariable')
            UVList(nn).UV = obj.Input{ii}{2};
            UVList(nn).Pos = ii;
            UVList(nn).StructPos = [];
            nn = nn + 1;
        end
    else
        UVPos = MultiVariateAnalysis.FindObjectInStructure(obj.Input{ii}{2},'UncertainVariable');
        for jj = 1:length(UVPos)
            UVList(nn).UV = getfield(obj.Input{ii}{2},UVPos{jj}{:});
            UVList(nn).Pos = ii;
            UVList(nn).StructPos = UVPos{jj};
            nn = nn+1;
        end
    end
end
if isempty(UVList); error('No uncertain variables defined in the input'); end

obj.UVInputList = UVList;
end
