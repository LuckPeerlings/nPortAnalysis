function [ output_args ] = WriteToTextFile_PGF(DataFile,Headers, Data, Comment )

 fileID = fopen(DataFile,'w');
 
%Write the comment
for ii = 1:length(Comment)
    fprintf(fileID,['# ',Comment{ii},'\n']);
end
if length(Headers) ~= size(Data,2)
    error('The number of headers is not corresponding to the number of columns in the data file')
end
for ii = 1:length(Headers)
   fprintf(fileID,Headers{ii});
   if ii ~= length(Headers)
    fprintf(fileID,'\t');
   else
    fprintf(fileID,'\n');
   end
end
for ii = 1:size(Data,1)
    for jj = 1:size(Data,2)
        fprintf(fileID,'%e',Data(ii,jj));
        if jj ~= size(Data,2)
            fprintf(fileID,'\t');
        elseif ii ~= size(Data,1)
            fprintf(fileID,'\n');
        end
    end
end

fclose(fileID);


end

