function Output = UpdateSolution(obj,Input)
%Setting the input parameter for the 
for ii = 1:length(Input)
    set(obj.ClassHandle,Input{ii}{1},Input{ii}{2});
end

%Execute the MethodHandle(s) to update the solution
for ii = 1:length(obj.MethodHandles)
    obj.ClassHandle.(obj.MethodHandles{ii});
end
%And obtain the base output
for ii = 1:length(obj.OutputProperties)
    Output{ii} = get(obj.ClassHandle,obj.OutputProperties{ii});
end
end