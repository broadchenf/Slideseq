function BeadMappingFile=FindMostRecentMapping(PuckDirectory)
    %Given a puck directory, this just finds the directory of the most
    %recent puck mapping.
        directorymap=dir(PuckDirectory);

        for k=find([directorymap.isdir] & string({directorymap.name})~='.'& string({directorymap.name})~='..' )
            datearray=[];
            filenamesplit=split(directorymap(k).name,'_');
            filenamesplitsplit=split(filenamesplit(2),'-');
            try
                datearray(1)=str2num(char(filenamesplitsplit(1)));
                datearray(2)=str2num(char(filenamesplitsplit(2)));
                datearray(3)=str2num(char(filenamesplit(3)));
                RunTime(k)=2400*31*datearray(1)+2400*datearray(2)+datearray(3); %this will fail at the year boundary, though
            catch
                RunTime(k)=-Inf; %Files that don't end in properly formatted dates will be ignored
            end                
        end
        latestrun=find(RunTime==max(RunTime));
        BeadMappingFile=directorymap(latestrun).name;