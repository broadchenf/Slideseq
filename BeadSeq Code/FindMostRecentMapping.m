function BeadMappingFile=FindMostRecentMapping(PuckDirectory)
    %Given a puck directory, this just finds the directory of the most
    %recent puck mapping.
        directorymap=dir(PuckDirectory);

        for k=find([directorymap.isdir] & string({directorymap.name})~='.'& string({directorymap.name})~='..' )
            datearray=[];
            filenamesplit=split(directorymap(k).name,'_');
            if filenamesplit{1}~="BeadMapping"
                continue
            end
            try
                yr=char(filenamesplit(2));
                if yr(1:2)=="20" %This cannot be true if filenamesplit(2) is a month
                    filenamesplitsplit=split(filenamesplit(3),'-');
                    datearray(1)=str2num(char(filenamesplitsplit(1)));
                    datearray(2)=str2num(char(filenamesplitsplit(2)));
                    datearray(3)=str2num(char(filenamesplit(4)));
                    datearray(4)=str2num(char(filenamesplit(2)));
                else
                    filenamesplitsplit=split(filenamesplit(2),'-');
                    datearray(1)=str2num(char(filenamesplitsplit(1)));
                    datearray(2)=str2num(char(filenamesplitsplit(2)));
                    datearray(3)=str2num(char(filenamesplit(3)));
                    datearray(4)=2018;
                end
                RunTime(k)=2400*31*12*datearray(4)+2400*31*datearray(1)+2400*datearray(2)+datearray(3);
            catch
                RunTime(k)=-Inf; %Files that don't end in properly formatted dates will be ignored
            end                
        end
        latestrun=find(RunTime==max(RunTime));
        BeadMappingFile=directorymap(latestrun).name;