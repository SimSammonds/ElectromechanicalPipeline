function [gSWdata,iSWdata] = organSWs(SWdata,g_list)
% organSWs is a function to segregate organised SWdata by the organ of
% origin

clsts = unique(SWdata(:,1));
gSWdata = [];
iSWdata = [];

for clst = clsts'
   idx = SWdata(:,1) == clst;
   if ismember(clst,g_list)
       gSWdata = [gSWdata; SWdata(idx,:)];
   else
       iSWdata = [iSWdata; SWdata(idx,:)];
   end
end

end