function list = rd_list(file)

%RD_LIST Reads a flat list file.
%       LIST = RD_LIST(FILE)
%  LIST         - Cell array where every line in FILE is an element of the array
%  FILE - Ascii file, usually a listing of files but could be anything
%  Lines can be commented out by using a % at the start of the line. These lines will
%       not show up in LIST.
%
% Added extra line to initially assign list to an empty matrix, in case FILE is empty. John Newgard, June 2010

list = {};
fid = fopen(file);

ct = 0; notdone = 1;
while notdone
   lin = fgetl(fid);
   if lin == -1; break; end
   if ~strcmp('%',lin(1))
      ct = ct +1;
      list{ct} = lin;
   end   
end

fclose(fid);
