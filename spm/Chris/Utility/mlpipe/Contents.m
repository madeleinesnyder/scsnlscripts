% 
% mlpipe: A pipe I/O library for Matlab and Octave.
%     
%     (c) Copyright Jonas Larsson 2008
% 
%     This file is part of mlpipe: a pipe I/O library for Matlab and Octave.
% 
%     mlpipe is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     mlpipe is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with mlpipe.  If not, see <http://www.gnu.org/licenses/>.
% 
%     PIPE I/O COMMANDS
%     ===============================================================
%     General-purpose pipe manipulation commands. Compressed files are 
%     automatically treated as pipes and opened with the corresponding
%     file compression/decompression utility [gzip/gunzip (.gz), 
%     compress/uncompress (.Z), zip/unzip (.zip)] which must be installed
%     on the system (default on most modern Unix systems).
%
%     NB. THESE FUNCTIONS RELY ON UNIX PIPES AND ARE NOT SUPPORTED ON WINDOWS 
%     (but may work with Cygwin)
%
%     mlpopen    :  Opens a pipe, returns information about open pipe/s
%     mlpclose   :  Closes a pipe or all pipes
%     
%     mlpseek    :  Move position indicator of a read pipe
%     mlprewind  :  Resets a read pipe to the beginning-of-pipe
%     mlptell    :  Returns position indicator of an open pipe
% 
%     mlpread    :  Reads binary data from a pipe
%     mlpwrite   :  Writes binary data to a pipe
%
%     COMPRESSED FILE I/O COMMANDS
%     ================================================================
%     Convenience functions providing a unified I/O interface to read 
%     from and write to normal and compressed files. Call the appropriate
%     pread/pwrite etc or fread/fwrite etc functions, depending on file type.
%
%     zopen    :  Opens a pipe or file, returns information about open pipe/s & files
%     zclose   :  Closes one or more files and/or pipes
%     
%     zseek    :  Move position indicator of a read pipe or file
%     zrewind  :  Resets a read pipe to the beginning-of-pipe or -file
%     ztell    :  Returns position indicator of an open pipe or file
% 
%     zread    :  Reads binary data from a pipe or file
%     zwrite   :  Writes binary data to a pipe or file
%
% Note for Octave users:
%  Octave natively provides support for pipes (unlike Matlab) with the commands
%  popen and pclose, which returns file identifiers that can be accessed
%  using standard file I/O commands (fread etc).
%  The functions in this packages are mainly useful on Octave if you 
%  share code across the two environments; otherwise use the built-in 
%  Octave functions. In particular, note that the functions mlpopen and 
%  mlpclose are NOT compatible with the Octave functions popen and pclose.

