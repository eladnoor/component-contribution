function fullpath = getBasePath()

fullpath = which('getBasePath.m');
fullpath = regexprep(fullpath,'getBasePath.m','');
fullpath = [fullpath '../'];
