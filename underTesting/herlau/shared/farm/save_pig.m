function save_pig(pig)
%name = pig.name;
%parfile = pig.parfile;
%runs = pig.runs;
%a = pig.a; 
%com = pig.com;
%submitted = pig.submitted;


save(pig.name,'-struct','pig');
%save(pig.name, 'name', 'parfile', 'runs', 'a', 'com', 'submitted');
end