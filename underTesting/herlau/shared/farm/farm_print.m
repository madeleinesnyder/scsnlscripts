function farm_print(pigfile, message)
if nargin < 2, message = 'default'; end
pig = load(pigfile);
pig.com.t = tic();
pig.com.message = message;
save_pig(pig);
end