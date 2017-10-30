function write_txt (fname, colname, rowname, data_matrix)

nhead_col = length(colname);
ndata_col = size(data_matrix, 2);
ndata_row = size(data_matrix, 1);

%-Open the file for writing
txtfid = fopen(fname, 'w');

%-Write first line (row names)
for headcolcnt = 1:nhead_col-1
  fprintf(txtfid, '%s\t', colname{headcolcnt});
end
fprintf(txtfid, '%s\n', colname{nhead_col});

%-Write the data matrix
for datarow_cnt = 1:ndata_row
  fprintf(txtfid, '%s\t', rowname{datarow_cnt});
  for datacol_cnt = 1:ndata_col-1
    fprintf(txtfid, '%.5f\t', data_matrix(datarow_cnt, datacol_cnt));
  end
  fprintf(txtfid, '%.5f\n', data_matrix(datarow_cnt, ndata_col));
end

%-Close the file
fclose(txtfid);

end