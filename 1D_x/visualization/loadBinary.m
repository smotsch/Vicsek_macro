function l=loadBinary(nameFile)
  %% read the binary data in the file 'nameFile' in 'l'
  fid = fopen(nameFile, 'r');
  fseek(fid, 4, 'bof');
  l = fread(fid, 'double')';
  fclose(fid);
end
