function fid = writeFastScapeVelocity(outfile,x,t,Vx,Vy)

fid = fopen(outfile,'w');

nx = length(x);
nt = length(t);

fwrite(fid,nx,'int32');
fwrite(fid,nt,'int32');
fwrite(fid,x,'double');
fwrite(fid,t,'double');
fwrite(fid,Vx','double');
fwrite(fid,Vy','double');

fid = fclose(fid);

