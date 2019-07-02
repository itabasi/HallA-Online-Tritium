#! /bin/csh
foreach f (./rawdata/cyric*.kdf)
echo ${f}
echo ${f} | sed s@rawdata@data@ | sed s/\.kdf/\_02.root/
#./bin/offline ${f} -w `echo ${f} | sed s@rawdata@data@ | sed s/\.kdf/\_02.root/` -p param/cyric02.param
#./bin/offline rawdata/cyric030.kdf -w ./data/cyric030_01.root -p param/cyric02.param
end
