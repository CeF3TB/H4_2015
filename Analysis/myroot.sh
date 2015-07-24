#!/bin/bash

### SOURCE ROOT
echo "Sourcing ROOT v5.34.18 py26"
source /swshare/ROOT/root_v5.34.18_slc5_amd64_py26_pythia6/bin/thisroot.sh

### LINK correct mysql
#g++ --shared -fPIC -Wl,-as-needed -o libRMySQL.so /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/TMySQLResult.o  /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/TMySQLRow.o  /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/TMySQLServer.o  /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/TMySQLStatement.o /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/G__MySQL.o  -Llib -rdynamic -L/usr/lib/mysql -lmysqlclient -lz -lcrypt -lnsl -lm -L/usr/lib -L/usr/lib64 -lssl -lcrypto -Llib -lCore -lCint -ldl -L/usr/lib64/mysql
echo "Making dir $PWD/lib"
[ -d $PWD/lib ] || mkdir -p $PWD/lib
[ -f lib/libRMySQL.so ] || g++ --shared -fPIC -Wl,-as-needed -o lib/libRMySQL.so /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/TMySQLResult.o  /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/TMySQLRow.o  /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/TMySQLServer.o  /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/TMySQLStatement.o /swshare/ROOT/root_v5.34.18.source/sql/mysql/src/G__MySQL.o  -Llib -rdynamic -L/usr/lib/mysql -lmysqlclient -lz -lcrypt -lnsl -lm -L/usr/lib -L/usr/lib64 -lssl -lcrypto -Llib -lCore -lCint -ldl -L/usr/lib64/mysql `root-config --libs`

LD_LIBRARY_PATH="$PWD/lib":$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH is $LD_LIBRARY_PATH"
