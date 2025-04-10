gcc -c sp_ienv.c
gcc -c -I../../../SuperLU/SRC/ dSLUsolve.c 
ar -cr libslu.a sp_ienv.o dSLUsolve.o
mv libslu.a /fermion2/greene/tonzani/HHG_CPC/make/../src/../lib
