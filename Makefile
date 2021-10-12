SHELL=/bin/ksh
SRCS=	gsw.f prepare.f cons.f bogusuv.f \
	gausl3.f getrdy.f outflds.f \
	lgndr.f sortml.f tranrs.f \
	transr.f fftx.f make_list.f \
	joinrs.f joinsr.f ujoinrs.f ujoinsr.f mpe_transpose_rs.f \
	mpe_transpose_sr.f mpe_unify.f mpe_init.f mpe_finalize.f \
        intgrt.f trngra.f tracking.f hdiffu.f \
        tranuv.f tranrs1.f mpe_transpose_rs1.f 
 
OBJS=	gsw.o prepare.o cons.o bogusuv.o \
	gausl3.o getrdy.o outflds.o \
        lgndr.o sortml.o tranrs.o \
        transr.o fftx.o make_list.o \
        joinrs.o joinsr.o ujoinrs.o ujoinsr.o mpe_transpose_rs.o \
	mpe_transpose_sr.o mpe_unify.o mpe_init.o mpe_finalize.o \
        intgrt.o trngra.o tracking.o hdiffu.o \
        tranuv.o tranrs1.o mpe_transpose_rs1.o 

#FC =     /package/dms/dms.v4/dmsmpif77
FC =    /users/xa09/dms_compile/fx10/dmsmpif77
FC90 =          mpifrtpx

CMD =           vor720.exe

#FFLAGS =  -O3 -q64 -qsave -qrealsize=8 -qarch=auto -qfixed=132
FFLAGS = -c -Kfast -CcdRR8 -Fwide
#LDFLAGS = -lmass -lmassvp4 -lessl -bmaxdata:0x12000000000 -bmaxstack:0x10000000000
LDFLAGS = -L/package/fx10/dms/dms.v4/lib -lrdms -lgdbm -L/users/xa09/operlib/lib -lnwp -SSL2

all:            $(CMD)

$(CMD):         $(OBJS)
#cc -c +DD64 +DSitanium2 timelib.c
#	xlf -c -q64 writes.f
#	xlf -c -q64 reads.f
	frtpx -c writes.f reads.f
	$(FC) -o $(@) $(LDFLAGS)  $(OBJS) writes.o reads.o

clean:
	-rm -f $(OBJS) *.lst *.exe
 
clobber:        clean
	-rm -f $(CMD)
 
void:   clobber
	-rm -f $(SRCS)
