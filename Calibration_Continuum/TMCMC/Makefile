CC = mpicc
CPP = mpic++
MPIF90 = mpif90


CFLAGS = -O3 -std=c99
CFLAGS += -Wall -D_XOPEN_SOURCE=700 -D_BSD_SOURCE


GSLROOT=/usr/local
CFLAGS  += -I${GSLROOT}/include


LDFLAGS += -L${GSLROOT}/lib -lgsl -lgslcblas

CFLAGS += `torc_cflags`
LDFLAGS += `torc_libs` -lpthread -lm

TARGETS = sample

all: $(TARGETS)



sample: engine_tmcmc.c engine_tmcmc.h spawner.h tmcmc_stats.c tmcmc_db.c auxil.c fitfun.h
	$(CC) $(CFLAGS) -o sample engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c spawner.c fitfun.c $(LDFLAGS)




clean:
	rm -f sample sample_theta sample_theta_fast sample_psi  *.o

clear:
	rm -rf *.txt
