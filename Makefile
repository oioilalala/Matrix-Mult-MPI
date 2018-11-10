targets := mmm_mpi mandelbrot_mpi 
mobjs	:= matrix_checksum.o mmm_mpi.o
objs	:= mandelbrot_mpi.o 
CC	:= mpicc

CFLAGS	:= -Wall -Werror -O2
CFLAGS	+= -g

ifneq ($(V),1)
Q = @
endif

all: $(targets)

mmm_mpi: $(mobjs)
	@echo "CC $@"
	$(Q)$(CC) $(CFLAGS) -o $@ $^ -lm

mandelbrot_mpi: $(objs)
	@echo "CC $@"
	$(Q)$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c
	@echo "CC $@"
	$(Q)$(CC) $(CFLAGS) -c -o $@ $<

%mpi.o: %mpi.c
	@echo "CC $@"
	$(Q)$(CC) $(CFLAGS) -c -o $@ $<

clean:
	@echo "clean"
	$(Q)rm -f $(targets) $(mobjs) $(objs) 
