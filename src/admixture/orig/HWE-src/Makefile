#***************************************************************************
# Makefile
#***************************************************************************

NAME   = hwe

CFLAGS = -O
CC     = gcc
LIB    = -lm

OBJS = $(NAME).o cal_const.o cal_n.o cal_prob.o\
		 check_file.o do_switch.o\
		 ln_p_value.o log_factorial.o print_data.o\
		random_choose.o  read_data.o select_index.o\
		 stamp_time.o test_switch.o 

$(OBJS): hwe.h func.h
$(NAME): $(OBJS)
	$(CC) $(CFLAGS) -o $(NAME) $(OBJS) $(LIB)

lixc:
	lint $(NAME).c cal_const.c cal_n.c cal_prob.c\
	check_file.c do_switch.c\ 
	ln_p_value.c log_factorial.c print_data.c\
	random_choose.c  read_data.c select_index.c stamp_time.c\
	test_switch.c 

print: *.h *.o
	pr $? | lpr -Pstat
		 touch print

clean:
	rm -f *.o

