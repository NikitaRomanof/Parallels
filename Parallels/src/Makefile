.PHONY: all clean check tests gcov_report dvi install uninstall dist
XX=g++

FILES_A= Ant/s21_graph_algorithms.cc \
	   Ant/console_interface.cc \
	   loader.cc \
	   Ant/mainant.cc \

FILES_W= Winograd/mainwinograd.cc \
	   Winograd/winograd.cc \
	   loader.cc \
	   Winograd/interfaceWinograd.cc \

FILES_G= Gauss/console_interface.cc \
		Gauss/gauss.cc \
		Gauss/main.cc

LDFLAGS= -Wall -Wextra -Werror -std=c++17 -g
OS := $(shell uname -s)
TEST_FLAGS= -lgtest -lgtest_main -lpthread
CFLAGS= -fprofile-arcs -ftest-coverage
BD= build/

ifeq ($(OS), Darwin)
	CC += -D OS_MAC
else
	CC += -D OS_LINUX
endif

all: clean
	$(XX) $(LDFLAGS) $(FILES_A) -o ant.out
	$(XX) $(LDFLAGS) $(FILES_W) -o winograd.out
	$(XX) $(LDFLAGS) $(FILES_G) -o gauss.out

ant: clean
	$(XX) $(LDFLAGS) $(FILES_A) -o ant.out

winograd: clean
	$(XX) $(LDFLAGS) $(FILES_W) -o winograd.out

gauss: clean
	$(XX) $(LDFLAGS) $(FILES_G) -o gauss.out

clean:
	rm -rf *.out *.a *.o *.cfg *.gcda *.gcno *.html *.info *.dSYM report

check: 
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	cp ../materials/linters/suppressions.cfg suppressions.cfg
	python3 ../materials/linters/cpplint.py *.cc *.h Ant/*.cc Ant/*.h Gauss/*.cc Gauss//*.h Winograd/*.cc  Winograd/*.h
	cppcheck -q --std=c++17 --enable=all --language=c++ --suppressions-list=suppressions.cfg *.cc *.h Ant/*.cc Ant/*.h Gauss/*.cc Gauss//*.h Winograd/*.cc  Winograd/*.h
