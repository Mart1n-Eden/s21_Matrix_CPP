CC = g++
EFLAGS = -Wall -Werror -Wextra -std=c++17 
CFLAGS = $(EFLAGS) -c

TEST_FLAGS = -lgtest -lstdc++
GCOVFLAGS = $(TEST_FLAGS) -L. --coverage

OS = $(shell uname)
ifeq ($(OS), Linux)
 TEST_FLAGS += $(shell pkg-config --libs gtest)
endif


all: s21_matrix_oop.a

s21_matrix_oop.a:
	$(CC) $(CFLAGS) -lstdc++ s21_matrix_oop.cc -o s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

test: clean s21_matrix_oop.a 
	$(CC) $(EFLAGS) $(TEST_FLAGS) s21_matrix_oop.a test.cc -o test
	./test

gcov_report: clean
	$(CC) $(GCOVFLAGS) test.cc s21_matrix_oop.cc -o gcov_report
	./gcov_report
	lcov -t "stest" -o s21_test.info -c -d . --ignore-errors mismatch 
	genhtml -o report s21_test.info
	open ./report/index.html

clean:
	rm -rf *.o *.out s21_matrix_oop.a test
	rm -rf *.gcda *.gcno report gcov_report.* gcov_report *.info gcovreport

leaks: test
	leaks -atExit -- ./test

clang:
	cp ../materials/linters/.clang-format ../src/.clang-format
	clang-format -style=Google -n *.cc *.h
	rm .clang-format