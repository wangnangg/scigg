curr_dir:=${shell pwd}
cc_compiler:=g++
c_compiler:=gcc
gtest_dir:=${curr_dir}/googletest/googletest
src_dir:=${curr_dir}/src
test_dir:=${curr_dir}/test
example_dir:=${curr_dir}/example
build_base?=${curr_dir}/build
flags = -I${src_dir} -std=c++1z -Wall -Wfloat-conversion -Wsign-conversion -Werror -MMD
link_flags =-lstdc++ -lm -pthread -lblas
ifeq ($(release), 1)
	flags += -O3 -DNDEBUG
	link_flags += -O3 -DNDEBUG
	build_dir:=${build_base}/release
else
  ifeq ($(profile), 1)
	  flags += -g -O3 -pg -no-pie
	  link_flags += -g -O3 -no-pie
	  build_dir:=${build_base}/profile
  else
    ifeq ($(coverage), 1)
      flags += -g -ftest-coverage -fprofile-arcs
	  link_flags += -g -fprofile-arcs
	  build_dir:=${build_base}/coverage
    else
	  flags += -g
	  link_flags += -g
	  build_dir:=${build_base}/debug
    endif
  endif
endif
test_build_dir:=${build_dir}/test
utest_exe:=${test_build_dir}/utest

example_build_dir:=${build_dir}/example
example_exe:=${example_build_dir}/example

lib_file:=${build_dir}/libmarkovgg.a

cpp_files:=${wildcard ${src_dir}/*.cpp}
obj_files:=${addprefix ${build_dir}/,${notdir ${cpp_files:.cpp=.o}}}

test_cpp_files:=${wildcard ${test_dir}/*.cpp}
test_obj_files:=${addprefix ${test_build_dir}/,${notdir ${test_cpp_files:.cpp=.o}}}

example_cpp_files:=${wildcard ${example_dir}/*.cpp}
example_obj_files:=${addprefix ${example_build_dir}/,${notdir ${example_cpp_files:.cpp=.o}}}

deps=${obj_files:.o=.d} ${test_obj_files:.o=.d} ${example_obj_files:.o=.d}

gtest_obj:=${test_build_dir}/gtest-all.o

gtest_flags:=-isystem ${gtest_dir}/include -I${gtest_dir}


.PHONY: run_utest example lib
run_utest: ${utest_exe}
	${utest_exe} ${args}

example: ${example_exe}
	${example_exe} ${args}

lib: ${lib_file}

${lib_file}: ${obj_files}
	ar crf $@ $^

${utest_exe}: ${obj_files} ${test_obj_files} ${gtest_obj}
	${cc_compiler} $^ ${link_flags} -o $@

${example_exe}: ${obj_files} ${example_obj_files} ${gtest_obj}
	${cc_compiler} $^ ${link_flags} -o $@

.PHONY: build
build: ${utest_exe}

${build_dir}:
	mkdir -p $@

${test_build_dir}:
	mkdir -p $@

${example_build_dir}:
	mkdir -p $@

${lib_build_dir}:
	mkdir -p $@


.PHONY: clean
clean:
	rm -rf build

${test_build_dir}/%.o: ${test_dir}/%.cpp | ${test_build_dir}
	${cc_compiler} ${gtest_flags} ${flags} -c $< -o $@

${example_build_dir}/%.o: ${example_dir}/%.cpp | ${example_build_dir}
	${cc_compiler} ${gtest_flags} ${flags} -c $< -o $@

${build_dir}/%.o: src/%.cpp | ${build_dir}
	${cc_compiler} ${flags} -c $< -o $@

${gtest_obj}: ${gtest_dir}/src/gtest-all.cc | ${test_build_dir}
	${cc_compiler} ${gtest_flags} -c $^ -o $@

-include ${deps}
