curr_dir:=${shell pwd}
compiler:=g++
gtest_dir:=${curr_dir}/googletest/googletest
src_dir:=${curr_dir}/src
test_dir:=${curr_dir}/test
example_dir:=${curr_dir}/example
ifeq ($(release), 1)
	flags:=-I${src_dir} -std=c++1z -Werror -MMD -O3 -DNDEBUG
	link_flags:=-lstdc++ -lm -pthread -O3 -DNDEBUG
	build_dir:=${curr_dir}/build/release
else
  ifeq ($(profile), 1)
	  flags:=-I${src_dir} -std=c++1z -Werror -MMD -g -O3 -pg -no-pie
	  link_flags:=-lstdc++ -lm -pthread -g -O3 -no-pie
	  build_dir:=${curr_dir}/build/profile
  else
	  flags:=-I${src_dir} -std=c++1z -Werror -MMD -g
	  link_flags:=-lstdc++ -lm -pthread -g
	  build_dir:=${curr_dir}/build/debug
  endif
endif
test_build_dir:=${build_dir}/test
utest_exe:=${test_build_dir}/utest

example_build_dir:=${build_dir}/example
example_exe:=${example_build_dir}/example

cpp_files:=${wildcard ${src_dir}/*.cpp}
obj_files:=${addprefix ${build_dir}/,${notdir ${cpp_files:.cpp=.o}}}
test_cpp_files:=${wildcard ${test_dir}/*.cpp}
test_obj_files:=${addprefix ${test_build_dir}/,${notdir ${test_cpp_files:.cpp=.o}}}
example_cpp_files:=${wildcard ${example_dir}/*.cpp}
example_obj_files:=${addprefix ${example_build_dir}/,${notdir ${example_cpp_files:.cpp=.o}}}

deps=${obj_files:.o=.d} ${test_obj_files:.o=.d} ${example_obj_files:.o=.d}

gtest_obj:=${test_build_dir}/gtest-all.o

gtest_flags:=-isystem ${gtest_dir}/include -I${gtest_dir}


.PHONY: run_utest example
run_utest: ${utest_exe}
	${utest_exe} ${args}

example: ${example_exe}
	${example_exe} ${args}

${utest_exe}: ${obj_files} ${test_obj_files} ${gtest_obj}
	${compiler} $^ ${link_flags} -o $@

${example_exe}: ${obj_files} ${example_obj_files} ${gtest_obj}
	${compiler} $^ ${link_flags} -o $@

${gtest_obj}: ${gtest_dir}/src/gtest-all.cc | ${test_build_dir}
	${compiler} ${gtest_flags} -c $^ -o $@

.PHONY: build
build: ${utest_exe}

${build_dir}:
	mkdir -p $@

${test_build_dir}:
	mkdir -p $@

${example_build_dir}:
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf build

${test_build_dir}/%.o: ${test_dir}/%.cpp | ${test_build_dir}
	${compiler} ${gtest_flags} ${flags} -c $< -o $@

${example_build_dir}/%.o: ${example_dir}/%.cpp | ${example_build_dir}
	${compiler} ${gtest_flags} ${flags} -c $< -o $@

${build_dir}/%.o: src/%.cpp | ${build_dir}
	${compiler} ${flags} -c $< -o $@

-include ${deps}