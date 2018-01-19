
config?=debug
build_dir:=build/${config}
flags:= -Isrc -std=c++1z -Wall -Wfloat-conversion -Wsign-conversion -Werror -MMD
gtest_dir:= googletest/googletest
gtest_flags:=-isystem ${gtest_dir}/include -I${gtest_dir}
link_flags:=-lstdc++ -lm -pthread -lblas
ifeq ($(config), release)
  flags += -O3 -DNDEBUG
  link_flags += -O3 -DNDEBUG
else
  ifeq ($(config), profile)
    flags += -g -O3 -pg -no-pie
    link_flags += -g -O3 -no-pie
  else
    ifeq ($(config), coverage)
      flags += -g -ftest-coverage -fprofile-arcs
      link_flags += -g -fprofile-arcs
    else
      ifeq ($(config), debug)
        flags += -g
        link_flags += -g
      else
$(error Unknown config: $(config))
      endif
    endif
  endif
endif

.PHONY: utest lib
utest: ${build_dir}/utest
lib: ${build_dir}/libscigg.a
${build_dir}/src/graph/digraph_algo.o: src/graph/digraph_algo.cpp | ${build_dir}/src/graph 
	g++ ${flags} -c src/graph/digraph_algo.cpp -o ${build_dir}/src/graph/digraph_algo.o
${build_dir}/src/graph:
	mkdir -p $@
${build_dir}/src/matvec/vector.o: src/matvec/vector.cpp | ${build_dir}/src/matvec 
	g++ ${flags} -c src/matvec/vector.cpp -o ${build_dir}/src/matvec/vector.o
${build_dir}/src/matvec:
	mkdir -p $@
${build_dir}/src/matvec/matvec_oper.o: src/matvec/matvec_oper.cpp | ${build_dir}/src/matvec 
	g++ ${flags} -c src/matvec/matvec_oper.cpp -o ${build_dir}/src/matvec/matvec_oper.o
${build_dir}/src/matvec/blas.o: src/matvec/blas.cpp | ${build_dir}/src/matvec 
	g++ ${flags} -c src/matvec/blas.cpp -o ${build_dir}/src/matvec/blas.o
${build_dir}/src/matvec/matrix.o: src/matvec/matrix.cpp | ${build_dir}/src/matvec 
	g++ ${flags} -c src/matvec/matrix.cpp -o ${build_dir}/src/matvec/matrix.o
${build_dir}/src/spmatvec/spblas.o: src/spmatvec/spblas.cpp | ${build_dir}/src/spmatvec 
	g++ ${flags} -c src/spmatvec/spblas.cpp -o ${build_dir}/src/spmatvec/spblas.o
${build_dir}/src/spmatvec:
	mkdir -p $@
${build_dir}/src/spmatvec/spmatvec_oper.o: src/spmatvec/spmatvec_oper.cpp | ${build_dir}/src/spmatvec 
	g++ ${flags} -c src/spmatvec/spmatvec_oper.cpp -o ${build_dir}/src/spmatvec/spmatvec_oper.o
${build_dir}/src/spmatvec/spvector.o: src/spmatvec/spvector.cpp | ${build_dir}/src/spmatvec 
	g++ ${flags} -c src/spmatvec/spvector.cpp -o ${build_dir}/src/spmatvec/spvector.o
${build_dir}/src/spmatvec/spmatrix.o: src/spmatvec/spmatrix.cpp | ${build_dir}/src/spmatvec 
	g++ ${flags} -c src/spmatvec/spmatrix.cpp -o ${build_dir}/src/spmatvec/spmatrix.o
${build_dir}/src/debug/debug_utils.o: src/debug/debug_utils.cpp | ${build_dir}/src/debug 
	g++ ${flags} -c src/debug/debug_utils.cpp -o ${build_dir}/src/debug/debug_utils.o
${build_dir}/src/debug:
	mkdir -p $@
${build_dir}/src/splinalg/splinalg_eigen.o: src/splinalg/splinalg_eigen.cpp | ${build_dir}/src/splinalg 
	g++ ${flags} -c src/splinalg/splinalg_eigen.cpp -o ${build_dir}/src/splinalg/splinalg_eigen.o
${build_dir}/src/splinalg:
	mkdir -p $@
${build_dir}/src/splinalg/splinalg_decomp.o: src/splinalg/splinalg_decomp.cpp | ${build_dir}/src/splinalg 
	g++ ${flags} -c src/splinalg/splinalg_decomp.cpp -o ${build_dir}/src/splinalg/splinalg_decomp.o
${build_dir}/src/splinalg/splinalg_solve.o: src/splinalg/splinalg_solve.cpp | ${build_dir}/src/splinalg 
	g++ ${flags} -c src/splinalg/splinalg_solve.cpp -o ${build_dir}/src/splinalg/splinalg_solve.o
${build_dir}/src/optimize/optimize.o: src/optimize/optimize.cpp | ${build_dir}/src/optimize 
	g++ ${flags} -c src/optimize/optimize.cpp -o ${build_dir}/src/optimize/optimize.o
${build_dir}/src/optimize:
	mkdir -p $@
${build_dir}/src/linalg/linalg_least_square.o: src/linalg/linalg_least_square.cpp | ${build_dir}/src/linalg 
	g++ ${flags} -c src/linalg/linalg_least_square.cpp -o ${build_dir}/src/linalg/linalg_least_square.o
${build_dir}/src/linalg:
	mkdir -p $@
${build_dir}/src/linalg/linalg_decomp.o: src/linalg/linalg_decomp.cpp | ${build_dir}/src/linalg 
	g++ ${flags} -c src/linalg/linalg_decomp.cpp -o ${build_dir}/src/linalg/linalg_decomp.o
${build_dir}/src/linalg/linalg_solve.o: src/linalg/linalg_solve.cpp | ${build_dir}/src/linalg 
	g++ ${flags} -c src/linalg/linalg_solve.cpp -o ${build_dir}/src/linalg/linalg_solve.o
${build_dir}/test/test_digraph.o: test/test_digraph.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_digraph.cpp -o ${build_dir}/test/test_digraph.o
${build_dir}/test:
	mkdir -p $@
${build_dir}/test/test_spmatrix.o: test/test_spmatrix.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_spmatrix.cpp -o ${build_dir}/test/test_spmatrix.o
${build_dir}/test/test_optimize.o: test/test_optimize.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_optimize.cpp -o ${build_dir}/test/test_optimize.o
${build_dir}/test/test_main.o: test/test_main.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_main.cpp -o ${build_dir}/test/test_main.o
${build_dir}/test/test_linalg.o: test/test_linalg.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_linalg.cpp -o ${build_dir}/test/test_linalg.o
${build_dir}/test/test_spblas.o: test/test_spblas.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_spblas.cpp -o ${build_dir}/test/test_spblas.o
${build_dir}/test/test_matvec_oper.o: test/test_matvec_oper.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_matvec_oper.cpp -o ${build_dir}/test/test_matvec_oper.o
${build_dir}/test/test_common.o: test/test_common.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_common.cpp -o ${build_dir}/test/test_common.o
${build_dir}/test/test_blas.o: test/test_blas.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_blas.cpp -o ${build_dir}/test/test_blas.o
${build_dir}/test/test_splinalg.o: test/test_splinalg.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_splinalg.cpp -o ${build_dir}/test/test_splinalg.o
${build_dir}/test/test_large.o: test/test_large.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_large.cpp -o ${build_dir}/test/test_large.o
${build_dir}/test/common.o: test/common.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/common.cpp -o ${build_dir}/test/common.o
${build_dir}/test/test_matrix.o: test/test_matrix.cpp | ${build_dir}/test 
	g++ ${gtest_flags} ${flags} -c test/test_matrix.cpp -o ${build_dir}/test/test_matrix.o
${build_dir}/${gtest_dir}/src/gtest-all.o: ${gtest_dir}/src/gtest-all.cc | ${build_dir}/${gtest_dir}/src 
	g++ ${gtest_flags} -c ${gtest_dir}/src/gtest-all.cc -o ${build_dir}/${gtest_dir}/src/gtest-all.o
${build_dir}/${gtest_dir}/src:
	mkdir -p $@
${build_dir}/utest: ${build_dir}/src/graph/digraph_algo.o ${build_dir}/src/matvec/vector.o ${build_dir}/src/matvec/matvec_oper.o ${build_dir}/src/matvec/blas.o ${build_dir}/src/matvec/matrix.o ${build_dir}/src/spmatvec/spblas.o ${build_dir}/src/spmatvec/spmatvec_oper.o ${build_dir}/src/spmatvec/spvector.o ${build_dir}/src/spmatvec/spmatrix.o ${build_dir}/src/debug/debug_utils.o ${build_dir}/src/splinalg/splinalg_eigen.o ${build_dir}/src/splinalg/splinalg_decomp.o ${build_dir}/src/splinalg/splinalg_solve.o ${build_dir}/src/optimize/optimize.o ${build_dir}/src/linalg/linalg_least_square.o ${build_dir}/src/linalg/linalg_decomp.o ${build_dir}/src/linalg/linalg_solve.o ${build_dir}/test/test_digraph.o ${build_dir}/test/test_spmatrix.o ${build_dir}/test/test_optimize.o ${build_dir}/test/test_main.o ${build_dir}/test/test_linalg.o ${build_dir}/test/test_spblas.o ${build_dir}/test/test_matvec_oper.o ${build_dir}/test/test_common.o ${build_dir}/test/test_blas.o ${build_dir}/test/test_splinalg.o ${build_dir}/test/test_large.o ${build_dir}/test/common.o ${build_dir}/test/test_matrix.o ${build_dir}/${gtest_dir}/src/gtest-all.o  | ${build_dir}
	g++ ${build_dir}/src/graph/digraph_algo.o ${build_dir}/src/matvec/vector.o ${build_dir}/src/matvec/matvec_oper.o ${build_dir}/src/matvec/blas.o ${build_dir}/src/matvec/matrix.o ${build_dir}/src/spmatvec/spblas.o ${build_dir}/src/spmatvec/spmatvec_oper.o ${build_dir}/src/spmatvec/spvector.o ${build_dir}/src/spmatvec/spmatrix.o ${build_dir}/src/debug/debug_utils.o ${build_dir}/src/splinalg/splinalg_eigen.o ${build_dir}/src/splinalg/splinalg_decomp.o ${build_dir}/src/splinalg/splinalg_solve.o ${build_dir}/src/optimize/optimize.o ${build_dir}/src/linalg/linalg_least_square.o ${build_dir}/src/linalg/linalg_decomp.o ${build_dir}/src/linalg/linalg_solve.o ${build_dir}/test/test_digraph.o ${build_dir}/test/test_spmatrix.o ${build_dir}/test/test_optimize.o ${build_dir}/test/test_main.o ${build_dir}/test/test_linalg.o ${build_dir}/test/test_spblas.o ${build_dir}/test/test_matvec_oper.o ${build_dir}/test/test_common.o ${build_dir}/test/test_blas.o ${build_dir}/test/test_splinalg.o ${build_dir}/test/test_large.o ${build_dir}/test/common.o ${build_dir}/test/test_matrix.o ${build_dir}/${gtest_dir}/src/gtest-all.o  ${link_flags} -o ${build_dir}/utest
${build_dir}:
	mkdir -p $@
${build_dir}/libscigg.a: ${build_dir}/src/graph/digraph_algo.o ${build_dir}/src/matvec/vector.o ${build_dir}/src/matvec/matvec_oper.o ${build_dir}/src/matvec/blas.o ${build_dir}/src/matvec/matrix.o ${build_dir}/src/spmatvec/spblas.o ${build_dir}/src/spmatvec/spmatvec_oper.o ${build_dir}/src/spmatvec/spvector.o ${build_dir}/src/spmatvec/spmatrix.o ${build_dir}/src/debug/debug_utils.o ${build_dir}/src/splinalg/splinalg_eigen.o ${build_dir}/src/splinalg/splinalg_decomp.o ${build_dir}/src/splinalg/splinalg_solve.o ${build_dir}/src/optimize/optimize.o ${build_dir}/src/linalg/linalg_least_square.o ${build_dir}/src/linalg/linalg_decomp.o ${build_dir}/src/linalg/linalg_solve.o  | ${build_dir}
	ar crf ${build_dir}/libscigg.a ${build_dir}/src/graph/digraph_algo.o ${build_dir}/src/matvec/vector.o ${build_dir}/src/matvec/matvec_oper.o ${build_dir}/src/matvec/blas.o ${build_dir}/src/matvec/matrix.o ${build_dir}/src/spmatvec/spblas.o ${build_dir}/src/spmatvec/spmatvec_oper.o ${build_dir}/src/spmatvec/spvector.o ${build_dir}/src/spmatvec/spmatrix.o ${build_dir}/src/debug/debug_utils.o ${build_dir}/src/splinalg/splinalg_eigen.o ${build_dir}/src/splinalg/splinalg_decomp.o ${build_dir}/src/splinalg/splinalg_solve.o ${build_dir}/src/optimize/optimize.o ${build_dir}/src/linalg/linalg_least_square.o ${build_dir}/src/linalg/linalg_decomp.o ${build_dir}/src/linalg/linalg_solve.o 
deps:=${build_dir}/src/graph/digraph_algo.d ${build_dir}/src/matvec/vector.d ${build_dir}/src/matvec/matvec_oper.d ${build_dir}/src/matvec/blas.d ${build_dir}/src/matvec/matrix.d ${build_dir}/src/spmatvec/spblas.d ${build_dir}/src/spmatvec/spmatvec_oper.d ${build_dir}/src/spmatvec/spvector.d ${build_dir}/src/spmatvec/spmatrix.d ${build_dir}/src/debug/debug_utils.d ${build_dir}/src/splinalg/splinalg_eigen.d ${build_dir}/src/splinalg/splinalg_decomp.d ${build_dir}/src/splinalg/splinalg_solve.d ${build_dir}/src/optimize/optimize.d ${build_dir}/src/linalg/linalg_least_square.d ${build_dir}/src/linalg/linalg_decomp.d ${build_dir}/src/linalg/linalg_solve.d ${build_dir}/test/test_digraph.d ${build_dir}/test/test_spmatrix.d ${build_dir}/test/test_optimize.d ${build_dir}/test/test_main.d ${build_dir}/test/test_linalg.d ${build_dir}/test/test_spblas.d ${build_dir}/test/test_matvec_oper.d ${build_dir}/test/test_common.d ${build_dir}/test/test_blas.d ${build_dir}/test/test_splinalg.d ${build_dir}/test/test_large.d ${build_dir}/test/common.d ${build_dir}/test/test_matrix.d ${build_dir}/${gtest_dir}/src/gtest-all.d 

.PHONY: clean
clean:
	  rm build -rf
-include ${deps}
